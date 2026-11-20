"""Download every file in data/processed/manifest.parquet, with resume and MD5 checks.

- Streams each file to <name>.part, resuming with an HTTP Range request when possible,
  checks the MD5 against the manifest, then renames it into place.
- Files already present with the right MD5 are skipped, so the script can be re-run safely.
- After downloading: unzips the hg19 FASTA and builds its .fai index (pyfaidx), and
  unzips the liftover chain (rtracklayer::import.chain needs it uncompressed).
- Writes data/processed/download_log.tsv; exits non-zero if any file failed.

Usage:
  .venv/Scripts/python -I preprocess/python/02_download.py --dry-run
  .venv/Scripts/python -I preprocess/python/02_download.py [--category kd_counts ...] [--workers 4]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "data" / "processed" / "manifest.parquet"
LOG = ROOT / "data" / "processed" / "download_log.tsv"
CHUNK = 1 << 20
RETRIES = 4
PAUSE = 0.5  # seconds between requests per worker, to stay polite to the ENCODE server
BLOCK_LIMIT = 10  # consecutive 403/429 responses before the run stops (server is refusing us)
_local = threading.local()
_block = {"consecutive": 0, "stop": threading.Event(), "lock": threading.Lock()}


class Blocked(Exception):
    """Server refused the request (403/429), most likely rate limiting."""


def _note_status(code: int) -> None:
    with _block["lock"]:
        if code in (403, 429):
            _block["consecutive"] += 1
            if _block["consecutive"] >= BLOCK_LIMIT:
                _block["stop"].set()
        else:
            _block["consecutive"] = 0


def session() -> requests.Session:
    if not hasattr(_local, "s"):
        _local.s = requests.Session()
        _local.s.headers["User-Agent"] = "encode-eclip-atlas downloader (python-requests)"
    return _local.s


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def fetch(rec: dict, bar: tqdm) -> dict:
    dest = ROOT / rec["local_path"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    out = {"accession": rec["accession"], "category": rec["category"], "local_path": rec["local_path"],
           "status": "", "md5_ok": False, "bytes": 0, "error": ""}

    if dest.exists() and dest.stat().st_size == rec["file_size"] and md5_of(dest) == rec["md5sum"]:
        bar.update(rec["file_size"])
        return {**out, "status": "skipped (verified)", "md5_ok": True, "bytes": dest.stat().st_size}

    part = dest.with_name(dest.name + ".part")
    for attempt in range(RETRIES):
        if _block["stop"].is_set():
            return {**out, "status": "not attempted (server blocking)", "error": "stopped after repeated 403/429"}
        try:
            time.sleep(PAUSE)
            have = part.stat().st_size if part.exists() else 0
            if rec["file_size"] and have > rec["file_size"]:
                part.unlink(); have = 0
            headers = {"Range": f"bytes={have}-"} if have else {}
            with session().get(rec["href"], headers=headers, stream=True, timeout=(30, 60)) as r:  # 60 s read timeout: retry stalled connections quickly
                _note_status(r.status_code)
                if r.status_code in (403, 429):
                    raise Blocked(f"HTTP {r.status_code}")
                if r.status_code == 416:  # .part already holds the whole file
                    bar.update(have)
                else:
                    r.raise_for_status()
                    resumed = bool(have) and r.status_code == 206  # else server ignored Range: restart
                    if resumed:
                        bar.update(have)
                    with open(part, "ab" if resumed else "wb") as fh:
                        for block in r.iter_content(CHUNK):
                            fh.write(block); bar.update(len(block))
            got = md5_of(part)
            if got == rec["md5sum"]:
                part.replace(dest)
                return {**out, "status": "downloaded", "md5_ok": True, "bytes": dest.stat().st_size}
            bar.update(-part.stat().st_size)
            part.unlink()
            out["error"] = f"md5 mismatch (got {got})"
        except Blocked as e:
            out["error"] = str(e)
            time.sleep(30 * (attempt + 1))  # long back-off when refused
            continue
        except (requests.RequestException, OSError) as e:
            out["error"] = f"{type(e).__name__}: {e}"
        time.sleep(2 ** attempt)
    if _block["stop"].is_set():
        return {**out, "status": "not attempted (server blocking)"}
    return {**out, "status": "failed"}


def gunzip(src: Path, dest: Path) -> bool:
    """Decompress src to dest unless dest is already newer; returns True if work was done."""
    if dest.exists() and dest.stat().st_mtime >= src.stat().st_mtime and dest.stat().st_size > 0:
        return False
    tmp = dest.with_name(dest.name + ".tmp")
    with gzip.open(src, "rb") as fi, open(tmp, "wb") as fo:
        shutil.copyfileobj(fi, fo, CHUNK)
    tmp.replace(dest)
    return True


def postprocess(man: pd.DataFrame) -> None:
    for cat, out_name in [("genome", "hg19.fa"), ("chain", "hg38ToHg19.over.chain")]:
        rows = man[man.category == cat]
        if rows.empty:
            continue
        src = ROOT / rows.local_path.iloc[0]
        if not src.exists():
            continue
        dest = src.parent / out_name
        print(f"Unzipping {src.name} -> {dest.name} ...", flush=True)
        changed = gunzip(src, dest)
        if cat == "genome":
            fai = dest.with_name(dest.name + ".fai")
            if changed and fai.exists():
                fai.unlink()
            if not fai.exists():
                print("Indexing hg19.fa with pyfaidx ...", flush=True)
                from pyfaidx import Fasta
                Fasta(str(dest)).close()
            with open(fai) as fh:
                n = sum(1 for _ in fh)
            print(f"  {fai.name}: {n} sequences")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print what would be downloaded and exit")
    ap.add_argument("--category", nargs="+", help="only these manifest categories")
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()

    man = pd.read_parquet(MANIFEST)
    if args.category:
        unknown = set(args.category) - set(man.category)
        if unknown:
            raise SystemExit(f"Unknown categories {sorted(unknown)}; choose from {sorted(man.category.unique())}")
        man = man[man.category.isin(args.category)]

    summary = man.groupby("category").agg(files=("accession", "size"), gb=("file_size", lambda s: s.sum() / 1e9))
    print(summary.to_string(float_format=lambda v: f"{v:.3f}"))
    print(f"TOTAL  {len(man)} files  {man.file_size.sum() / 1e9:.2f} GB")
    free = shutil.disk_usage(ROOT).free / 1e9
    print(f"Free disk: {free:.0f} GB")
    if args.dry_run:
        return 0
    if free < man.file_size.sum() / 1e9 * 2 + 5:
        raise SystemExit("Not enough free disk space.")

    recs = man.to_dict("records")
    results = []
    with tqdm(total=int(man.file_size.sum()), unit="B", unit_scale=True, unit_divisor=1024,
              desc="download", mininterval=2) as bar, ThreadPoolExecutor(args.workers) as pool:
        # Large files first so they overlap with the many small ones.
        futs = [pool.submit(fetch, r, bar) for r in sorted(recs, key=lambda r: -r["file_size"])]
        for f in as_completed(futs):
            res = f.result()
            results.append(res)
            if res["status"] == "failed":
                tqdm.write(f"FAILED {res['accession']}: {res['error']}")

    log = pd.DataFrame(results).sort_values(["category", "accession"])
    if LOG.exists() and args.category:  # keep log rows for categories not run this time
        old = pd.read_csv(LOG, sep="\t")
        log = pd.concat([old[~old.category.isin(args.category)], log])
    log.to_csv(LOG, sep="\t", index=False)
    print("\n" + log.groupby(["category", "status"]).size().to_string())

    if _block["stop"].is_set():
        print(f"\nStopped: the server refused {BLOCK_LIMIT} requests in a row (HTTP 403/429), "
              "most likely temporary rate limiting. Wait 30-60 minutes and re-run; "
              "verified files are skipped.", file=sys.stderr)
        return 2
    failed = (log.status == "failed").sum()
    if failed:
        print(f"\n{failed} file(s) failed; re-run to retry (completed files are skipped).", file=sys.stderr)
        return 1
    postprocess(man)
    print(f"\nAll {len(log)} files verified. Log: {LOG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
