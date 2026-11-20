"""Consistency checks for the shared pre-processing outputs (steps 01-04, 05a).

Prints PASS/FAIL per check and exits non-zero if any check fails. Read-only.

Usage:  .venv/Scripts/python -I preprocess/checks/check_preprocess.py [--quick]
        --quick  skip re-hashing every downloaded file (MD5, ~5.6 GB, ~1-2 min)
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data/processed"
TAB = ROOT / "preprocess/results/tables"
CANONICAL = [f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"]]
REGIONS = ["CDS", "3UTR", "5UTR", "small_ncRNA", "noncoding_exon", "proximal_intron", "distal_intron", "intergenic"]
NP_COLS = ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
results: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""), flush=True)


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    # ---- manifest ---------------------------------------------------------------------------
    spec = importlib.util.spec_from_file_location("bm", ROOT / "preprocess/python/01_build_manifest.py")
    bm = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "preprocess/python"))
    spec.loader.exec_module(bm)
    man = pd.read_parquet(PROC / "manifest.parquet")
    counts = man.category.value_counts()
    off = {c: (int(counts.get(c, 0)), n) for c, n in bm.EXPECTED.items() if int(counts.get(c, 0)) != n}
    check("manifest category counts equal the 2026-11-16 survey", not off, f"differences: {off}" if off else f"{len(man)} files")
    check("manifest: unique local paths, every file has an MD5",
          not man.local_path.duplicated().any() and man.md5sum.notna().all())
    missing = [p for p in man.local_path if not (ROOT / p).exists()]
    check("every manifest file exists locally", not missing, f"{len(missing)} missing" if missing else "")
    if not args.quick:
        bad = [r.local_path for r in man.itertuples() if (ROOT / r.local_path).exists()
               and md5_of(ROOT / r.local_path) != r.md5sum]
        check("every downloaded file matches its MD5", not bad, f"{len(bad)} mismatches: {bad[:3]}" if bad else "")
    for f in ["hg19.fa", "hg19.fa.fai", "hg38ToHg19.over.chain", "gencode_v19.sqlite"]:
        check(f"reference file present: data/ref/{f}", (ROOT / "data/ref" / f).exists())

    # ---- liftover ---------------------------------------------------------------------------
    lift = pd.read_csv(TAB / "liftover_stats.tsv", sep="\t")
    check("liftover: all GRCh38-only experiments lifted", len(lift) == int(man.needs_liftover.sum()), f"{len(lift)} sets")
    check("liftover: no experiment loses > 5% of peaks", (lift.pct_dropped <= 5).all(), f"max {lift.pct_dropped.max():.2f}%")
    check("liftover: output files exist", all((ROOT / p).exists() for p in lift.output))

    # ---- annotated peaks vs raw BED files ---------------------------------------------------
    peaks = pd.read_parquet(PROC / "peaks_annotated.parquet")
    ecl = man[man.category.str.startswith("eclip_peaks")].copy()
    ecl["peak_path"] = np.where(ecl.needs_liftover,
                                "data/raw/eclip_peaks_hg19/" + ecl.accession + "_lifted.bed.gz", ecl.local_path)
    n_ann = peaks.groupby("experiment").size()
    mism, raw_min = [], {"signalValue": np.inf, "pValue": np.inf}
    for r in ecl.itertuples():
        b = pd.read_csv(ROOT / r.peak_path, sep="\t", header=None, names=NP_COLS, usecols=range(10))
        b = b[b.chrom.isin(CANONICAL)]
        if len(b) != n_ann.get(r.experiment, 0):
            mism.append((r.experiment, len(b), int(n_ann.get(r.experiment, 0))))
        for c in raw_min:
            raw_min[c] = min(raw_min[c], b[c].min())
    check("every eCLIP experiment is annotated", set(n_ann.index) == set(ecl.experiment), f"{len(n_ann)} experiments")
    check("annotated peak count = raw canonical-chromosome BED rows, per experiment", not mism,
          f"mismatches: {mism[:3]}" if mism else f"{len(peaks)} peaks")
    check("ENCODE thresholds hold: log2FC >= 3 and -log10 p >= 3",
          peaks.signalValue.min() >= 3 and peaks.pValue.min() >= 3 and min(raw_min.values()) >= 3,
          f"min log2FC {peaks.signalValue.min():.2f}, min -log10p {peaks.pValue.min():.2f}")
    check("peaks: canonical chromosomes, strand +/-, start < end",
          peaks.chrom.isin(CANONICAL).all() and peaks.strand.isin(["+", "-"]).all() and (peaks.start < peaks.end).all())
    check("peaks: valid region labels", peaks.region.astype(str).isin(REGIONS).all())
    inter = peaks.region.astype(str) == "intergenic"
    check("peaks: gene_id missing exactly for intergenic peaks", (peaks.gene_id.isna() == inter).all(),
          f"{100 * inter.mean():.1f}% intergenic")
    check("peaks: no duplicate coordinates within an experiment",
          not peaks.duplicated(["experiment", "chrom", "start", "end", "strand"]).any())
    qc = pd.read_csv(TAB / "peak_qc.tsv", sep="\t")
    check("peak_qc.tsv agrees with peaks_annotated", qc.n_peaks.sum() == len(peaks) and len(qc) == len(n_ann))

    # ---- expression -------------------------------------------------------------------------
    expr = pd.read_parquet(PROC / "expression_by_cell.parquet")
    check("expression: both cell lines present", set(expr.cell_line) == {"K562", "HepG2"})
    ids = {"GAPDH": "ENSG00000111640", "ACTB": "ENSG00000075624", "ALB": "ENSG00000163631",
           "HBG1": "ENSG00000213934", "HBG2": "ENSG00000196565"}
    expr["gid"] = expr.gene_id.str.split(".").str[0]
    ex = expr.set_index(["gid", "cell_line"]).expressed

    def on(g, c):
        return bool(ex.get((ids[g], c), False))
    check("housekeeping genes expressed in both lines (GAPDH, ACTB)",
          all(on(g, c) for g in ["GAPDH", "ACTB"] for c in ["K562", "HepG2"]))
    check("cell identity: ALB expressed in HepG2 only", on("ALB", "HepG2") and not on("ALB", "K562"))
    check("cell identity: HBG1/2 expressed in K562 only",
          all(on(g, "K562") and not on(g, "HepG2") for g in ["HBG1", "HBG2"]))

    # ---- experiment metadata ----------------------------------------------------------------
    meta = pd.read_parquet(PROC / "experiment_meta.parquet")
    check("experiment_meta covers every eCLIP experiment", set(ecl.experiment) <= set(meta.experiment))

    n_fail = sum(not ok for _, ok, _ in results)
    print(f"\n{len(results) - n_fail}/{len(results)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
