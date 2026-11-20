"""Build the download manifest for the ENCODE eCLIP atlas (hg19 build).

Outputs (data/processed/):
  manifest.parquet  one row per file to download
  kd_pairs.parquet  knockdown experiment <-> matched control <-> RBP <-> cell line
  excluded.tsv      experiments dropped, with the reason

Usage:  .venv/Scripts/python -I preprocess/python/01_build_manifest.py [--force]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # needed under python -I
import pandas as pd

import encode_api as api

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "processed"
CELL_LINES = ["K562", "HepG2"]
KD_ASSAYS = ["shRNA RNA-seq", "CRISPRi RNA-seq"]

CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
CHAIN_MD5_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/md5sum.txt"
REFERENCES = {  # ENCODE file id -> (category, local folder)
    "ENCFF590SSX": ("genome", "data/ref"),  # hg19.fa.gz
    "gencode.v19.annotation": ("annotation", "data/ref"),  # GENCODE V19 GTF
}

# Expected counts from the 2026-11-16 API survey; used to catch large changes at ENCODE.
EXPECTED = {
    "eclip_peaks_hg19": 223,
    "eclip_peaks_grch38": 27,
    "kd_counts": 387,
    "ctrl_counts": 94,
    "genome": 1,
    "annotation": 1,
    "chain": 1,
    "rbns_enrichment": 168,
}

FILE_FIELDS = ["accession", "dataset", "assembly", "file_type", "file_format", "output_type",
               "technical_replicates", "biological_replicates", "file_size", "md5sum", "href",
               "date_created", "genome_annotation", "status"]


def exp_id(f: dict) -> str:
    return f["dataset"].strip("/").split("/")[-1]


def pick_idr(files: list[dict]) -> dict[str, dict]:
    """Per experiment, the replicate-merged IDR peak file (2 technical replicates); latest if several."""
    best: dict[str, dict] = {}
    for f in files:
        if len(f.get("technical_replicates", [])) < 2:
            continue  # per-replicate, unfiltered peaks
        e = exp_id(f)
        if e not in best or f["date_created"] > best[e]["date_created"]:
            best[e] = f
    return best


def row(f: dict, category: str, folder: str, **extra) -> dict:
    href = f["href"]
    name = href.rsplit("/", 1)[-1]
    return {
        "accession": f.get("accession") or f.get("title") or name.split(".")[0],
        "category": category,
        "experiment": exp_id(f) if f.get("dataset", "").startswith("/experiments/") else None,
        "assembly": f.get("assembly"),
        "href": href if href.startswith("http") else api.BASE + href,
        "file_size": int(f.get("file_size") or 0),
        "md5sum": f.get("md5sum"),
        "local_path": f"{folder}/{name}",
        **extra,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="write the manifest even if counts differ from expected")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows, excluded = [], []

    # ---- 1. eCLIP experiments and peak files --------------------------------------
    print("Querying eCLIP experiments ...")
    eclip = api.search({"type": "Experiment", "assay_title": "eCLIP",
                        "biosample_ontology.term_name": CELL_LINES,
                        "field": ["accession", "target.label", "biosample_ontology.term_name",
                                  "date_released", "audit"]})
    emeta = {x["accession"]: x for x in eclip}
    peak_q = {"type": "File", "assay_title": "eCLIP", "file_type": "bed narrowPeak",
              "biosample_ontology.term_name": CELL_LINES, "field": FILE_FIELDS}
    hg19 = pick_idr(api.search({**peak_q, "assembly": "hg19"}))
    grch38 = pick_idr(api.search({**peak_q, "assembly": "GRCh38"}))

    for acc, x in sorted(emeta.items()):
        common = dict(target=x["target"]["label"], cell_line=x["biosample_ontology"]["term_name"],
                      assay_title="eCLIP", audit_flags=api.audit_summary(x))
        if acc in hg19:
            rows.append(row(hg19[acc], "eclip_peaks_hg19", "data/raw/eclip_peaks_hg19",
                            needs_liftover=False, coord_source="hg19_native", **common))
        elif acc in grch38:
            rows.append(row(grch38[acc], "eclip_peaks_grch38", "data/raw/eclip_peaks_grch38",
                            needs_liftover=True, coord_source="GRCh38_liftover", **common))
        else:
            excluded.append({"experiment": acc, "kind": "eCLIP", "target": common["target"],
                             "cell_line": common["cell_line"], "reason": "no IDR peak file in hg19 or GRCh38"})
    eclip_pairs = {(x["target"]["label"], x["biosample_ontology"]["term_name"]) for x in eclip}

    # ---- 2. Knockdown experiments, matched controls, gene counts ------------------
    print("Querying knockdown RNA-seq experiments ...")
    kd = api.search({"type": "Experiment", "assay_title": KD_ASSAYS,
                     "biosample_ontology.term_name": CELL_LINES,
                     "field": ["accession", "assay_title", "target.label", "biosample_ontology.term_name",
                               "possible_controls.accession", "audit"]})
    kd = [x for x in kd if "target" in x
          and (x["target"]["label"], x["biosample_ontology"]["term_name"]) in eclip_pairs]
    ctrl_ids = {c["accession"] for x in kd for c in x.get("possible_controls", [])}
    ctrl_meta = {x["accession"]: x for x in api.search(
        {"type": "Experiment", "accession": sorted(ctrl_ids),
         "field": ["accession", "assay_title", "biosample_ontology.term_name", "audit"]})} if ctrl_ids else {}

    print("Querying hg19 / GENCODE V19 gene quantification files ...")
    gq = api.search_by_dataset([x["accession"] for x in kd] + sorted(ctrl_ids),
                               {"output_type": "gene quantifications", "assembly": "hg19",
                                "genome_annotation": "V19", "field": FILE_FIELDS})
    by_exp: dict[str, list[dict]] = {}
    for f in gq:
        by_exp.setdefault(exp_id(f), []).append(f)

    pairs, used_ctrls = [], set()
    for x in sorted(kd, key=lambda z: z["accession"]):
        acc, tgt, cell = x["accession"], x["target"]["label"], x["biosample_ontology"]["term_name"]
        ctrls = [c["accession"] for c in x.get("possible_controls", [])]
        reason = None
        if len(ctrls) != 1:
            reason = f"expected 1 matched control, found {len(ctrls)}"
        elif acc not in by_exp:
            reason = "no hg19/V19 gene quantifications for knockdown"
        elif ctrls[0] not in by_exp:
            reason = f"no hg19/V19 gene quantifications for control {ctrls[0]}"
        if reason:
            excluded.append({"experiment": acc, "kind": x["assay_title"], "target": tgt,
                             "cell_line": cell, "reason": reason})
            continue
        ctl = ctrls[0]
        used_ctrls.add(ctl)
        pairs.append({"kd_experiment": acc, "control_experiment": ctl, "target": tgt, "cell_line": cell,
                      "assay_title": x["assay_title"], "kd_audit": api.audit_summary(x),
                      "ctrl_audit": api.audit_summary(ctrl_meta.get(ctl, {}))})
        for f in by_exp[acc]:
            rows.append(row(f, "kd_counts", "data/raw/kd_counts", target=tgt, cell_line=cell,
                            assay_title=x["assay_title"], control_experiment=ctl,
                            replicate=",".join(map(str, f.get("biological_replicates", []))),
                            audit_flags=api.audit_summary(x)))
    for ctl in sorted(ctrl_ids - used_ctrls):
        if ctl not in by_exp:
            cm = ctrl_meta.get(ctl, {})
            excluded.append({"experiment": ctl, "kind": "control", "target": None,
                             "cell_line": cm.get("biosample_ontology", {}).get("term_name"),
                             "reason": "no hg19/V19 gene quantifications"})
    for ctl in sorted(used_ctrls):
        cm = ctrl_meta.get(ctl, {})
        for f in by_exp[ctl]:
            rows.append(row(f, "ctrl_counts", "data/raw/ctrl_counts",
                            cell_line=cm.get("biosample_ontology", {}).get("term_name"),
                            assay_title=cm.get("assay_title"),
                            replicate=",".join(map(str, f.get("biological_replicates", []))),
                            audit_flags=api.audit_summary(cm)))

    # ---- 2b. RNA Bind-n-Seq (in vitro) enrichment for RBPs that also have eCLIP -----
    print("Querying RNA Bind-n-Seq experiments ...")
    eclip_rbps = {t for t, _ in eclip_pairs}
    rbns = [x for x in api.search({"type": "Experiment", "assay_title": "RNA Bind-n-Seq",
                                   "field": ["accession", "target.label"]})
            if x.get("target", {}).get("label") in eclip_rbps]
    rbns_target = {x["accession"]: x["target"]["label"] for x in rbns}
    for f in api.search_by_dataset(list(rbns_target), {"output_type": "enrichment",
                                                       "field": FILE_FIELDS + ["submitted_file_name"]}):
        rows.append(row(f, "rbns_enrichment", "data/raw/rbns", assay_title="RNA Bind-n-Seq",
                        target=rbns_target[exp_id(f)],
                        file_label=f.get("submitted_file_name", "").split("/")[-1]))

    # ---- 3. Reference files --------------------------------------------------------
    print("Adding reference files ...")
    for fid, (cat, folder) in REFERENCES.items():
        f = api.get(f"/files/{fid}/")
        rows.append(row({**f, "dataset": ""}, cat, folder, accession=fid))
    import requests  # UCSC chain file; md5 from UCSC's md5sum.txt
    md5_txt = requests.get(CHAIN_MD5_URL, timeout=60).text
    chain_md5 = next(l.split()[0] for l in md5_txt.splitlines() if l.endswith("hg38ToHg19.over.chain.gz"))
    head = requests.head(CHAIN_URL, timeout=60, allow_redirects=True)
    rows.append({"accession": "hg38ToHg19.over.chain", "category": "chain", "href": CHAIN_URL,
                 "file_size": int(head.headers.get("content-length", 0)), "md5sum": chain_md5,
                 "local_path": "data/ref/hg38ToHg19.over.chain.gz", "assembly": "GRCh38->hg19"})

    # ---- 4. Write outputs and summary ----------------------------------------------
    man = pd.DataFrame(rows)
    cols = ["accession", "category", "experiment", "assay_title", "target", "cell_line", "replicate",
            "control_experiment", "assembly", "needs_liftover", "coord_source", "href", "file_size",
            "md5sum", "audit_flags", "file_label", "local_path"]
    man = man.reindex(columns=cols)
    man["needs_liftover"] = man["needs_liftover"].astype("boolean").fillna(False)
    dup = man["local_path"].duplicated()
    if dup.any():
        raise SystemExit(f"Duplicate local paths in manifest: {man.loc[dup, 'local_path'].tolist()[:5]}")
    if man["md5sum"].isna().any():
        raise SystemExit(f"Files without md5: {man.loc[man.md5sum.isna(), 'accession'].tolist()}")

    summary = man.groupby("category").agg(files=("accession", "size"), gb=("file_size", lambda s: s.sum() / 1e9))
    print("\nManifest summary")
    print(summary.to_string(float_format=lambda v: f"{v:.3f}"))
    print(f"TOTAL  {len(man)} files  {man.file_size.sum() / 1e9:.2f} GB")
    pr = pd.DataFrame(pairs)
    print(f"\neCLIP experiments: {man.category.str.startswith('eclip').sum()} "
          f"({man.target[man.category.str.startswith('eclip')].nunique()} RBPs)")
    print(f"Knockdown pairs: {len(pr)} experiments, {pr[['target', 'cell_line']].drop_duplicates().shape[0]} "
          f"RBP x cell-line pairs, {pr.control_experiment.nunique()} controls")
    print(f"Excluded: {len(excluded)} experiments (see data/processed/excluded.tsv)")

    bad = {c: (int(summary.files.get(c, 0)), n) for c, n in EXPECTED.items()
           if abs(int(summary.files.get(c, 0)) - n) > max(2, 0.1 * n)}
    if bad and not args.force:
        print(f"\nCounts differ substantially from expected (found, expected): {bad}. "
              "Re-run with --force to accept.", file=sys.stderr)
        return 1

    man.to_parquet(OUT / "manifest.parquet", index=False)
    pr.to_parquet(OUT / "kd_pairs.parquet", index=False)
    pd.DataFrame(excluded, columns=["experiment", "kind", "target", "cell_line", "reason"]).to_csv(
        OUT / "excluded.tsv", sep="\t", index=False)
    print(f"\nWrote {OUT / 'manifest.parquet'}, kd_pairs.parquet, excluded.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
