"""Fetch release date and lab for every eCLIP experiment in the manifest.

Output: data/processed/experiment_meta.parquet
Usage:  .venv/Scripts/python -I preprocess/python/05a_experiment_meta.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # needed under python -I
import pandas as pd

import encode_api as api

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    man = pd.read_parquet(ROOT / "data/processed/manifest.parquet")
    exps = sorted(man.loc[man.category.str.startswith("eclip"), "experiment"].unique())
    hits = api.search({"type": "Experiment", "accession": exps,
                       "field": ["accession", "date_released", "lab.title", "award.rfa"]})
    meta = pd.DataFrame([{"experiment": h["accession"], "date_released": h.get("date_released"),
                          "lab": h.get("lab", {}).get("title"), "award_rfa": h.get("award", {}).get("rfa")}
                         for h in hits])
    missing = set(exps) - set(meta.experiment)
    if missing:
        print(f"Missing metadata for {len(missing)} experiments: {sorted(missing)[:5]}", file=sys.stderr)
        return 1
    meta["release_year"] = pd.to_datetime(meta.date_released).dt.year
    meta.to_parquet(ROOT / "data/processed/experiment_meta.parquet", index=False)
    print(meta.groupby(["lab", "release_year"]).size().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
