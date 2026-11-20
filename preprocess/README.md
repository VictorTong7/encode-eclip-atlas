# Shared pre-processing

Download, liftover, region annotation and peak annotation. Every downstream analysis reads these outputs. After any change, re-run the affected steps and both checks below.

## Steps (run from the project root)

| Step | Command | Output (`data/processed/` unless noted) | ~Time |
|---|---|---|---|
| 01 | `.venv/Scripts/python -I preprocess/python/01_build_manifest.py` | `manifest.parquet`, `kd_pairs.parquet`, `excluded.tsv` | 2 min |
| 02 | `.venv/Scripts/python -I preprocess/python/02_download.py` | `data/raw/`, `data/ref/` (902 files, 5.6 GB, MD5-checked) | 15 min |
| 02b | `Rscript preprocess/R/02b_liftover_peaks.R` | `data/raw/eclip_peaks_hg19/*_lifted.bed.gz`, `results/tables/liftover_stats.tsv` | 1 min |
| 03 | `Rscript preprocess/R/03_build_annotation.R` | `annot_features.rds`, `annot_flat.rds`, `expression_by_cell.parquet`, `data/ref/gencode_v19.sqlite`, `results/tables/region_lengths.tsv` | 2 min |
| 04 | `Rscript preprocess/R/04_annotate_peaks.R` | `peaks_annotated.parquet`, `results/tables/{region_counts,region_enrichment,peak_qc}.tsv` | 1 min |
| 05a | `.venv/Scripts/python -I preprocess/python/05a_experiment_meta.py` | `experiment_meta.parquet` | <1 min |
| check | `.venv/Scripts/python -I preprocess/checks/check_preprocess.py` (`--quick` skips MD5) and `Rscript preprocess/checks/check_annotation.R` | PASS/FAIL list, non-zero exit on failure | 2 min |

`results/tables/` above means `preprocess/results/tables/`. The shared helpers are `R/lib_regions.R` (regions, priority, `pick_best()`, `build_landmarks()`), `python/encode_api.py` and `python/rbns.py` (RNA Bind-n-Seq loader).

## Peak filtering parameters

| Stage | Rule | Where |
|---|---|---|
| Peak file | Replicate-merged IDR narrowPeak (≥ 2 technical replicates); latest if several; hg19 native, else GRCh38 lifted over | `01` `pick_idr()` |
| ENCODE thresholds (upstream) | log2 fold enrichment over size-matched input (`signalValue`) ≥ 3 and −log10 p (`pValue`) ≥ 3. Both minimums are exactly 3.00. p is capped at 400 (2.6% of peaks). No summits (`peak` = −1), so peaks are placed at their midpoint | ENCODE |
| Liftover (27 sets) | One-to-one mapping only, canonical chromosome, same strand, width within 10% of the original. Median drop < 1%, max 2.65% (APEX1) | `02b` |
| Chromosomes | chr1–22, X, Y, M | `lib_regions.R` `CANONICAL` |
| Region | Midpoint, strand-aware. Priority CDS > 3′UTR > 5′UTR > small ncRNA > noncoding exon > proximal intron (≤ 500 nt from a splice site) > distal intron; otherwise intergenic. Gene ties go to the higher control TPM | `04`, `lib_regions.R` |
| Transcript set | GENCODE V19 `basic` (protein-coding genes: coding transcripts only), all small-ncRNA transcripts, no pseudogenes | `03` |
| Expressed gene | Median control TPM ≥ 1 per cell line; small ncRNAs are always counted as expressed. Used for region enrichment only, never to drop peaks | `03` |

**No other peak filtering is applied here.** Downstream analyses set their own minimum depth; the binding-landscape and cell-line analyses use experiments with ≥ 100 peaks (231 of 250).

## Current outputs (2026-11-19)
- 902 files; 250 eCLIP experiments (168 RBPs, 223 hg19-native + 27 lifted); 935,151 peaks.
- Peak share by region: CDS 30.8%, distal intron 28.6%, proximal intron 19.7%, 3′UTR 11.6%, 5′UTR 4.2%, intergenic 3.1%, noncoding exon 1.1%, small ncRNA 0.9%.
- No ENCODE audit flags on any eCLIP experiment. 6 (RBP, cell line) pairs have two experiments.

## Known limitations (kept by decision, 2026-11-19)
- **Fragmented peaks:** 16.1% of peaks overlap or directly abut another peak in the same experiment (25.3% are within 50 nt), and 8.9% are narrower than 20 nt. ENCODE's IDR output splits some binding sites into several peaks. Peaks are kept as given; the cell-line comparison treats peaks within 50 nt as the same site.
- **chrY in K562:** 43 peaks, although K562 is female (likely mapping artefacts). **chrM:** 1,638 peaks on mitochondrial transcripts. Both kept.
- The IDR threshold is a hard cut: evidence just below it is discarded, and peak strength above it is not used as a weight.
