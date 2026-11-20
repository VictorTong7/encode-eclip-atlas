# Assign every eCLIP peak to a transcript region and gene; compute region enrichment.
#
# Inputs : data/processed/manifest.parquet, eCLIP peaks (hg19 native + *_lifted),
#          data/processed/annot_features.rds, expression_by_cell.parquet,
#          preprocess/results/tables/region_lengths.tsv
# Outputs: data/processed/peaks_annotated.parquet   one row per peak
#          preprocess/results/tables/region_counts.tsv          experiment x region counts (all peaks)
#          preprocess/results/tables/region_enrichment.tsv      log2 enrichment (peaks on expressed genes)
#          preprocess/results/tables/peak_qc.tsv
#
# Each peak is placed at its midpoint (strand-aware) and given the highest-priority region it
# overlaps (preprocess/R/lib_regions.R); ties between genes go to the gene with the higher control TPM.
#
# Usage: Rscript preprocess/R/04_annotate_peaks.R   (from the project root)

suppressPackageStartupMessages({ library(GenomicRanges); library(dplyr); library(tidyr) })
source("preprocess/R/lib_regions.R")

feat <- readRDS("data/processed/annot_features.rds")
expr <- arrow::read_parquet("data/processed/expression_by_cell.parquet")
lens <- read.delim("preprocess/results/tables/region_lengths.tsv")
man <- arrow::read_parquet("data/processed/manifest.parquet") |>
  filter(startsWith(category, "eclip_peaks")) |>
  mutate(peak_path = ifelse(needs_liftover,
                            file.path("data/raw/eclip_peaks_hg19", paste0(accession, "_lifted.bed.gz")),
                            local_path))
stopifnot(all(file.exists(man$peak_path)))

np_cols <- c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
annotate_one <- function(i) {
  r <- man[i, ]
  b <- read.delim(r$peak_path, header = FALSE, col.names = np_cols)
  b <- b[b$chrom %in% CANONICAL, ]
  mid <- b$start + (b$end - b$start) %/% 2 + 1L                 # 1-based midpoint
  gr <- GRanges(b$chrom, IRanges(mid, width = 1L), strand = b$strand)
  e <- expr[expr$cell_line == r$cell_line, ]
  tpm_f <- e$median_tpm[match(feat$gene_id, e$gene_id)]
  best <- pick_best(findOverlaps(gr, feat), feat$priority, tpm_f)
  region <- rep("intergenic", length(gr)); gene_id <- rep(NA_character_, length(gr))
  region[best$query] <- feat$region[best$feature]
  gene_id[best$query] <- feat$gene_id[best$feature]
  data.frame(experiment = r$experiment, accession = r$accession, target = r$target,
             cell_line = r$cell_line, coord_source = r$coord_source,
             chrom = b$chrom, start = b$start, end = b$end, strand = b$strand,
             signalValue = b$signalValue, pValue = b$pValue, region = region,
             subtype = ifelse(region == "small_ncRNA", feat$subtype[match(gene_id, feat$gene_id)], NA),
             gene_id = gene_id,
             gene_name = feat$gene_name[match(gene_id, feat$gene_id)],
             gene_type = feat$gene_type[match(gene_id, feat$gene_id)],
             gene_expressed = !is.na(gene_id) & gene_id %in% e$gene_id[e$expressed])
}
message("Annotating ", nrow(man), " peak sets ...")
peaks <- bind_rows(lapply(seq_len(nrow(man)), annotate_one)) |>
  mutate(region = factor(region, REGION_LEVELS))
arrow::write_parquet(peaks, "data/processed/peaks_annotated.parquet")
message(sprintf("  %d peaks annotated", nrow(peaks)))

# ---- counts, QC and enrichment -----------------------------------------------------------
meta <- distinct(peaks, experiment, target, cell_line, coord_source)
counts <- peaks |> count(experiment, region, .drop = FALSE) |>
  pivot_wider(names_from = region, values_from = n, values_fill = 0)
counts <- right_join(meta, counts, by = "experiment")
write.table(counts, "preprocess/results/tables/region_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

qc <- peaks |> group_by(experiment, target, cell_line, coord_source) |>
  summarise(n_peaks = n(), pct_intergenic = 100 * mean(region == "intergenic"),
            pct_expressed_gene = 100 * mean(gene_expressed),
            n_genes = n_distinct(gene_id[!is.na(gene_id)]), .groups = "drop")
write.table(qc, "preprocess/results/tables/peak_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Enrichment uses only peaks on expressed genes, matching the region-length background.
enr <- peaks |> filter(gene_expressed) |>
  mutate(experiment = factor(experiment, meta$experiment)) |>
  count(experiment, region, .drop = FALSE) |>
  mutate(experiment = as.character(experiment), region = as.character(region)) |>
  filter(region != "intergenic") |>
  group_by(experiment) |> mutate(n_total = sum(n)) |> ungroup() |>
  left_join(meta |> select(experiment, cell_line), by = "experiment") |>
  left_join(lens |> group_by(cell_line) |>
              mutate(frac_len = length_nt / sum(length_nt), frac_len_tpm = length_x_tpm / sum(length_x_tpm)) |>
              ungroup() |> select(cell_line, region, frac_len, frac_len_tpm),
            by = c("cell_line", "region")) |>
  mutate(frac_peaks = (n + 0.5) / (n_total + 0.5 * length(REGIONS)),
         log2_enrich = log2(frac_peaks / frac_len),
         log2_enrich_tpm = log2(frac_peaks / frac_len_tpm)) |>
  left_join(meta |> select(experiment, target, coord_source), by = "experiment") |>
  # Peak share relative to the average RBP in the same cell line. Easiest to interpret:
  # length-normalised enrichment penalises introns (short-lived RNA) and inflates small regions.
  group_by(cell_line, region) |> mutate(log2_vs_all_rbps = log2(frac_peaks / mean(frac_peaks))) |> ungroup() |>
  select(experiment, target, cell_line, coord_source, region, n, n_total, frac_peaks,
         log2_vs_all_rbps, log2_enrich, log2_enrich_tpm)
write.table(enr, "preprocess/results/tables/region_enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nPeaks by region (all experiments):\n")
print(round(100 * prop.table(table(peaks$region)), 1))
cat(sprintf("\nMedian per experiment: intergenic %.1f%%, on expressed genes %.1f%%\n",
            median(qc$pct_intergenic), median(qc$pct_expressed_gene)))
cat("Small ncRNA subtypes:\n"); print(sort(table(peaks$subtype), decreasing = TRUE))
