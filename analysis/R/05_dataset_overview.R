# Aim 1, part A: dataset overview (COMPUTE ONLY; figures in 05_dataset_overview_figures.R).
#
# Questions: what was profiled; how many peaks per experiment; what peaks look like; is depth driven
# by batch; how many genes each RBP binds; do many RBPs bind the same genes.
#
# Inputs : data/processed/{manifest,kd_pairs,experiment_meta,peaks_annotated,expression_by_cell}.parquet,
#          preprocess/results/tables/peak_qc.tsv
# Outputs: data/analysis/aim1_overview_plotdata.rds   (everything the figures need)
#          analysis/results/tables/dataset_overview.tsv  (one row per eCLIP experiment)
#          analysis/results/tables/rbp_availability.tsv  (one row per RBP: which data exist)
#          analysis/results/tables/gene_rbp_counts.tsv   (expressed gene x cell line: RBPs and peaks)
#          analysis/results/tables/batch_tests.tsv
#
# Usage: Rscript analysis/R/05_dataset_overview.R   (from the project root)

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
MIN_PEAKS_AIM1 <- 100      # experiments used for binding profiles (Aim 1)
MIN_PEAKS_MODEL <- 1000    # experiments used for sequence models (Aim 3)

man <- arrow::read_parquet("data/processed/manifest.parquet")
kd <- arrow::read_parquet("data/processed/kd_pairs.parquet")
meta <- arrow::read_parquet("data/processed/experiment_meta.parquet")
qc <- read.delim("preprocess/results/tables/peak_qc.tsv")
expr <- arrow::read_parquet("data/processed/expression_by_cell.parquet")
peaks <- arrow::read_parquet("data/processed/peaks_annotated.parquet",
                             col_select = c("experiment", "target", "cell_line", "start", "end", "signalValue",
                                            "pValue", "region", "gene_id", "gene_name"))

# ---- per experiment ------------------------------------------------------------------------------
pk_stats <- peaks |> group_by(experiment) |>
  summarise(median_width = median(end - start), median_log2fc = median(signalValue),
            median_neglog10p = median(pValue), .groups = "drop")
exps <- qc |>
  left_join(meta |> select(experiment, release_year, phase = award_rfa), by = "experiment") |>
  left_join(pk_stats, by = "experiment") |>
  group_by(target, cell_line) |> mutate(n_same_rbp_cell = n()) |> ungroup() |>
  mutate(used_aim1 = n_peaks >= MIN_PEAKS_AIM1, used_models = n_peaks >= MIN_PEAKS_MODEL,
         lifted = coord_source == "GRCh38_liftover") |>
  arrange(desc(n_peaks))
stopifnot(nrow(exps) == sum(startsWith(man$category, "eclip_peaks")), !anyNA(exps$phase))

# ---- what was profiled: data availability per RBP ------------------------------------------------------
ecl <- man |> filter(startsWith(category, "eclip_peaks"))
shrna <- kd |> filter(assay_title == "shRNA RNA-seq")
rbns <- unique(man$target[man$category == "rbns_enrichment"])
avail <- tibble(target = sort(unique(ecl$target))) |>
  mutate(eCLIP_K562 = target %in% ecl$target[ecl$cell_line == "K562"],
         eCLIP_HepG2 = target %in% ecl$target[ecl$cell_line == "HepG2"],
         knockdown_K562 = target %in% shrna$target[shrna$cell_line == "K562"],
         knockdown_HepG2 = target %in% shrna$target[shrna$cell_line == "HepG2"],
         RBNS = target %in% rbns)

# ---- batch: does depth / peak strength differ by ENCODE phase? -----------------------------------------
batch_test <- function(d, metric) {
  w <- wilcox.test(d[[metric]][d$phase == "ENCODE3"], d[[metric]][d$phase == "ENCODE4"])
  tibble(metric = metric, n_ENCODE3 = sum(d$phase == "ENCODE3"), n_ENCODE4 = sum(d$phase == "ENCODE4"),
         median_ENCODE3 = median(d[[metric]][d$phase == "ENCODE3"]),
         median_ENCODE4 = median(d[[metric]][d$phase == "ENCODE4"]), wilcoxon_p = w$p.value)
}
batch <- bind_rows(lapply(c("all", "K562", "HepG2"), function(cl) {
  d <- if (cl == "all") exps else filter(exps, cell_line == cl)
  bind_rows(batch_test(d, "n_peaks"), batch_test(d, "median_log2fc")) |> mutate(cell_line = cl, .before = 1)
}))

# ---- genes: how many RBPs bind each expressed gene ------------------------------------------------------
gene_names <- peaks |> filter(!is.na(gene_id)) |> distinct(gene_id, gene_name)
feat <- readRDS("data/processed/annot_features.rds")             # gene span = first to last annotated base
gene_len <- tibble(gene_id = feat$gene_id, s = BiocGenerics::start(feat), e = BiocGenerics::end(feat)) |>
  group_by(gene_id) |> summarise(gene_span = max(e) - min(s) + 1, .groups = "drop")
bound <- peaks |> filter(!is.na(gene_id)) |>
  group_by(cell_line, gene_id) |>
  summarise(n_rbps = n_distinct(target), n_peaks = n(), .groups = "drop")
n_rbps_cell <- ecl |> group_by(cell_line) |> summarise(n_rbps_profiled = n_distinct(target))
genes <- expr |> filter(expressed, !small_ncRNA) |>                       # measurable, expressed genes
  select(cell_line, gene_id, median_tpm) |>
  left_join(bound, by = c("cell_line", "gene_id")) |>
  mutate(n_rbps = coalesce(n_rbps, 0L), n_peaks = coalesce(n_peaks, 0L)) |>
  left_join(gene_names, by = "gene_id") |>
  inner_join(gene_len, by = "gene_id") |>                       # annotated genes only (peaks can only land there)
  left_join(n_rbps_cell, by = "cell_line") |>
  mutate(frac_rbps = n_rbps / n_rbps_profiled)
rho <- genes |> group_by(cell_line) |>
  summarise(spearman_rbps_vs_tpm = cor(n_rbps, median_tpm, method = "spearman"),
            spearman_rbps_vs_span = cor(n_rbps, gene_span, method = "spearman", use = "complete.obs"),
            # partial association: RBPs ~ log TPM + log span (standardised coefficients)
            beta_log_tpm = coef(lm(scale(n_rbps) ~ scale(log10(median_tpm)) + scale(log10(gene_span))))[2],
            beta_log_span = coef(lm(scale(n_rbps) ~ scale(log10(median_tpm)) + scale(log10(gene_span))))[3],
            pct_unbound = 100 * mean(n_rbps == 0), median_rbps = median(n_rbps), .groups = "drop")

# ---- peak-level values for the distribution panels (kept raw so binning stays a figure choice) -----------
peak_vals <- peaks |> transmute(cell_line, width = end - start, log2fc = signalValue, neglog10p = pValue)

# ---- write -----------------------------------------------------------------------------------------------
w <- function(x, f) write.table(x, file.path("analysis/results/tables", f), sep = "\t", quote = FALSE, row.names = FALSE)
w(exps |> mutate(across(where(is.double), ~ round(.x, 3))), "dataset_overview.tsv")
w(avail, "rbp_availability.tsv")
w(genes |> select(cell_line, gene_id, gene_name, median_tpm, gene_span, n_rbps, n_rbps_profiled, frac_rbps, n_peaks) |>
    mutate(across(c(median_tpm, frac_rbps), ~ round(.x, 4))) |> arrange(cell_line, desc(n_rbps)), "gene_rbp_counts.tsv")
w(batch |> mutate(across(where(is.double), ~ signif(.x, 4))), "batch_tests.tsv")
saveRDS(list(exps = exps, avail = avail, batch = batch, genes = genes, gene_rho = rho, peak_vals = peak_vals,
             thresholds = c(aim1 = MIN_PEAKS_AIM1, models = MIN_PEAKS_MODEL)),
        "data/analysis/aim1_overview_plotdata.rds")

# ---- console summary -------------------------------------------------------------------------------------------
cat(sprintf("%d experiments, %d RBPs, %d peaks; %d lifted over; %d RBP x cell pairs with >1 experiment\n",
            nrow(exps), n_distinct(exps$target), sum(exps$n_peaks), sum(exps$lifted),
            n_distinct(paste(exps$target, exps$cell_line)[exps$n_same_rbp_cell > 1])))
cat(sprintf("Below %d peaks: %d; below %d: %d\n", MIN_PEAKS_AIM1, sum(!exps$used_aim1), MIN_PEAKS_MODEL, sum(!exps$used_models)))
cat("RBP availability (n RBPs):\n"); print(colSums(avail[, -1]))
cat("Batch tests:\n"); print(as.data.frame(batch), digits = 3)
cat("Genes:\n"); print(as.data.frame(rho), digits = 3)
cat("Top hub genes:\n")
print(genes |> group_by(cell_line) |> slice_max(n_rbps, n = 5, with_ties = FALSE) |>
        select(cell_line, gene_name, n_rbps, n_rbps_profiled, median_tpm) |> as.data.frame())
