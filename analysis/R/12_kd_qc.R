# Binding vs expression, step 1: knockdown RNA-seq QC, nuclear TPM, and the analysis set (COMPUTE).
#
# Sample QC  : every shRNA knockdown and control replicate (RSEM gene files, ENSG genes only):
#              assigned counts, genes with >= 10 counts, chrM share of counts, nuclear TPM (chrM genes
#              removed, renormalised to 1e6) and genes at nuclear TPM >= 1, share of TPM in genes < 500 nt,
#              lowest replicate correlation (log1p nuclear TPM). Flag: the small-RNA-rich control batches
#              (> SMALL_RNA_FLAG of TPM in genes < 500 nt) and their knockdowns.
# Expression : median nuclear TPM per gene x cell line over all control samples (same samples and layout
#              as preprocess/R/03's expression_by_cell.parquet); expressed = median >= 1 or small ncRNA.
# Analysis set (decisions 2026-11-20):
#   - one shRNA knockdown per RBP x cell line (analysis/R/10);
#   - binding set (171): the knockdown <-> eCLIP pairing and >= MIN_BOUND bound-gene rule of analysis/R/11, computed
#     here (section 3a): eCLIP experiment = the one with most gene-assigned peaks for the same RBP x cell line; bound
#     genes counted among genes tested by DESeq2 for the knockdown and expressed by the whole-cell TPM rule
#     (expression_by_cell.parquet), as in analysis/R/11 (column n_bound_aim3, name kept for the table schema);
#   - group = Aim 1 main binding category; an ncRNA-category RBP with >= REASSIGN_GENES bound mRNA genes
#     and < REASSIGN_NC_PCT % of peaks on ncRNA moves to its second category (HLTF, SRSF9 HepG2 -> CDS);
#   - ncRNA binders reported separately; weak knockdowns (efficiency < MIN_KD) excluded;
#   - bound-gene tiers (genes bound in the tested universe, nuclear-TPM expressed): < 300, 300-999, >= 1000.
#   Exclusion order: too few bound genes -> ncRNA binder -> weak knockdown.
#
# Outputs: analysis/results/tables/kd_sample_qc.tsv   one row per sample
#          analysis/results/tables/kd_groups.tsv      one row per selected knockdown (184)
#          data/analysis/expression_nuclear.parquet   per gene x cell line
#          data/analysis/kd_qc_plotdata.rds           for analysis/R/12_kd_qc_figures.R
# Usage  : Rscript analysis/R/12_kd_qc.R   (~2 min; needs analysis/R/05 and 10 outputs)

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
MIN_KD <- 0.5; TIERS <- c(300, 1000); REASSIGN_GENES <- 300; REASSIGN_NC_PCT <- 10
MIN_BOUND <- 30                                                                         # as analysis/R/11
SHORT_NT <- 500; SMALL_RNA_FLAG <- 0.4
CLASSES <- c("3UTR", "CDS", "5UTR", "proximal_intron", "distal_intron", "noncoding")   # as analysis/R/11
NC_REGIONS <- c("small_ncRNA", "noncoding_exon")

man   <- arrow::read_parquet("data/processed/manifest.parquet")
pairs <- arrow::read_parquet("data/processed/kd_pairs.parquet")
excl  <- read.delim("data/processed/excluded.tsv")
ksel  <- read.delim("analysis/results/tables/kd_selected.tsv")
ksum  <- read.delim("analysis/results/tables/kd_summary.tsv")
cats  <- read.delim("analysis/results/tables/rbp_categories.tsv")
expr0 <- arrow::read_parquet("data/processed/expression_by_cell.parquet")
feat  <- readRDS("data/processed/annot_features.rds")
chrM  <- unique(feat$gene_id[as.character(GenomicRanges::seqnames(feat)) == "chrM"])
stopifnot(length(chrM) == 37)

# ---- 1. per-sample QC ----------------------------------------------------------------------------------
files <- man |> filter(category %in% c("kd_counts", "ctrl_counts")) |>
  dplyr::select(accession, category, experiment, target, cell_line, replicate, local_path)
shrna_kd <- pairs$kd_experiment[pairs$assay_title == "shRNA RNA-seq"]
ctrl_of <- setNames(pairs$control_experiment, pairs$kd_experiment)
shrna_ctl <- unique(ctrl_of[shrna_kd])

message("Reading ", nrow(files), " RSEM gene files ...")
ntpm <- list(); qc <- vector("list", nrow(files))
for (i in seq_len(nrow(files))) {
  d <- read.delim(files$local_path[i], colClasses = c(gene_id = "character"))[, c("gene_id", "length", "expected_count", "TPM")]
  d <- d[startsWith(d$gene_id, "ENSG"), ]
  mt <- d$gene_id %in% chrM
  nt <- ifelse(mt, NA, d$TPM); nt <- nt / sum(nt, na.rm = TRUE) * 1e6
  ntpm[[files$accession[i]]] <- setNames(nt, d$gene_id)
  qc[[i]] <- data.frame(accession = files$accession[i], counts_M = sum(d$expected_count) / 1e6,
                        nuclear_counts_M = sum(d$expected_count[!mt]) / 1e6,
                        genes_c10 = sum(d$expected_count >= 10), genes_ntpm1 = sum(nt >= 1, na.rm = TRUE),
                        chrM_share = sum(d$expected_count[mt]) / sum(d$expected_count),
                        short_tpm_share = sum(d$TPM[d$length < SHORT_NT]) / sum(d$TPM))
}
genes <- names(ntpm[[1]]); stopifnot(all(vapply(ntpm, function(x) identical(names(x), genes), TRUE)))
qc <- bind_rows(qc)
rep_r <- files |> group_by(experiment) |>
  summarise(rep_r_min = { m <- log1p(do.call(cbind, ntpm[accession])); m <- m[!is.na(m[, 1]), , drop = FALSE]
                          cc <- cor(m); min(cc[upper.tri(cc)]) }, .groups = "drop")
samples <- files |> left_join(qc, by = "accession") |> left_join(rep_r, by = "experiment") |>
  mutate(kind = ifelse(category == "kd_counts", "knockdown", "control"),
         control_batch = ifelse(kind == "control", experiment, unname(ctrl_of[experiment])),
         shrna = (kind == "knockdown" & experiment %in% shrna_kd) | (kind == "control" & experiment %in% shrna_ctl),
         selected = experiment %in% ksel$kd_experiment[ksel$selected] |
           (kind == "control" & experiment %in% ksel$control_experiment[ksel$selected]))
odd <- samples |> filter(kind == "control", shrna) |> group_by(control_batch) |>
  summarise(s = median(short_tpm_share), .groups = "drop") |> filter(s > SMALL_RNA_FLAG) |> pull(control_batch)
samples <- samples |> mutate(qc_flag = ifelse(control_batch %in% odd, "small-RNA-rich batch", ""))

# ---- 2. nuclear-TPM expression (same control samples as preprocess/R/03) -------------------------------
ctl <- files |> filter(category == "ctrl_counts")
expr <- bind_rows(lapply(split(ctl, ctl$cell_line), function(x) {
  m <- do.call(cbind, ntpm[x$accession])
  data.frame(gene_id = genes, cell_line = x$cell_line[1], median_tpm = apply(m, 1, median), n_samples = ncol(m))
})) |> filter(!is.na(median_tpm)) |>
  left_join(expr0 |> dplyr::select(gene_id, cell_line, small_ncRNA, expressed_old = expressed), by = c("gene_id", "cell_line")) |>
  mutate(small_ncRNA = coalesce(small_ncRNA, FALSE), expressed = median_tpm >= 1 | small_ncRNA)
expr_change <- expr |> group_by(cell_line) |>
  summarise(gained = sum(expressed & !expressed_old, na.rm = TRUE), lost = sum(!expressed & expressed_old, na.rm = TRUE),
            gained_median_tpm = median(median_tpm[expressed & !expressed_old]),   # computed before the counts below
            n_expressed_old = sum(expressed_old, na.rm = TRUE), n_expressed_new = sum(expressed), .groups = "drop")
arrow::write_parquet(expr |> dplyr::select(gene_id, cell_line, median_tpm, n_samples, expressed, small_ncRNA),
                     "data/analysis/expression_nuclear.parquet")

# ---- 3. per-knockdown table: efficiency, binding depth, group, tier, exclusion ---------------------------
med_ntpm <- function(exp, gid) median(vapply(files$accession[files$experiment == exp],
                                             function(a) unname(ntpm[[a]][gid]), 0))
kd <- ksum |> dplyr::select(kd_experiment, target, cell_line, control_experiment, kd_efficiency, n_genes, n_de) |>
  left_join(ksel |> dplyr::select(kd_experiment, target_gene_id), by = "kd_experiment") |>
  rowwise() |>
  mutate(ntpm_ratio = (med_ntpm(kd_experiment, target_gene_id) + 0.5) / (med_ntpm(control_experiment, target_gene_id) + 0.5)) |>
  ungroup() |>
  mutate(kd_efficiency_ntpm = 1 - ntpm_ratio, n_kd_same_control = as.integer(table(control_experiment)[control_experiment]))

# ---- 3a. knockdown <-> eCLIP pairing and binding set (rule of analysis/R/11) ------------------------------
peaks_all <- arrow::read_parquet("data/processed/peaks_annotated.parquet",
                                 col_select = c("experiment", "target", "cell_line", "gene_id", "region")) |>
  mutate(region = as.character(region))
eclip <- peaks_all |> filter(!is.na(gene_id)) |> count(experiment, target, cell_line) |> group_by(target, cell_line) |>
  slice_max(n, n = 1, with_ties = FALSE) |> ungroup()
de_genes <- arrow::read_parquet("data/analysis/de_results.parquet", col_select = c("kd_experiment", "gene_id"))
pair <- kd |> dplyr::select(kd_experiment, target, cell_line) |>
  inner_join(eclip |> dplyr::select(target, cell_line, eclip_experiment = experiment), by = c("target", "cell_line"))
n_bound_wc <- de_genes |> inner_join(pair |> dplyr::select(kd_experiment, cell_line, eclip_experiment), by = "kd_experiment") |>
  semi_join(expr0 |> filter(expressed), by = c("gene_id", "cell_line")) |>              # whole-cell TPM rule, as R/11
  semi_join(peaks_all |> filter(!is.na(gene_id)) |> distinct(experiment, gene_id),
            by = c(eclip_experiment = "experiment", "gene_id")) |>
  count(kd_experiment, name = "n_bound_aim3")
pair <- pair |> inner_join(n_bound_wc, by = "kd_experiment") |> filter(n_bound_aim3 >= MIN_BOUND) |>
  dplyr::select(kd_experiment, eclip_experiment, n_bound_aim3)
stopifnot(!anyDuplicated(pair$kd_experiment))

peaks <- peaks_all |> filter(experiment %in% pair$eclip_experiment) |> dplyr::select(experiment, gene_id, region)
pk <- peaks |> group_by(experiment) |>
  summarise(n_peaks = n(), pct_peaks_ncRNA = 100 * mean(region %in% NC_REGIONS), .groups = "drop")
gene_bind <- peaks |> filter(!is.na(gene_id)) |>
  mutate(cls = ifelse(region %in% NC_REGIONS, "noncoding", region)) |>
  group_by(experiment, gene_id) |> summarise(class = CLASSES[min(match(cls, CLASSES))], .groups = "drop")
univ <- de_genes |> inner_join(kd |> dplyr::select(kd_experiment, cell_line), by = "kd_experiment") |>
  semi_join(expr |> filter(expressed), by = c("gene_id", "cell_line"))
bound <- univ |> inner_join(pair |> dplyr::select(kd_experiment, eclip_experiment), by = "kd_experiment") |>
  inner_join(gene_bind, by = c(eclip_experiment = "experiment", "gene_id")) |>
  group_by(kd_experiment) |>
  summarise(n_bound = n(), n_bound_mrna = sum(class != "noncoding"), .groups = "drop")

cat_of <- cats |> dplyr::select(experiment, aim1_category = main_category, second_group)
tier_lab <- c(sprintf("< %d", TIERS[1]), sprintf("%d–%d", TIERS[1], TIERS[2] - 1), sprintf("≥ %s", format(TIERS[2], big.mark = ",")))
kd <- kd |>
  left_join(pair, by = "kd_experiment") |>
  left_join(pk, by = c(eclip_experiment = "experiment")) |>
  left_join(bound, by = "kd_experiment") |>
  left_join(cat_of, by = c(eclip_experiment = "experiment")) |>
  mutate(in_binding_set = !is.na(eclip_experiment),
         reassigned = in_binding_set & coalesce(aim1_category == "ncRNA", FALSE) &
           n_bound_mrna >= REASSIGN_GENES & pct_peaks_ncRNA < REASSIGN_NC_PCT,
         group = case_when(!in_binding_set ~ NA_character_, reassigned ~ second_group,
                           is.na(aim1_category) ~ "uncategorised", TRUE ~ aim1_category),
         tier = ifelse(in_binding_set, as.character(cut(n_bound, c(0, TIERS, Inf), right = FALSE, labels = tier_lab)), NA),
         weak_kd = kd_efficiency < MIN_KD,
         excluded_reason = case_when(!in_binding_set ~ "too few bound genes", group == "ncRNA" ~ "ncRNA binder",
                                     weak_kd ~ "weak knockdown", TRUE ~ "none"),
         analysed = excluded_reason == "none",
         flag = ifelse(control_experiment %in% odd, "small-RNA-rich batch", "")) |>
  arrange(cell_line, target)
stopifnot(sum(kd$in_binding_set) == nrow(pair), !any(kd$reassigned & kd$group == "ncRNA"))

# flow of the analysis set (all counts computed)
n_found <- sum(pairs$assay_title == "shRNA RNA-seq") + sum(excl$kind == "shRNA RNA-seq")
flow <- tibble(
  step = c("shRNA knockdowns found", "with hg19 gene quantifications", "one per RBP × cell line",
           "with enough eCLIP binding", "not ncRNA binders", "knockdown ≥ 50%"),
  n = c(n_found, sum(pairs$assay_title == "shRNA RNA-seq"), nrow(kd), sum(kd$in_binding_set),
        sum(kd$in_binding_set & kd$group != "ncRNA"), sum(kd$analysed)),
  dropped_reason = c(NA, "no hg19 / GENCODE V19 quantifications", "weaker duplicate knockdowns",
                     sprintf("< 30 bound genes (Aim 2 rule)"), "ncRNA binders (reported separately)",
                     "weak knockdown (< 50% of target RNA removed)")) |>
  mutate(dropped = lag(n) - n)

# ---- 4. write --------------------------------------------------------------------------------------------
r4 <- function(x) round(x, 4)
write.table(samples |> dplyr::select(accession, experiment, kind, target, cell_line, replicate, control_batch, shrna, selected,
                                     counts_M, nuclear_counts_M, genes_c10, genes_ntpm1, chrM_share, short_tpm_share,
                                     rep_r_min, qc_flag) |> mutate(across(where(is.double), r4)),
            "analysis/results/tables/kd_sample_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(kd |> dplyr::select(kd_experiment, target, cell_line, control_experiment, n_kd_same_control, kd_efficiency,
                                kd_efficiency_ntpm, n_genes, n_de, eclip_experiment, n_peaks, pct_peaks_ncRNA, n_bound_aim3,
                                n_bound, n_bound_mrna, aim1_category, second_group, reassigned, group, tier, weak_kd,
                                excluded_reason, analysed, flag) |> mutate(across(where(is.double), r4)),
            "analysis/results/tables/kd_groups.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
saveRDS(list(samples = samples, kd = kd, flow = flow, expr_change = expr_change, odd_batches = odd,
             params = list(MIN_KD = MIN_KD, TIERS = TIERS, tier_lab = tier_lab, REASSIGN_GENES = REASSIGN_GENES,
                           REASSIGN_NC_PCT = REASSIGN_NC_PCT, SHORT_NT = SHORT_NT)),
        "data/analysis/kd_qc_plotdata.rds")

# ---- 5. summary -------------------------------------------------------------------------------------------
s <- samples |> filter(shrna)
cat(sprintf("\nshRNA samples: %d (%d knockdown, %d control replicates)\n", nrow(s), sum(s$kind == "knockdown"), sum(s$kind == "control")))
print(s |> group_by(cell_line) |> summarise(reads_M_min = min(counts_M), reads_M_med = median(counts_M), reads_M_max = max(counts_M),
                                            genes_c10_min = min(genes_c10), genes_c10_max = max(genes_c10),
                                            genes_ntpm1_min = min(genes_ntpm1), chrM_med = median(chrM_share),
                                            rep_r_min = min(rep_r_min)) |> as.data.frame(), digits = 3)
cat("Small-RNA-rich control batches:", paste(odd, collapse = ", "), "\n")
cat("Their knockdowns:", paste(unique(s$experiment[s$qc_flag != "" & s$kind == "knockdown"]), collapse = ", "),
    "| selected:", paste(unique(s$experiment[s$qc_flag != "" & s$kind == "knockdown" & s$selected]), collapse = ", "), "\n")
cat("\nExpressed genes, old (TPM) vs nuclear TPM:\n"); print(as.data.frame(expr_change))
cat(sprintf("\nKnockdown efficiency: DESeq2 vs nuclear TPM, Spearman %.3f; weak (< %.0f%%): %d of %d\n",
            cor(kd$kd_efficiency, kd$kd_efficiency_ntpm, method = "spearman"), 100 * MIN_KD, sum(kd$weak_kd), nrow(kd)))
cat("Reassigned:", paste(kd$target[kd$reassigned], kd$cell_line[kd$reassigned], "->", kd$group[kd$reassigned], collapse = "; "), "\n")
cat(sprintf("Bound genes (nuclear universe) vs Aim 2 count: Spearman %.4f, median change %+.0f\n",
            cor(kd$n_bound, kd$n_bound_aim3, method = "spearman", use = "complete.obs"),
            median(kd$n_bound - kd$n_bound_aim3, na.rm = TRUE)))
cat("\nFlow:\n"); print(as.data.frame(flow))
cat("\nAnalysed knockdowns by group x tier:\n")
print(with(kd |> filter(analysed), table(group, factor(tier, tier_lab))))
print(with(kd |> filter(analysed), table(cell_line, factor(tier, tier_lab))))
message("Wrote kd_sample_qc.tsv, kd_groups.tsv, expression_nuclear.parquet, kd_qc_plotdata.rds")
