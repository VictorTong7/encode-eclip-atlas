# Aim 1, part B: binding landscape (COMPUTE ONLY; figures in 05_binding_landscape_figures.R).
#
# Profile = per-region log2 enrichment (peaks on expressed genes vs region length), z-scored per region
# across experiments. Z-scoring removes the region-wide length term, so profiles compare RBPs with each
# other (raw length enrichment penalises introns, whose RNA is short-lived).
#
# Grouping (2026-11-19): each experiment gets a MAIN BINDING CATEGORY = the group of the region with the
# highest peak share relative to the average RBP (log2_vs_all_rbps). If the best group beats the next group
# by < MARGIN log2, the experiment is "mixed" (no clear preference). The earlier Ward clusters (k chosen by
# silhouette) are still computed and written, unchanged, because Aim 1b (05b) and Aim 3 (07) use them as a
# "same-cluster" null; they are no longer used for Aim 1 figures (weak structure: see cluster_silhouette.tsv).
#
# Inputs : preprocess/results/tables/{region_enrichment,peak_qc}.tsv
# Outputs: data/analysis/aim1_landscape_plotdata.rds
#          analysis/results/tables/rbp_clusters.tsv      (experiment, cluster, category, PCs, z-scores)
#          analysis/results/tables/rbp_categories.tsv    (category, margin, top regions)
#          analysis/results/tables/cluster_silhouette.tsv, concordance.tsv, category_concordance.tsv,
#          pca_variance.tsv, pca_loadings.tsv
#
# Usage: Rscript analysis/R/05_binding_landscape.R   (from the project root)

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
source("preprocess/R/lib_regions.R")
set.seed(42)
MIN_PEAKS <- 100
MARGIN <- 0.25                                   # log2; below this the top group is not a clear preference
CATEGORY_OF <- c(CDS = "CDS", `3UTR` = "3UTR", `5UTR` = "5UTR", small_ncRNA = "ncRNA", noncoding_exon = "ncRNA",
                 proximal_intron = "splice_site", distal_intron = "deep_intron")
CATEGORIES <- c("CDS", "3UTR", "5UTR", "ncRNA", "splice_site", "deep_intron", "mixed")
w <- function(x, f) write.table(x, file.path("analysis/results/tables", f), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- profile matrix ---------------------------------------------------------------------------------
qc <- read.delim("preprocess/results/tables/peak_qc.tsv")
enr <- read.delim("preprocess/results/tables/region_enrichment.tsv") |>
  inner_join(qc |> select(experiment, n_peaks), by = "experiment") |>
  filter(n_peaks >= MIN_PEAKS)
meta <- distinct(enr, experiment, target, cell_line, coord_source, n_peaks) |>
  mutate(label = paste(target, cell_line))
wide <- enr |> select(experiment, region, log2_enrich) |>
  pivot_wider(names_from = region, values_from = log2_enrich)
X <- as.matrix(wide[, REGIONS]); rownames(X) <- wide$experiment
Z <- scale(X)
meta <- meta[match(rownames(Z), meta$experiment), ]
message(sprintf("%d experiments with >= %d peaks (%d RBPs)", nrow(Z), MIN_PEAKS, n_distinct(meta$target)))

# ---- Ward clusters (unchanged; kept for the 05b / 07 same-cluster nulls) ------------------------------------
D <- dist(Z)
hc <- hclust(D, method = "ward.D2")
sil <- sapply(4:12, function(k) mean(cluster::silhouette(cutree(hc, k), D)[, "sil_width"]))
sil_tab <- data.frame(k = 4:12, mean_silhouette = round(sil, 4))
K <- sil_tab$k[which.max(sil_tab$mean_silhouette)]
cl <- cutree(hc, K)
top2 <- sapply(split(as.data.frame(Z), cl), function(d) {
  r <- REGIONS[order(-colMeans(d))[1:2]]
  paste(c(CDS = "CDS", `3UTR` = "3′UTR", `5UTR` = "5′UTR", small_ncRNA = "Small ncRNA", noncoding_exon = "Noncoding exon",
          proximal_intron = "Proximal intron", distal_intron = "Distal intron")[r], collapse = " + ")
})
size_order <- order(-table(cl))
meta$cluster <- match(cl, size_order)
meta$cluster_name <- setNames(paste0("C", seq_len(K), " ", top2[size_order]), seq_len(K))[as.character(meta$cluster)]

# ---- main binding category ------------------------------------------------------------------------------------
rel <- enr |> select(experiment, region, log2_vs_all_rbps) |>
  pivot_wider(names_from = region, values_from = log2_vs_all_rbps)
R <- as.matrix(rel[, REGIONS]); rownames(R) <- rel$experiment; R <- R[rownames(Z), ]
grp_score <- sapply(setdiff(CATEGORIES, "mixed"), function(g) apply(R[, names(CATEGORY_OF)[CATEGORY_OF == g], drop = FALSE], 1, max))
ord <- t(apply(grp_score, 1, order, decreasing = TRUE))
best <- colnames(grp_score)[ord[, 1]]; second <- colnames(grp_score)[ord[, 2]]
margin <- grp_score[cbind(seq_len(nrow(R)), ord[, 1])] - grp_score[cbind(seq_len(nrow(R)), ord[, 2])]
meta$top_group <- best; meta$second_group <- second
meta$category_margin <- margin
meta$main_category <- factor(ifelse(margin < MARGIN, "mixed", best), CATEGORIES)
meta$top_region <- REGIONS[apply(R, 1, which.max)]

# ---- PCA (replaces the UMAP: linear, so axes and region loadings are interpretable) -------------------------------
pca <- prcomp(Z, center = TRUE, scale. = FALSE)
meta$pc1 <- pca$x[, 1]; meta$pc2 <- pca$x[, 2]
pca_var <- data.frame(pc = colnames(pca$x), sd = pca$sdev, pct_variance = 100 * pca$sdev^2 / sum(pca$sdev^2))
pca_load <- data.frame(region = rownames(pca$rotation), pca$rotation, row.names = NULL)

# ---- experiment x experiment similarity (Pearson correlation of z-score profiles) ----------------------------------
S <- cor(t(Z))
within_order <- unlist(lapply(levels(meta$main_category), function(g) {
  i <- which(meta$main_category == g)
  if (length(i) < 3) return(i)
  i[hclust(as.dist(1 - S[i, i]), method = "average")$order]
}))

# ---- K562 vs HepG2 concordance (unchanged test) + category agreement --------------------------------------------------
both <- intersect(meta$target[meta$cell_line == "K562"], meta$target[meta$cell_line == "HepG2"])
k_rows <- which(meta$cell_line == "K562")[match(both, meta$target[meta$cell_line == "K562"])]
h_rows <- which(meta$cell_line == "HepG2")[match(both, meta$target[meta$cell_line == "HepG2"])]
Zk <- Z[k_rows, , drop = FALSE]; Zh <- Z[h_rows, , drop = FALSE]
cross <- sqrt(outer(rowSums(Zk^2), rowSums(Zh^2), "+") - 2 * Zk %*% t(Zh))
same <- diag(cross); diff_d <- cross[row(cross) != col(cross)]
obs <- mean(same)
perm <- replicate(10000, mean(cross[cbind(seq_along(both), sample(seq_along(both)))]))
p_perm <- (sum(perm <= obs) + 1) / (length(perm) + 1)
pct <- sapply(seq_along(both), function(i) mean(cross[i, -i] <= cross[i, i]))
conc <- data.frame(target = both, distance_same = round(same, 3), percentile_vs_other_rbps = round(pct, 3),
                   cluster_K562 = meta$cluster[k_rows], cluster_HepG2 = meta$cluster[h_rows]) |>
  arrange(distance_same)
cat_conc <- data.frame(target = both, category_K562 = meta$main_category[k_rows], category_HepG2 = meta$main_category[h_rows],
                       corr_same_rbp = round(S[cbind(k_rows, h_rows)], 3))
S_kh <- S[k_rows, h_rows]                        # rows K562, cols HepG2, same RBP order
kh_order <- order(meta$main_category[k_rows], -meta$category_margin[k_rows])

# ---- region composition per category (mean peak share; for the stacked bars) -----------------------------------------
peak_frac <- enr |> select(experiment, region, frac_peaks) |>
  inner_join(meta |> select(experiment, main_category), by = "experiment") |>
  group_by(main_category, region) |> summarise(frac = mean(frac_peaks), .groups = "drop") |>
  group_by(main_category) |> mutate(frac = frac / sum(frac)) |> ungroup() |>
  left_join(meta |> count(main_category, name = "n"), by = "main_category")

# ---- write ---------------------------------------------------------------------------------------------------------------
w(sil_tab, "cluster_silhouette.tsv")
w(meta |> select(experiment, target, cell_line, coord_source, n_peaks, cluster, cluster_name, main_category,
                 category_margin, pc1, pc2) |>
    mutate(across(c(category_margin, pc1, pc2), ~ round(.x, 3))) |>
    bind_cols(as.data.frame(round(Z, 3)) |> rename_with(~ paste0("z_", .x))), "rbp_clusters.tsv")
w(meta |> select(experiment, target, cell_line, n_peaks, main_category, top_group, second_group, category_margin, top_region) |>
    mutate(category_margin = round(category_margin, 3)) |> arrange(main_category, desc(category_margin)), "rbp_categories.tsv")
w(conc, "concordance.tsv")
w(cat_conc, "category_concordance.tsv")
w(pca_var |> mutate(across(where(is.double), ~ round(.x, 3))), "pca_variance.tsv")
w(pca_load |> mutate(across(where(is.double), ~ round(.x, 4))), "pca_loadings.tsv")
saveRDS(list(meta = meta, Z = Z, S = S, within_order = within_order, sil = sil_tab, K = K,
             pca = list(var = pca_var, load = pca_load), peak_frac = peak_frac,
             conc = list(both = both, same = same, diff = diff_d, obs = obs, p_perm = p_perm, table = conc,
                         category = cat_conc, S_kh = S_kh, kh_order = kh_order),
             margin = MARGIN),
        "data/analysis/aim1_landscape_plotdata.rds")

# ---- console summary -------------------------------------------------------------------------------------------------------
message(sprintf("Ward k = %d (mean silhouette %.3f; range over k = 4-12: %.3f-%.3f)", K, max(sil), min(sil), max(sil)))
cat("\nMain binding category (margin >= ", MARGIN, " log2):\n", sep = ""); print(table(meta$main_category))
cat(sprintf("Before the mixed rule: %s\n", paste(names(table(best)), table(best), sep = " ", collapse = ", ")))
cat(sprintf("Same category in K562 and HepG2: %d/%d (Ward cluster: %d/%d)\n",
            sum(cat_conc$category_K562 == cat_conc$category_HepG2), nrow(cat_conc),
            sum(conc$cluster_K562 == conc$cluster_HepG2), nrow(conc)))
cat(sprintf("PCA variance: %s\n", paste(sprintf("%s %.1f%%", pca_var$pc, pca_var$pct_variance), collapse = ", ")))
cat(sprintf("Concordance: mean same-RBP distance %.2f vs %.2f different (permutation p = %.4g)\n", obs, mean(diff_d), p_perm))
cat(sprintf("Profile correlation, same RBP across lines: median %.2f; different RBPs: median %.2f\n",
            median(S[cbind(k_rows, h_rows)]), median(S_kh[row(S_kh) != col(S_kh)])))
