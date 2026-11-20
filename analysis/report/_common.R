# Shared helpers for the analysis website (sourced by every page's setup chunk).
# Every number shown on the site is computed here from saved result tables; nothing is typed by hand.

suppressPackageStartupMessages(library(dplyr))

# project root = the folder holding _quarto.yml (analysis/), found from the page's working directory
.find_root <- function(d = getwd()) {
  while (!file.exists(file.path(d, "_quarto.yml"))) {
    up <- dirname(d); if (up == d) stop("_quarto.yml not found above ", getwd()); d <- up
  }
  d
}
ROOT <- .find_root()

tab <- function(f) {  # analysis tables, falling back to the shared pre-processing tables (e.g. peak_qc)
  p <- file.path(ROOT, "results", "tables", paste0(f, ".tsv"))
  if (!file.exists(p)) p <- file.path(ROOT, "..", "preprocess", "results", "tables", paste0(f, ".tsv"))
  read.delim(p)
}
fig <- function(f) paste0("/results/figures/", f, ".png")   # "/" = project root; Quarto rewrites per page
n0 <- function(x) format(x, big.mark = ",")
r2 <- function(x) sprintf("%.2f", x)
pct <- function(x) sprintf("%.0f%%", 100 * x)
src <- function(path) sprintf("Script: `%s`.", path)        # figure captions name the script that made them

# ---- numbers per aim --------------------------------------------------------------------------------------
aim1_numbers <- function() {
  qc <- tab("peak_qc"); sil <- tab("cluster_silhouette"); ov <- tab("dataset_overview"); av <- tab("rbp_availability")
  bt <- tab("batch_tests"); gr <- tab("gene_rbp_counts"); rc <- tab("rbp_categories"); cc <- tab("category_concordance")
  pv <- tab("pca_variance")
  within(list(), {
    n_exp <- nrow(qc); n_rbp <- n_distinct(qc$target); n_peaks <- sum(qc$n_peaks)
    best_sil <- max(sil$mean_silhouette); min_sil <- min(sil$mean_silhouette)
    n_lift <- sum(ov$lifted); n_lt100 <- sum(ov$n_peaks < 100); n_lt1k <- sum(ov$n_peaks < 1000)
    n_kd_match <- sum((av$eCLIP_K562 & av$knockdown_K562) | (av$eCLIP_HepG2 & av$knockdown_HepG2))
    n_rbns <- sum(av$RBNS)
    k562_enc4 <- bt$median_ENCODE4[bt$cell_line == "K562" & bt$metric == "n_peaks"]
    k562_enc3 <- bt$median_ENCODE3[bt$cell_line == "K562" & bt$metric == "n_peaks"]
    k562_p <- bt$wilcoxon_p[bt$cell_line == "K562" & bt$metric == "n_peaks"]
    rho_tpm <- mean(sapply(split(gr, gr$cell_line), function(g) cor(g$n_rbps, g$median_tpm, method = "spearman")))
    pct_unbound <- 100 * mean(gr$n_rbps == 0)
    n_prof <- nrow(rc); n_mixed <- sum(rc$main_category == "mixed")
    n_both <- nrow(cc); same_cat <- sum(cc$category_K562 == cc$category_HepG2); corr_same <- median(cc$corr_same_rbp)
    pc12 <- sum(pv$pct_variance[1:2]); pc1 <- pv$pct_variance[1]; pc2 <- pv$pct_variance[2]; pc3 <- pv$pct_variance[3]
    median_width <- median(ov$median_width)
    depth_fold <- max(ov$n_peaks) / min(ov$n_peaks)
    # extended-legend details
    peaks_min <- min(ov$n_peaks); peaks_max <- max(ov$n_peaks); peaks_med <- median(ov$n_peaks)
    set_n <- colSums(av[, c("eCLIP_K562", "eCLIP_HepG2", "knockdown_K562", "knockdown_HepG2", "RBNS")])
    n_dup <- sum(!duplicated(paste(ov$target, ov$cell_line)[ov$n_same_rbp_cell > 1]))
    ph <- function(cl, ph) sum(ov$cell_line == cl & ov$phase == ph)
    n_k3 <- ph("K562", "ENCODE3"); n_k4 <- ph("K562", "ENCODE4"); n_h3 <- ph("HepG2", "ENCODE3"); n_h4 <- ph("HepG2", "ENCODE4")
    bp <- function(cl, m) bt$wilcoxon_p[bt$cell_line == cl & bt$metric == m]
    hepg2_p <- bp("HepG2", "n_peaks"); str_p_k <- bp("K562", "median_log2fc"); str_p_h <- bp("HepG2", "median_log2fc")
    med_genes <- median(ov$n_genes); med_ppg <- median(ov$n_peaks / ov$n_genes)
    med_expr <- median(ov$pct_expressed_gene); med_inter <- median(ov$pct_intergenic)
    n_genes_k <- sum(gr$cell_line == "K562"); n_genes_h <- sum(gr$cell_line == "HepG2")
    n_prof_k <- gr$n_rbps_profiled[gr$cell_line == "K562"][1]; n_prof_h <- gr$n_rbps_profiled[gr$cell_line == "HepG2"][1]
    cat_n <- table(factor(rc$main_category, c("CDS", "3UTR", "5UTR", "ncRNA", "splice_site", "deep_intron", "mixed")))
    # K562 vs HepG2 profile distances, as in analysis/R/05_binding_landscape.R (first experiment per RBP x cell line)
    cl <- tab("rbp_clusters"); zc <- grep("^z_", names(cl), value = TRUE)
    both <- intersect(cl$target[cl$cell_line == "K562"], cl$target[cl$cell_line == "HepG2"])
    zk <- as.matrix(cl[cl$cell_line == "K562", zc][match(both, cl$target[cl$cell_line == "K562"]), ])
    zh <- as.matrix(cl[cl$cell_line == "HepG2", zc][match(both, cl$target[cl$cell_line == "HepG2"]), ])
    dd <- sqrt(outer(rowSums(zk^2), rowSums(zh^2), "+") - 2 * zk %*% t(zh))
    dist_same <- mean(diag(dd)); dist_diff <- mean(dd[row(dd) != col(dd)]); n_diff_pairs <- sum(row(dd) != col(dd))
    corr_diff <- median(cor(t(zk), t(zh))[row(dd) != col(dd)])
  })
}

aim1b_numbers <- function() {
  pairs <- bind_rows(tab("cross_cell_pairs"), tab("cross_cell_null_pairs"))
  dec <- tab("cross_cell_decomposition"); grp <- tab("cross_cell_groups"); gfe <- tab("cross_cell_group_features")
  med <- pairs |> group_by(type) |>
    summarise(site = median(site_matched, na.rm = TRUE), gene = median(gene_frac, na.rm = TRUE),
              map = median(map_r, na.rm = TRUE), .groups = "drop")
  dfrac <- dec |> group_by(type, cls) |> summarise(n = sum(n), .groups = "drop") |> group_by(type) |> mutate(f = n / sum(n))
  list(m = function(t, v) med[[v]][med$type == t],
       unexplained = function(t) dfrac$f[dfrac$type == t & dfrac$cls == "unexplained"],
       frac = function(t, c) dfrac$f[dfrac$type == t & dfrac$cls == c],
       map_r_of = function(rbp) pairs$map_r[pairs$type == "test" & pairs$target == rbp],
       n_map = sum(!is.na(pairs$map_r[pairs$type == "test"])),
       groups = table(grp$group),
       n_pairs = table(pairs$type),
       n_peaks_decomp = tapply(dec$n, dec$type, sum),
       depth_rho = gfe$spearman_with_site_sharing[gfe$feature == "log10_min_peaks"])
}
