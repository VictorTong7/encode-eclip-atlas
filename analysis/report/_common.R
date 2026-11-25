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

# ---- knockdown data (analysis/R/12_kd_qc.R) ---------------------------------------------------------------------
kd_numbers <- function() {
  s <- tab("kd_sample_qc") |> filter(shrna); g <- tab("kd_groups"); sel <- tab("kd_selected")
  a <- g |> filter(analysed)
  rho_de <- sapply(split(g, g$cell_line), function(x) cor(x$kd_efficiency, x$n_de, method = "spearman"))
  odd <- s |> filter(qc_flag != "")
  list(n_samples = nrow(s), n_kd_rep = sum(s$kind == "knockdown"), n_ctl_rep = sum(s$kind == "control"),
       reads_min = min(s$counts_M), reads_med = median(s$counts_M), reads_max = max(s$counts_M),
       genes10_min = min(s$genes_c10), genes10_max = max(s$genes_c10), rep_r_min = min(s$rep_r_min),
       chrm_min = min(s$chrM_share), chrm_max = max(s$chrM_share), chrm_med = median(s$chrM_share),
       n_odd_ctl = n_distinct(odd$control_batch), n_odd_kd = n_distinct(odd$experiment[odd$kind == "knockdown"]),
       n_shrna = sum(sel$assay_title == "shRNA RNA-seq"), n_crispri = sum(sel$assay_title != "shRNA RNA-seq"),
       n_dup = sum(sel$reason == "weaker duplicate shRNA knockdown"),
       n_sel = nrow(g), n_rbp = n_distinct(g$target),
       n_rbp_both = sum(table(unique(g[, c("target", "cell_line")])$target) == 2),
       n_ctl = n_distinct(g$control_experiment), max_per_ctl = max(table(g$control_experiment)),
       eff_med = median(g$kd_efficiency), eff_min = min(g$kd_efficiency), eff_max = max(g$kd_efficiency),
       n_weak = sum(g$weak_kd), n_weak_k = sum(g$weak_kd & g$cell_line == "K562"), n_weak_h = sum(g$weak_kd & g$cell_line == "HepG2"),
       eff_ntpm_rho = cor(g$kd_efficiency, g$kd_efficiency_ntpm, method = "spearman"),
       de_med = median(g$n_de), de_min = min(g$n_de), de_max = max(g$n_de), rho_de = rho_de,
       n_binding = sum(!is.na(g$eclip_experiment)), n_low_binding = sum(is.na(g$eclip_experiment)),
       n_ncrna = sum(g$excluded_reason == "ncRNA binder"), n_weak_excl = sum(g$excluded_reason == "weak knockdown"),
       n_analysed = nrow(a), n_an_k = sum(a$cell_line == "K562"), n_an_h = sum(a$cell_line == "HepG2"),
       tiers = table(factor(a$tier, c("< 300", "300–999", "≥ 1,000"))), bound_med = median(a$n_bound),
       reassigned = paste(g$target[g$reassigned %in% TRUE], g$cell_line[g$reassigned %in% TRUE], collapse = " and "),
       fxr2_eff = g$kd_efficiency[g$target == "FXR2" & g$cell_line == "K562"])
}

# ---- binding and RNA change (analysis/R/13, 17, 18, 19) --------------------------------------------------------------
bvu_numbers <- function() {
  base <- tab("bvu_baseline_expression"); cr <- tab("change_rate"); crs <- tab("change_rate_summary")
  ms <- tab("bvu_magnitude_summary"); mg <- tab("bvu_magnitude")
  hs <- tab("hub_attribution_summary"); ha <- tab("hub_attribution"); hg <- tab("hub_groups")
  ur <- tab("hub_unbound_rates") |> filter(n_kd >= 20); go <- tab("hub_go"); pr <- tab("hub_gene_profile")
  m <- cr |> filter(version == "main")
  row <- function(d, s) d[d$set == s, ]
  kd_val <- function(d, col, lab) setNames(d[[col]][match(lab, paste(d$target, d$cell_line))], lab)
  pos <- c("SLBP K562", "CPSF6 K562", "IGF2BP3 HepG2", "RBFOX2 K562", "RBFOX2 HepG2", "DDX6 HepG2")
  urw <- ur |> dplyr::select(cell_line, dec, hub_group, median_rate) |> tidyr::pivot_wider(names_from = hub_group, values_from = median_rate)
  urm <- ur |> group_by(cell_line, hub_group) |> summarise(r = median(median_rate), .groups = "drop")
  gom <- go |> filter(background == "matched")
  list(n_kd = nrow(base), n_more_expr = sum(base$tpm_ratio > 1), tpm_b = median(base$tpm_bound), tpm_u = median(base$tpm_unbound),
       tpm_ratio = median(base$tpm_ratio), bound_med = median(base$n_bound), genes_med = median(base$n_genes),
       cr = row(crs, "main"), cr_hub = row(crs, "hub"), cr_strict = row(crs, "strict"), cr_null = row(crs, "null"),
       cr_k = row(crs, "main K562"), cr_h = row(crs, "main HepG2"), n_rate_higher = sum(m$rate_bound > m$rate_unbound),
       or_pos = kd_val(m, "or_mh", pos), or_top = m |> filter(q_mh < 0.05) |> arrange(desc(or_mh)) |> head(5),
       or_low = m |> filter(q_mh < 0.05, or_mh < 1),
       mag = row(ms, "main"), mag_null = row(ms, "null (random bound sets)"), mag_pos = kd_val(mg, "diff_median", pos),
       mag_rho_or = cor(mg$diff_median, log(m$or_mh[match(mg$kd_experiment, m$kd_experiment)]), method = "spearman"),
       hub_n = table(hg$cell_line, hg$hub_group), hub_tpm = tapply(hg$ntpm, list(hg$cell_line, hg$hub_group), median),
       ur_more = sum(urw$`≥ 20` > urw$`< 5`, na.rm = TRUE), ur_cmp = sum(!is.na(urw$`≥ 20`) & !is.na(urw$`< 5`)), urm = urm,
       hub = row(hs, "main"), hub_null = row(hs, "null (random bound sets)"), hub_share_pos = kd_val(ha, "hub_share", pos),
       hub_cover = median(ha$n_genes / base$n_genes[match(ha$kd_experiment, base$kd_experiment)]),
       go_both = nrow(gom |> count(ontology, ID) |> filter(n == 2)),
       go_top = gom |> filter(ontology == "BP") |> group_by(cell_line) |> slice_min(p.adjust, n = 4, with_ties = FALSE) |> ungroup(),
       prof = pr)
}
