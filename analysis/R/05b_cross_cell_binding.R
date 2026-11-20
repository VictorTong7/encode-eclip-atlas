# Aim 1b: is RBP binding intrinsic to the protein? Compare each RBP's binding between K562 and
# HepG2 at the level of sites, target genes and positional maps; explain cell-specific peaks by
# expression; group RBPs into conserved / intermediate / cell-specific and characterise the groups.
#
# Comparisons (every metric):
#   test     same RBP, K562 vs HepG2 (primary experiment = most peaks per RBP x cell line)
#   ceiling  same RBP, two independent K562 experiments (6 RBPs)
#   null1    different RBPs, K562 vs HepG2
#   null2    different RBPs from the same Aim 1 binding cluster, K562 vs HepG2
#
# Inputs : data/processed/{peaks_annotated,expression_by_cell,experiment_meta}.parquet,
#          data/processed/annot_features.rds, data/ref/gencode_v19.sqlite,
#          analysis/results/tables/{rbp_clusters,region_enrichment}.tsv
# Outputs: analysis/results/tables/cross_cell_*.tsv, data/analysis/cross_cell.rds (for figures)
#
# Usage: Rscript analysis/R/05b_cross_cell_binding.R   (from the project root)

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomicFeatures); library(dplyr); library(tidyr)
})
source("preprocess/R/lib_regions.R")
set.seed(7)
MAXGAP <- 50L          # peaks within this distance count as the same site
MAP_HALF <- 300L; MAP_BIN <- 10L
MIN_MAP_HITS <- 30L

# ---- inputs --------------------------------------------------------------------------------------
cl <- read.delim("analysis/results/tables/rbp_clusters.tsv")
meta_x <- arrow::read_parquet("data/processed/experiment_meta.parquet")
expr <- arrow::read_parquet("data/processed/expression_by_cell.parquet")
feat <- readRDS("data/processed/annot_features.rds")
peaks <- arrow::read_parquet("data/processed/peaks_annotated.parquet",
                             col_select = c("experiment", "chrom", "start", "end", "strand",
                                            "signalValue", "region", "gene_id")) |>
  filter(experiment %in% cl$experiment)

both <- intersect(cl$target[cl$cell_line == "K562"], cl$target[cl$cell_line == "HepG2"])
primary <- cl |> filter(target %in% both) |> group_by(target, cell_line) |>
  slice_max(n_peaks, n = 1, with_ties = FALSE) |> ungroup()
pk <- primary |> select(target, cell_line, experiment) |> pivot_wider(names_from = cell_line, values_from = experiment)
dups <- cl |> filter(cell_line == "K562") |> group_by(target) |> filter(n() == 2) |>
  arrange(target, -n_peaks) |> summarise(A = experiment[1], B = experiment[2], .groups = "drop")
clus <- setNames(cl$cluster, cl$experiment)
message(sprintf("%d RBPs in both lines (>= 100 peaks); %d ceiling pairs", nrow(pk), nrow(dups)))

# ---- per-experiment peak midpoints ----------------------------------------------------------------
cell_of <- setNames(cl$cell_line, cl$experiment)
tpm_tab <- expr |> select(gene_id, cell_line, median_tpm, small_ncRNA) |>
  pivot_wider(names_from = cell_line, values_from = median_tpm)
small_ids <- tpm_tab$gene_id[tpm_tab$small_ncRNA]
peaks <- peaks |> group_by(experiment) |>
  mutate(strength_pct = rank(signalValue) / n()) |> ungroup()
mk_gr <- function(d) {
  mid <- d$start + (d$end - d$start) %/% 2 + 1L
  GRanges(d$chrom, IRanges(mid, width = 1L), strand = d$strand, signal = d$signalValue,
          gene_id = d$gene_id, region = d$region, strength_pct = d$strength_pct)
}
gr <- lapply(split(peaks, peaks$experiment), mk_gr)

universe_genes <- function(thr) {
  u <- tpm_tab$gene_id[(tpm_tab$K562 >= thr & tpm_tab$HepG2 >= thr) %in% TRUE]
  union(u, small_ids)
}

# ---- Level 1 + 2: sites and genes ------------------------------------------------------------------
pair_stats <- function(A, B, ga, gb, maxgap = MAXGAP) {
  if (length(A) == 0 || length(B) == 0) return(c(site_frac = NA, site_matched = NA, site_jacc = NA, site_rho = NA,
                                                  gene_frac = NA, gene_jacc = NA, gene_rho = NA))
  # depth-matched: keep only the strongest n peaks of the larger set (n = size of the smaller set),
  # then average the two directional shared fractions. Removes the bias from unequal peak counts.
  n <- min(length(A), length(B))
  Am <- A[order(-A$signal)[seq_len(n)]]; Bm <- B[order(-B$signal)[seq_len(n)]]
  hm <- findOverlaps(Am, Bm, maxgap = maxgap)
  site_matched <- (length(unique(queryHits(hm))) + length(unique(subjectHits(hm)))) / (2 * n)
  h <- findOverlaps(A, B, maxgap = maxgap)
  a_in <- unique(queryHits(h)); b_in <- unique(subjectHits(h))
  site_frac <- if (length(A) <= length(B)) length(a_in) / length(A) else length(b_in) / length(B)
  shared <- (length(a_in) + length(b_in)) / 2
  site_jacc <- shared / (length(A) + length(B) - shared)
  site_rho <- NA
  if (length(a_in) >= 10) {
    bm <- tapply(B$signal[subjectHits(h)], queryHits(h), max)
    site_rho <- suppressWarnings(cor(A$signal[as.integer(names(bm))], bm, method = "spearman"))
  }
  gi <- intersect(names(ga), names(gb)); gu <- union(names(ga), names(gb))
  gene_frac <- length(gi) / min(length(ga), length(gb))
  gene_jacc <- length(gi) / length(gu)
  va <- ga[gu]; vb <- gb[gu]; va[is.na(va)] <- 0; vb[is.na(vb)] <- 0
  gene_rho <- suppressWarnings(cor(va, vb, method = "spearman"))
  c(site_frac = site_frac, site_matched = site_matched, site_jacc = site_jacc, site_rho = site_rho,
    gene_frac = gene_frac, gene_jacc = gene_jacc, gene_rho = gene_rho)
}

build_pairs <- function() {
  test <- pk |> transmute(type = "test", target, target_b = target, a = K562, b = HepG2)
  ceil <- dups |> transmute(type = "ceiling", target, target_b = target, a = A, b = B)
  nul <- expand.grid(i = seq_len(nrow(pk)), j = seq_len(nrow(pk))) |> filter(i != j) |>
    transmute(type = "null1", target = pk$target[i], target_b = pk$target[j], a = pk$K562[i], b = pk$HepG2[j])
  nul$type[clus[nul$a] == clus[nul$b]] <- "null2"   # null2 is a subset of null1, kept separately
  bind_rows(test, ceil, nul)
}
pairs <- build_pairs()

run_level12 <- function(pairs, thr = 1, maxgap = MAXGAP) {
  ug <- universe_genes(thr)
  ex <- unique(c(pairs$a, pairs$b))
  gu <- lapply(setNames(ex, ex), function(e) { g <- gr[[e]]; g[g$gene_id %in% ug] })
  gs <- lapply(gu, function(g) tapply(g$signal, g$gene_id, sum))
  st <- t(vapply(seq_len(nrow(pairs)), function(i) {
    a <- pairs$a[i]; b <- pairs$b[i]
    pair_stats(gu[[a]], gu[[b]], gs[[a]], gs[[b]], maxgap)
  }, numeric(7)))
  bind_cols(pairs, as.data.frame(st))
}
CACHE <- "data/analysis/cross_cell_cache.rds"
if (file.exists(CACHE) && !("--recompute" %in% commandArgs(TRUE))) {
  message("Using cached overlap results (", CACHE, "); pass --recompute to rebuild")
  cache <- readRDS(CACHE)
} else cache <- list()
if (is.null(cache$lv12)) {
  message("Level 1-2: site and gene overlap for ", nrow(pairs), " pairs ...")
  cache$lv12 <- run_level12(pairs)
}
lv12 <- cache$lv12

# ---- Level 3: positional binding maps ----------------------------------------------------------------
message("Level 3: positional maps ...")
txdb <- AnnotationDbi::loadDb("data/ref/gencode_v19.sqlite")
tx_keep <- unique(feat$tx_id)
intr <- unlist(intronsByTranscript(txdb, use.names = TRUE)[tx_keep])
cds_r <- unlist(range(cdsBy(txdb, "tx", use.names = TRUE)))
cds_r <- cds_r[names(cds_r) %in% tx_keep]
lm_one <- function(x, fix) { y <- unique(granges(resize(x, 1L, fix = fix))); y[seqnames(y) %in% CANONICAL] }
landmarks <- list(`5′ splice site` = lm_one(intr, "start"), `3′ splice site` = lm_one(intr, "end"),
                  `Start codon` = lm_one(cds_r, "start"), `Stop codon` = lm_one(cds_r, "end"))
bins <- seq(-MAP_HALF, MAP_HALF, by = MAP_BIN)
exp_map <- function(e) {
  m <- gr[[e]]
  unlist(lapply(landmarks, function(L) {
    h <- findOverlaps(m, resize(L, 2L * MAP_HALF + 1L, fix = "center"))
    if (length(h) == 0) return(numeric(length(bins) - 1L))
    p <- start(m)[queryHits(h)]; l <- start(L)[subjectHits(h)]
    d <- ifelse(as.character(strand(L))[subjectHits(h)] == "-", l - p, p - l)
    as.numeric(table(cut(d, bins, include.lowest = TRUE))) / length(m)
  }))
}
exps_all <- unique(c(pairs$a, pairs$b))
maps <- sapply(exps_all, exp_map)                       # rows = landmark x bin, cols = experiment
map_hits <- colSums(maps) * sapply(exps_all, function(e) length(gr[[e]]))
map_cor <- function(a, b) {
  if (map_hits[a] < MIN_MAP_HITS || map_hits[b] < MIN_MAP_HITS) return(NA_real_)
  suppressWarnings(cor(maps[, a], maps[, b]))
}
lv12$map_r <- mapply(map_cor, lv12$a, lv12$b)

# ---- Level 4: why are some peaks found in only one line? -----------------------------------------------
message("Level 4: explaining one-line-only peaks ...")
directional <- function(a, b, type, target) {
  A <- gr[[a]]; A <- A[!is.na(A$gene_id)]
  B <- gr[[b]]
  sh <- overlapsAny(A, B, maxgap = MAXGAP)
  ca <- cell_of[[a]]; cb <- cell_of[[b]]
  i <- match(A$gene_id, tpm_tab$gene_id)
  t_this <- tpm_tab[[ca]][i]; t_other <- tpm_tab[[cb]][i]
  small <- A$gene_id %in% small_ids
  lr <- log2((t_other + 0.5) / (t_this + 0.5)); lr[small] <- 0
  lr <- pmin(pmax(lr, -6), 6)
  cls <- case_when(sh ~ "shared",
                   !small & t_other < EXPRESSED_TPM ~ "not expressed in other line",
                   !small & lr <= -1 ~ "expressed ≥2× lower",
                   A$strength_pct < 1 / 3 ~ "weak peak",
                   TRUE ~ "unexplained")
  data.frame(type = type, target = target, from = a, to = b, from_cell = ca, shared = sh, cls = cls,
             log2_tpm_ratio = lr, strength_pct = A$strength_pct, region = A$region)
}
pk4 <- bind_rows(
  lapply(seq_len(nrow(pk)), function(i) bind_rows(directional(pk$K562[i], pk$HepG2[i], "test", pk$target[i]),
                                                  directional(pk$HepG2[i], pk$K562[i], "test", pk$target[i]))),
  lapply(seq_len(nrow(dups)), function(i) bind_rows(directional(dups$A[i], dups$B[i], "ceiling", dups$target[i]),
                                                    directional(dups$B[i], dups$A[i], "ceiling", dups$target[i]))))
CLS <- c("shared", "not expressed in other line", "expressed ≥2× lower", "weak peak", "unexplained")
decomp <- pk4 |> mutate(cls = factor(cls, CLS)) |> count(type, target, cls, .drop = FALSE) |>
  group_by(type, target) |> filter(sum(n) > 0) |> mutate(frac = n / sum(n)) |> ungroup()

# Pooled model of sharing from expression + peak strength + region (test pairs only).
tst <- pk4 |> filter(type == "test") |> mutate(region = factor(region))
fit_full <- glm(shared ~ log2_tpm_ratio + strength_pct + region, binomial, data = tst)
fit_noexp <- glm(shared ~ strength_pct + region, binomial, data = tst)
fit_null <- glm(shared ~ 1, binomial, data = tst)
pr2 <- function(f) 1 - deviance(f) / deviance(fit_null)
auc <- function(y, p) { r <- rank(p); n1 <- as.numeric(sum(y)); n0 <- as.numeric(sum(!y)); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
tst$fitted <- fitted(fit_full)
model_summary <- data.frame(
  model = c("expression + strength + region", "strength + region (no expression)"),
  pseudo_r2 = round(c(pr2(fit_full), pr2(fit_noexp)), 4),
  auc = round(c(auc(tst$shared, fitted(fit_full)), auc(tst$shared, fitted(fit_noexp))), 4))
oe <- tst |> group_by(target) |>
  summarise(observed_shared = mean(shared), expected_shared = mean(fitted), .groups = "drop") |>
  mutate(obs_over_exp = observed_shared / expected_shared)

# ---- grouping ---------------------------------------------------------------------------------------------
classify <- function(l12) {
  ceil_min <- min(l12$site_matched[l12$type == "ceiling"], na.rm = TRUE)
  tt <- l12 |> filter(type == "test")
  bind_rows(lapply(seq_len(nrow(tt)), function(i) {
    x <- tt$target[i]
    own <- l12 |> filter(type %in% c("null1", "null2"), target == x | target_b == x)
    n2 <- own$site_matched[own$type == "null2"]
    nul <- if (sum(!is.na(n2)) >= 5) n2 else own$site_matched
    obs <- tt$site_matched[i]
    data.frame(target = x, site_matched = obs, null_used = if (sum(!is.na(n2)) >= 5) "null2" else "null1",
               null_n = sum(!is.na(nul)), null_median = median(nul, na.rm = TRUE),
               null_q95 = quantile(nul, 0.95, na.rm = TRUE, names = FALSE),
               pct_null_ge = mean(nul >= obs, na.rm = TRUE), ceiling_min = ceil_min,
               group = case_when(obs > quantile(nul, 0.95, na.rm = TRUE) & obs >= ceil_min ~ "conserved",
                                 obs <= median(nul, na.rm = TRUE) ~ "cell-specific",
                                 TRUE ~ "intermediate"))
  }))
}
groups <- classify(lv12)
print(table(groups$group))

message("Robustness: alternative overlap windows and expression thresholds ...")
if (is.null(cache$alt)) {
  cache$alt <- lapply(list(c(25, 1), c(100, 1), c(50, 5)), function(s)
    list(s = s, l12 = run_level12(pairs, thr = s[2], maxgap = as.integer(s[1]))))
  saveRDS(cache, CACHE)
}
robust <- bind_rows(lapply(cache$alt, function(z) {
  s <- z$s; g <- classify(z$l12)
  data.frame(setting = sprintf("maxgap %d nt, TPM >= %d", s[1], s[2]),
             agreement = mean(g$group == groups$group[match(g$target, groups$target)]),
             conserved = sum(g$group == "conserved"), intermediate = sum(g$group == "intermediate"),
             cell_specific = sum(g$group == "cell-specific"))
}))

# ---- features of each RBP ------------------------------------------------------------------------------------
message("Characterising groups ...")
rbp_gene <- feat |> as.data.frame() |> distinct(gene_name, gene_id) |> group_by(gene_name) |> slice(1) |> ungroup()
own_expr <- pk |> select(target) |> left_join(rbp_gene, by = c("target" = "gene_name")) |>
  left_join(tpm_tab |> select(gene_id, K562, HepG2), by = "gene_id") |>
  transmute(target, rbp_log2_tpm = log2((K562 + HepG2) / 2 + 1), rbp_abs_log2_ratio = abs(log2((K562 + 1) / (HepG2 + 1))))
enr <- read.delim("preprocess/results/tables/region_enrichment.tsv")
top_region <- enr |> filter(target %in% pk$target, experiment %in% c(pk$K562, pk$HepG2)) |>
  group_by(target, region) |> summarise(v = mean(log2_vs_all_rbps), .groups = "drop") |>
  group_by(target) |> slice_max(v, n = 1, with_ties = FALSE) |> ungroup() |> select(target, top_region = region)
tgt_div <- bind_rows(lapply(seq_len(nrow(pk)), function(i) {
  g <- unique(c(gr[[pk$K562[i]]]$gene_id, gr[[pk$HepG2[i]]]$gene_id)); g <- g[!is.na(g) & !(g %in% small_ids)]
  j <- match(g, tpm_tab$gene_id)
  data.frame(target = pk$target[i],
             target_expr_divergence = median(abs(log2((tpm_tab$K562[j] + 1) / (tpm_tab$HepG2[j] + 1))), na.rm = TRUE))
}))
go_terms <- c(nucleus = "GO:0005634", cytoplasm = "GO:0005737", rna_splicing = "GO:0008380",
              translation = "GO:0006412", mrna_stability = "GO:0043488", ribosome_biogenesis = "GO:0042254")
go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = pk$target,
                                             keytype = "SYMBOL", columns = "GOALL"))
go_flags <- sapply(go_terms, function(t) pk$target %in% go$SYMBOL[go$GOALL %in% t])
go_flags[!pk$target %in% go$SYMBOL[!is.na(go$GOALL)], ] <- NA
yr <- setNames(meta_x$release_year, meta_x$experiment)
nk <- setNames(cl$n_peaks, cl$experiment); cs <- setNames(cl$coord_source, cl$experiment)
features <- pk |> transmute(
  target,
  log10_min_peaks = log10(pmin(nk[K562], nk[HepG2])),
  abs_log2_peak_ratio = abs(log2(nk[K562] / nk[HepG2])),
  lifted_2022_batch = cs[K562] == "GRCh38_liftover" | cs[HepG2] == "GRCh38_liftover",
  same_release_year = yr[K562] == yr[HepG2],
  frac_peaks_shared_expressed = sapply(seq_len(nrow(pk)), function(i) {
    g <- c(gr[[K562[i]]]$gene_id, gr[[HepG2[i]]]$gene_id); mean(g %in% universe_genes(1)) })) |>
  left_join(own_expr, by = "target") |> left_join(tgt_div, by = "target") |> left_join(top_region, by = "target") |>
  mutate(intronic_binder = top_region %in% c("proximal_intron", "distal_intron"),
         utr3_or_cds_binder = top_region %in% c("3UTR", "CDS")) |>
  bind_cols(as.data.frame(go_flags))

# Effect sizes: conserved vs not conserved (intermediate + cell-specific). The planned
# conserved vs cell-specific contrast is not estimable: almost no RBP is cell-specific.
# Also a cut-off-free check: Spearman correlation of each feature with site sharing.
cliff <- function(x, y) mean(sign(outer(x, y, "-")))
boot_ci <- function(f, x, y, B = 2000) {
  v <- replicate(B, f(sample(x, replace = TRUE), sample(y, replace = TRUE)))
  quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
fg <- features |> left_join(groups |> select(target, group, site_matched), by = "target") |>
  mutate(conserved = group == "conserved")
feat_names <- setdiff(names(features), c("target", "top_region"))
eff <- bind_rows(lapply(feat_names, function(f) {
  v <- fg[[f]]; ok <- !is.na(v)
  x <- v[ok & fg$conserved]; y <- v[ok & !fg$conserved]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  if (is.logical(v)) {
    est <- mean(x) - mean(y); ci <- boot_ci(function(a, b) mean(a) - mean(b), x, y)
    p <- fisher.test(table(c(x, y), rep(c("c", "n"), c(length(x), length(y)))))$p.value
    kind <- "difference in proportion"
  } else {
    est <- cliff(x, y); ci <- boot_ci(cliff, x, y)
    p <- suppressWarnings(wilcox.test(x, y)$p.value); kind <- "Cliff's delta"
  }
  rho <- suppressWarnings(cor.test(as.numeric(v[ok]), fg$site_matched[ok], method = "spearman", exact = FALSE))
  data.frame(feature = f, effect = kind, estimate = est, ci_low = ci[1], ci_high = ci[2], p = p,
             conserved_n = length(x), not_conserved_n = length(y),
             spearman_with_site_sharing = unname(rho$estimate), spearman_p = rho$p.value)
})) |> mutate(p_bh = p.adjust(p, "BH"), spearman_p_bh = p.adjust(spearman_p, "BH")) |> arrange(p)

lm_dat <- fg |> mutate(y = conserved)
fit_grp <- glm(y ~ log10_min_peaks + intronic_binder + target_expr_divergence + nucleus, binomial, data = lm_dat)
grp_model <- data.frame(term = names(coef(fit_grp)), odds_ratio = exp(coef(fit_grp)),
                        ci_low = exp(confint.default(fit_grp)[, 1]), ci_high = exp(confint.default(fit_grp)[, 2]),
                        p = summary(fit_grp)$coefficients[, 4], row.names = NULL)

# ---- write -----------------------------------------------------------------------------------------------------
w <- function(x, f) write.table(x, file.path("analysis/results/tables", f), sep = "\t", quote = FALSE, row.names = FALSE)
lv12_out <- lv12 |> mutate(across(where(is.numeric), ~ round(.x, 4)))
w(lv12_out |> filter(type %in% c("test", "ceiling")), "cross_cell_pairs.tsv")
w(lv12_out |> filter(type %in% c("null1", "null2")), "cross_cell_null_pairs.tsv")
w(decomp |> mutate(frac = round(frac, 4)), "cross_cell_decomposition.tsv")
w(model_summary, "cross_cell_sharing_model.tsv")
w(groups |> left_join(oe, by = "target") |> left_join(lv12 |> filter(type == "test") |>
    select(target, site_frac, gene_frac, gene_rho, site_rho, map_r), by = "target") |>
    mutate(across(where(is.numeric), ~ round(.x, 4))) |> arrange(desc(site_matched)), "cross_cell_groups.tsv")
w(robust, "cross_cell_group_robustness.tsv")
w(fg |> mutate(across(where(is.numeric), ~ round(.x, 4))), "cross_cell_rbp_features.tsv")
w(eff |> mutate(across(where(is.numeric), ~ signif(.x, 3))), "cross_cell_group_features.tsv")
w(grp_model |> mutate(across(where(is.numeric), ~ signif(.x, 3))), "cross_cell_group_model.tsv")
saveRDS(list(lv12 = lv12, decomp = decomp, groups = groups, maps = maps, bins = bins, landmarks = names(landmarks),
             pk = pk, dups = dups, gr_targets = lapply(gr, function(g) data.frame(gene_id = g$gene_id, signal = g$signal)),
             universe = universe_genes(1), eff = eff, oe = oe, fg = fg),
        "data/analysis/cross_cell.rds")

# ---- console summary ---------------------------------------------------------------------------------------------
smry <- lv12 |> group_by(type) |>
  summarise(n = n(), site_frac = median(site_frac, na.rm = TRUE), site_matched = median(site_matched, na.rm = TRUE), gene_frac = median(gene_frac, na.rm = TRUE),
            gene_rho = median(gene_rho, na.rm = TRUE), map_r = median(map_r, na.rm = TRUE), .groups = "drop")
cat("\nMedians by comparison:\n"); print(as.data.frame(smry |> mutate(across(where(is.numeric), ~ round(.x, 3)))))
cat("\nOne-line-only peak explanations (pooled):\n")
print(as.data.frame(decomp |> group_by(type, cls) |> summarise(n = sum(n), .groups = "drop") |>
                      group_by(type) |> mutate(frac = round(n / sum(n), 3))))
cat("\nSharing model:\n"); print(model_summary)
cat("\nGroups:\n"); print(table(groups$group))
cat("\nRobustness:\n"); print(robust)
cat("\nGroup features (top by p):\n"); print(head(eff |> mutate(across(where(is.numeric), ~ signif(.x, 3))), 12))
cat("\nGroup model:\n"); print(grp_model |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
