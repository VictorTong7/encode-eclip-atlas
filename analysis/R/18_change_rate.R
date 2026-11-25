# Binding vs expression, figures b, c and e: conservative change-rate metric (COMPUTE). Uses the same genes and expression
# deciles as 13_bound_vs_unbound.R (figures a, d) via lib_bvu.R prep_kd_genes().
# Question: after an RBP is knocked down, do the genes it binds change at a higher rate than genes it does not bind?
# Only WHETHER a gene changes counts (not how much, not which direction).
#
# Per analysed knockdown (kd_groups.tsv, analysed == TRUE; 115):
#   genes    lib_bvu.R load_kd_data(): tested by DESeq2, expressed (control nuclear TPM >= 1), not the RBP's own gene;
#            genes with missing padj (DESeq2 independent filtering) are excluded (status unknown).
#   bound    >= 1 peak of the RBP x cell line eCLIP experiment.
#   changed  padj < 0.05 (either direction, no fold-change threshold); unchanged (control) = padj >= 0.05
#            (decision 2026-11-24).
#   2 x 2    a = bound changed, b = bound unchanged, c = unbound changed, d = unbound unchanged.
#            odds ratio OR = (a/b)/(c/d); rate ratio RR = (a/(a+b)) / (c/(c+d)).
#   main     adjusted for expression: 2 x 2 within each of N_DEC deciles of log1p nuclear TPM (within the knockdown),
#            combined with Mantel-Haenszel OR (95% CI, CMH chi-square test; mantelhaen.test) and MH rate ratio.
#            Unadjusted OR / Fisher test in the table. BH across knockdowns per version.
# Versions: main (expression deciles); hub (expression decile x hub group < 5 / 5-19 / >= 20 other RBPs; genes without a
#           hub count dropped); strict (unchanged = padj > 0.5 only; 0.05 <= padj <= 0.5 excluded); null (random
#           'bound' sets of the same size, expression deciles; OR should be ~1).
# Across knockdowns: median OR (bootstrap CI), Wilcoxon signed-rank of log OR, counts q < 0.05; by cell line;
# without FXR2 K562.
#
# Outputs: analysis/results/tables/change_rate.tsv, change_rate_summary.tsv; data/analysis/change_rate_plotdata.rds
# Figures: b (rates) and c (main OR) and e (hub / strict checks), drawn by 13_bound_vs_unbound_figures.R.
# Usage  : Rscript analysis/R/18_change_rate.R   (~2 min; needs 12_kd_qc.R and gene_rbp_counts.tsv from 05)

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/R/lib_bvu.R")
N_DEC <- 10; HUB_BREAKS <- c(-Inf, 4, 19, Inf); N_BOOT <- 5000
POS <- c("SLBP K562", "CPSF6 K562", "IGF2BP3 HepG2", "DDX6 HepG2", "RBFOX2 K562", "RBFOX2 HepG2")

kd <- read.delim("analysis/results/tables/kd_groups.tsv") |> filter(analysed)
message("Loading gene universes for ", nrow(kd), " knockdowns ...")
kdl_all <- load_kd_data(kd)                                 # full universe: main / strict / null versions
kdl_hub <- label_hub_categories(kdl_all, 20)               # adds n_other (hub count excl. this RBP); drops genes without it

# one knockdown, one version: x has columns bound, chg (changed: TRUE / FALSE / NA = excluded) and stratum
cr_one <- function(x) {
  x <- x[!is.na(x$chg), ]
  tab <- function(z) c(a = sum(z$bound & z$chg), b = sum(z$bound & !z$chg), c = sum(!z$bound & z$chg), d = sum(!z$bound & !z$chg))
  t0 <- tab(x)
  st <- split(x, x$stratum)
  st <- st[vapply(st, function(z) any(z$bound) && any(!z$bound) && nrow(z) >= 2, TRUE)]   # informative strata only
  arr <- vapply(st, function(z) { t <- tab(z); matrix(c(t[["a"]], t[["c"]], t[["b"]], t[["d"]]), 2) }, matrix(0, 2, 2))
  mh_num <- sum(vapply(st, function(z) { t <- tab(z); t[["a"]] * t[["d"]] / sum(t) }, 0))
  mh_den <- sum(vapply(st, function(z) { t <- tab(z); t[["b"]] * t[["c"]] / sum(t) }, 0))
  rr_num <- sum(vapply(st, function(z) { t <- tab(z); t[["a"]] * (t[["c"]] + t[["d"]]) / sum(t) }, 0))
  rr_den <- sum(vapply(st, function(z) { t <- tab(z); t[["c"]] * (t[["a"]] + t[["b"]]) / sum(t) }, 0))
  cmh <- tryCatch(mantelhaen.test(arr, exact = FALSE), error = function(e) NULL)
  fis <- fisher.test(matrix(c(t0[["a"]], t0[["c"]], t0[["b"]], t0[["d"]]), 2))
  data.frame(n_genes = nrow(x), a = t0[["a"]], b = t0[["b"]], c = t0[["c"]], d = t0[["d"]],
             rate_bound = t0[["a"]] / (t0[["a"]] + t0[["b"]]), rate_unbound = t0[["c"]] / (t0[["c"]] + t0[["d"]]),
             or_raw = unname(fis$estimate), p_raw = fis$p.value,
             rr_raw = (t0[["a"]] / (t0[["a"]] + t0[["b"]])) / (t0[["c"]] / (t0[["c"]] + t0[["d"]])),
             n_strata = length(st), or_mh = mh_num / mh_den,
             or_lo = if (is.null(cmh)) NA else cmh$conf.int[1], or_hi = if (is.null(cmh)) NA else cmh$conf.int[2],
             p_mh = if (is.null(cmh)) NA else cmh$p.value, rr_mh = rr_num / rr_den)
}

prep <- function(d, version, seed = 1) {
  d <- prep_kd_genes(d, N_DEC)          # shared setup with 13_bound_vs_unbound.R (same genes, same deciles)
  d$chg <- if (version == "strict") ifelse(d$changed, TRUE, ifelse(d$padj > 0.5, FALSE, NA)) else d$changed
  if (version == "hub") d <- d |> mutate(stratum = paste(dec, cut(n_other, HUB_BREAKS)))
  else d <- d |> mutate(stratum = as.character(dec))
  if (version == "null") { set.seed(seed); d$bound <- seq_len(nrow(d)) %in% sample.int(nrow(d), sum(d$bound)) }
  d
}

run_version <- function(version) {
  src <- if (version == "hub") kdl_hub else kdl_all
  bind_rows(lapply(seq_along(src), function(i) {
    x <- src[[i]]
    cbind(x$info |> dplyr::select(kd_experiment, target, cell_line, group, tier, flag), version = version,
          cr_one(prep(x$d, version, seed = i)))
  })) |> mutate(q_mh = p.adjust(p_mh, "BH"), q_raw = p.adjust(p_raw, "BH"))
}

message("Computing change-rate tables (main, hub, strict, null) ...")
res <- bind_rows(lapply(c("main", "hub", "strict", "null"), run_version))

# ---- across knockdowns -------------------------------------------------------------------------------------------
set.seed(18)
boot_med <- function(v) quantile(replicate(N_BOOT, median(sample(v, replace = TRUE))), c(0.025, 0.975), names = FALSE)
summ_one <- function(r, label) {
  v <- log(r$or_mh); v <- v[is.finite(v)]
  ci <- exp(boot_med(v))
  data.frame(set = label, n_kd = nrow(r), median_rate_bound = median(r$rate_bound), median_rate_unbound = median(r$rate_unbound),
             median_or_raw = median(r$or_raw), median_or = exp(median(v)), or_ci_lo = ci[1], or_ci_hi = ci[2],
             median_rr = median(r$rr_mh[is.finite(r$rr_mh)]), signed_rank_p = wilcox.test(v, exact = FALSE)$p.value,
             n_or_gt1 = sum(r$or_mh > 1, na.rm = TRUE), n_sig_up = sum(r$q_mh < 0.05 & r$or_mh > 1, na.rm = TRUE),
             n_sig_down = sum(r$q_mh < 0.05 & r$or_mh < 1, na.rm = TRUE), n_undefined = sum(!is.finite(r$or_mh)))
}
summ <- bind_rows(
  lapply(c("main", "hub", "strict", "null"), function(v) summ_one(res[res$version == v, ], v)),
  lapply(c("K562", "HepG2"), function(cl) summ_one(res[res$version == "main" & res$cell_line == cl, ], paste("main", cl))),
  list(summ_one(res[res$version == "main" & res$flag %in% c("", NA), ], "main without FXR2 K562")))

# ---- write ---------------------------------------------------------------------------------------------------------
r5 <- function(x) signif(x, 5)
write.table(res |> mutate(across(where(is.double), r5)), "analysis/results/tables/change_rate.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(summ |> mutate(across(where(is.double), r5)), "analysis/results/tables/change_rate_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
saveRDS(list(res = res, summ = summ, params = list(N_DEC = N_DEC, POS = POS)), "data/analysis/change_rate_plotdata.rds")

# ---- summary --------------------------------------------------------------------------------------------------------
m <- res |> filter(version == "main")
cat(sprintf("\nChecks: 2x2 totals = genes with padj in all knockdowns: %s; NA-padj genes excluded (median %.1f%%)\n",
            all(m$a + m$b + m$c + m$d == m$n_genes),
            100 * median(1 - m$n_genes / vapply(kdl_all, function(x) nrow(x$d), 0L)[m$kd_experiment])))
cat("\nAcross knockdowns:\n"); print(summ, digits = 3, row.names = FALSE)
cat("\nPositive controls (main):\n")
print(m |> filter(paste(target, cell_line) %in% POS) |>
        dplyr::select(target, cell_line, a, b, c, d, rate_bound, rate_unbound, or_raw, or_mh, or_lo, or_hi, rr_mh, q_mh) |>
        as.data.frame(), digits = 3)
cat("\nStrongest effects (main, q < 0.05):\n")
print(m |> filter(q_mh < 0.05) |> arrange(desc(or_mh)) |> head(12) |>
        dplyr::select(target, cell_line, rate_bound, rate_unbound, or_mh, rr_mh, q_mh) |> as.data.frame(), digits = 3)
cat("\nKnockdowns with OR < 1 (q < 0.05):\n")
print(m |> filter(q_mh < 0.05, or_mh < 1) |> dplyr::select(target, cell_line, rate_bound, rate_unbound, or_mh, q_mh) |> as.data.frame(), digits = 3)
message("Wrote change_rate.tsv, change_rate_summary.tsv, change_rate_plotdata.rds")
