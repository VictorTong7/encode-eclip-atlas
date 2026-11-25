# Binding vs expression, figures a and d (COMPUTE). Figures b and c come from 18_change_rate.R; all four use the same
# genes and expression deciles (lib_bvu.R prep_kd_genes()).
#
# Setup: per analysed knockdown (kd_groups.tsv, analysed == TRUE; 115), the gene universe (tested by DESeq2, expressed at
# control nuclear TPM >= 1, not the RBP's own gene) minus genes with missing padj; bound = >= 1 peak of the RBP x cell
# line eCLIP experiment; 10 deciles of log1p control nuclear TPM within the knockdown.
#
#   a  Baseline expression: median control nuclear TPM of bound vs unbound genes (ratio, Wilcoxon p).
#   d  Magnitude: do bound genes change more than unbound genes of the same expression? Change = |shrunken log2 FC|
#      (both directions), all genes, changed or not.
#      Statistic: within each decile, median |log2 FC| bound - unbound; the knockdown's value is the average of the decile
#      differences weighted by the decile's bound-gene count (deciles need >= MIN_PER_GROUP bound and unbound genes).
#      Bias correction: medians of small, skewed groups give a positive value under no effect, so the value expected
#      under no effect (mean of B permutations of the bound labels within deciles) is subtracted from the effect and CI.
#      95% CI: bootstrap of genes within deciles (B); p: two-sided permutation p (B); BH across knockdowns.
#      Gene-level resampling is optimistic (genes in a knockdown are not independent): across-knockdown summaries
#      (median with bootstrap CI, Wilcoxon signed-rank vs 0) are the main evidence.
#      Also in the table: unadjusted difference of medians, uncorrected decile statistic, decile-weighted difference of means.
#      Null check (summary table): random 'bound' sets of the same size.
#
# Outputs: analysis/results/tables/bvu_baseline_expression.tsv, bvu_magnitude.tsv, bvu_magnitude_summary.tsv
#          data/analysis/bound_vs_unbound_plotdata.rds
# Usage  : Rscript analysis/R/13_bound_vs_unbound.R   (~7 min; needs 12_kd_qc.R; run after 18_change_rate.R for the
#          gene-count check)

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/R/lib_bvu.R")
N_DEC <- 10; MIN_PER_GROUP <- 3; B <- 1000; N_BOOT_KD <- 5000
POS <- c("SLBP K562", "CPSF6 K562", "IGF2BP3 HepG2", "DDX6 HepG2", "RBFOX2 K562", "RBFOX2 HepG2")   # positive controls

kd <- read.delim("analysis/results/tables/kd_groups.tsv") |> filter(analysed)
message("Loading gene universes for ", nrow(kd), " knockdowns ...")
kdl <- lapply(load_kd_data(kd), function(x) { x$d <- prep_kd_genes(x$d, N_DEC); x })
info <- function(x) x$info |> dplyr::select(kd_experiment, target, cell_line, group, tier, flag)

# ---- a: baseline expression ----------------------------------------------------------------------------------------
base <- bind_rows(lapply(kdl, function(x) { d <- x$d
  cbind(info(x), data.frame(n_genes = nrow(d), n_bound = sum(d$bound), n_unbound = sum(!d$bound),
                            tpm_bound = median(d$ntpm[d$bound]), tpm_unbound = median(d$ntpm[!d$bound]),
                            wilcox_p = wilcox.test(d$ntpm[d$bound], d$ntpm[!d$bound], exact = FALSE)$p.value))
})) |> mutate(tpm_ratio = tpm_bound / tpm_unbound)

# ---- d: magnitude ------------------------------------------------------------------------------------------------------
# decile-weighted difference of medians (or means); v = |LFC|, b = bound flag, dec = decile
mag_stat <- function(v, b, dec, fun = median) {
  diffs <- vapply(split(seq_along(v), dec), function(i) {
    bi <- i[b[i]]; ui <- i[!b[i]]
    if (length(bi) < MIN_PER_GROUP || length(ui) < MIN_PER_GROUP) return(c(NA, 0))
    c(fun(v[bi]) - fun(v[ui]), length(bi))
  }, numeric(2))
  ok <- !is.na(diffs[1, ])
  if (!any(ok)) return(c(est = NA, n_dec = 0))
  c(est = sum(diffs[1, ok] * diffs[2, ok]) / sum(diffs[2, ok]), n_dec = sum(ok))
}
mag_one <- function(d, seed = 1) {
  v <- d$abs_lfc; b <- d$bound; dec <- d$dec
  est <- mag_stat(v, b, dec)
  if (is.na(est[["est"]])) return(NULL)
  set.seed(seed)
  # usable deciles, split once into bound / unbound values; resampling then works on these short vectors
  grp <- lapply(split(seq_along(v), dec), function(i) list(b = v[i[b[i]]], u = v[i[!b[i]]]))
  grp <- grp[vapply(grp, function(g) length(g$b) >= MIN_PER_GROUP && length(g$u) >= MIN_PER_GROUP, TRUE)]
  w <- vapply(grp, function(g) length(g$b), 0L)
  perm <- replicate(B, sum(w * vapply(grp, function(g) { pool <- c(g$b, g$u); s <- sample.int(length(pool), length(g$b))
    median(pool[s]) - median(pool[-s]) }, 0)) / sum(w))
  boot <- replicate(B, sum(w * vapply(grp, function(g)
    median(g$b[sample.int(length(g$b), replace = TRUE)]) - median(g$u[sample.int(length(g$u), replace = TRUE)]), 0)) / sum(w))
  p <- min(1, 2 * min(mean(c(perm, est[["est"]]) >= est[["est"]]), mean(c(perm, est[["est"]]) <= est[["est"]])))
  pm <- mean(perm)                                                    # value expected under no effect
  data.frame(n_genes = length(v), n_bound = sum(b), n_unbound = sum(!b), n_deciles_used = est[["n_dec"]],
             med_bound = median(v[b]), med_unbound = median(v[!b]), diff_raw = median(v[b]) - median(v[!b]),
             diff_uncorrected = est[["est"]], perm_mean = pm, diff_median = est[["est"]] - pm,
             ci_lo = quantile(boot, 0.025, names = FALSE) - pm, ci_hi = quantile(boot, 0.975, names = FALSE) - pm,
             p = p, diff_mean = mag_stat(v, b, dec, mean)[["est"]])
}
run_mag <- function(kdl, null = FALSE) {
  bind_rows(lapply(seq_along(kdl), function(k) {
    x <- kdl[[k]]; d <- x$d
    if (null) { set.seed(1000 + k); d$bound <- seq_len(nrow(d)) %in% sample.int(nrow(d), sum(d$bound)) }
    r <- mag_one(d, seed = k); if (is.null(r)) return(NULL)
    cbind(info(x), r)
  })) |> mutate(q = p.adjust(p, "BH"))
}
message("Magnitude, and the random-set null (", B, " permutations / bootstraps each) ...")
t0 <- Sys.time()
mag <- run_mag(kdl)
null_mag <- run_mag(kdl, null = TRUE)
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

set.seed(13)
boot_med <- function(v) quantile(replicate(N_BOOT_KD, median(sample(v, replace = TRUE))), c(0.025, 0.975), names = FALSE)
summ_one <- function(r, label) {
  ci <- boot_med(r$diff_median)
  data.frame(set = label, n_kd = nrow(r), median_diff = median(r$diff_median), ci_lo = ci[1], ci_hi = ci[2],
             median_diff_raw = median(r$diff_raw), median_diff_mean = median(r$diff_mean, na.rm = TRUE),
             signed_rank_p = wilcox.test(r$diff_median, exact = FALSE)$p.value, n_gt0 = sum(r$diff_median > 0),
             n_sig_up = sum(r$q < 0.05 & r$diff_median > 0), n_sig_down = sum(r$q < 0.05 & r$diff_median < 0))
}
summ <- bind_rows(summ_one(mag, "main"),
                  summ_one(mag[mag$cell_line == "K562", ], "K562"), summ_one(mag[mag$cell_line == "HepG2", ], "HepG2"),
                  summ_one(mag[mag$flag %in% c("", NA), ], "without FXR2 K562"),
                  summ_one(null_mag, "null (random bound sets)"))

# ---- write ---------------------------------------------------------------------------------------------------------------
r5 <- function(x) signif(x, 5)
w <- function(df, f) write.table(df |> mutate(across(where(is.double), r5)), file.path("analysis/results/tables", f),
                                 sep = "\t", quote = FALSE, row.names = FALSE, na = "")
w(base, "bvu_baseline_expression.tsv"); w(mag, "bvu_magnitude.tsv"); w(summ, "bvu_magnitude_summary.tsv")
saveRDS(list(base = base, mag = mag, summ = summ, params = list(N_DEC = N_DEC, B = B, POS = POS)),
        "data/analysis/bound_vs_unbound_plotdata.rds")

# ---- checks and summary ------------------------------------------------------------------------------------------------------
cr <- "analysis/results/tables/change_rate.tsv"
if (file.exists(cr)) {
  crm <- read.delim(cr) |> filter(version == "main") |> dplyr::select(kd_experiment, n_genes_cr = n_genes, or_mh)
  chk <- base |> inner_join(crm, by = "kd_experiment") |> inner_join(mag |> dplyr::select(kd_experiment, diff_median), by = "kd_experiment")
  stopifnot(nrow(chk) == nrow(base), all(chk$n_genes == chk$n_genes_cr))
  cat(sprintf("Check: same genes as the change-rate metric in all %d knockdowns. Magnitude vs log odds ratio: Spearman %.2f\n",
              nrow(chk), cor(chk$diff_median, log(chk$or_mh), method = "spearman")))
}
cat(sprintf("\na. Bound genes more expressed in %d of %d knockdowns; median TPM %.1f vs %.1f\n",
            sum(base$tpm_ratio > 1), nrow(base), median(base$tpm_bound), median(base$tpm_unbound)))
cat("\nd. Magnitude across knockdowns:\n"); print(summ, digits = 3, row.names = FALSE)
cat("\nPositive controls:\n")
print(mag |> filter(paste(target, cell_line) %in% POS) |>
        dplyr::select(target, cell_line, n_bound, diff_raw, diff_median, ci_lo, ci_hi, q) |> as.data.frame(), digits = 3)
message("Wrote bvu_baseline_expression.tsv, bvu_magnitude.tsv, bvu_magnitude_summary.tsv, bound_vs_unbound_plotdata.rds")
