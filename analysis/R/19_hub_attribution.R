# Binding vs expression, hub block (figures e-g; COMPUTE). How many changing targets are explained by being a hub gene
# (bound by many RBPs) rather than by binding of the knocked-down RBP?
#
# Same setup as figures a-d (lib_bvu.R prep_kd_genes(): tested, expressed, not the RBP's own gene, padj present;
# changed = padj < 0.05; 10 expression deciles within the knockdown), restricted to genes with an RBP count in
# gene_rbp_counts.tsv (label_hub_categories(); ~92% of genes, the rest mostly small RNAs).
#   hub group  number of OTHER profiled RBPs binding the gene (n_rbps minus the knocked-down RBP's own binding):
#              < 5 / 5-19 / >= 20.
#   e  Hub groups per cell line (gene-level: all profiled RBPs), genes tested in >= 50% of the cell line's knockdowns:
#      counts and control expression. Also the gene list for 17_hub_go.R.
#   f  Change rate of UNBOUND genes by hub group within expression deciles, per knockdown (hub effect without binding).
#   g  Attribution, per knockdown:
#        O       bound genes that change
#        E_expr  expected from expression only: sum over bound genes of the change rate of unbound genes in the same decile
#        E_hub   expected from expression + hub status: same, using unbound genes of the same decile AND hub group
#                (decile rate if that cell has < MIN_CELL unbound genes; counted)
#        extra changing targets O - E_expr = explained by hub status (E_hub - E_expr) + attributable to binding (O - E_hub)
#        hub share = (E_hub - E_expr) / (O - E_expr), where O - E_expr > 0
#      95% CIs: bootstrap of genes within deciles (B; optimistic). Across knockdowns: totals with bootstrap CIs over
#      knockdowns, median hub share, Wilcoxon signed-rank of each excess. Null: random 'bound' sets of the same size.
#
# Outputs: analysis/results/tables/hub_groups.tsv (gene x cell line; read by 17_hub_go.R), hub_attribution.tsv, hub_unbound_rates.tsv,
#          hub_attribution_summary.tsv; data/analysis/hub_attribution_plotdata.rds
# Usage  : Rscript analysis/R/19_hub_attribution.R   (~5 min; needs 12_kd_qc.R and gene_rbp_counts.tsv from 05)

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/R/lib_bvu.R")
N_DEC <- 10; HUB_BREAKS <- c(-Inf, 4, 19, Inf); GRP_LAB <- c("< 5", "5–19", "≥ 20"); MIN_CELL <- 10
B <- 1000; N_BOOT_KD <- 5000
POS <- c("SLBP K562", "CPSF6 K562", "IGF2BP3 HepG2", "DDX6 HepG2", "RBFOX2 K562", "RBFOX2 HepG2")   # positive controls

kd <- read.delim("analysis/results/tables/kd_groups.tsv") |> filter(analysed)
message("Loading gene universes for ", nrow(kd), " knockdowns ...")
kdl <- lapply(load_kd_data(kd), function(x) { x$d <- prep_kd_genes(x$d, N_DEC); x })
n_setup <- vapply(kdl, function(x) nrow(x$d), 0L)                    # genes in the figure a-d setup, per knockdown
kdl <- label_hub_categories(kdl, 20)                                  # keeps genes with an RBP count; adds n_other
kdl <- lapply(kdl, function(x) { x$d$hub <- as.integer(cut(x$d$n_other, HUB_BREAKS)); x })   # 1, 2, 3
info <- function(x) x$info |> dplyr::select(kd_experiment, target, cell_line, group, tier, flag)

# ---- e: hub groups (gene level: all profiled RBPs, as in gene_rbp_counts.tsv) ----------------------------------------
rbp <- read.delim("analysis/results/tables/gene_rbp_counts.tsv") |>
  dplyr::select(cell_line, gene_id, gene_name, gene_span, n_rbps, frac_rbps)
tested <- bind_rows(lapply(kdl, function(x) data.frame(cell_line = x$info$cell_line, gene_id = x$d$gene_id, ntpm = x$d$ntpm))) |>
  group_by(cell_line, gene_id) |> summarise(n_tested = n(), ntpm = first(ntpm), .groups = "drop") |>
  inner_join(kd |> count(cell_line, name = "n_kd"), by = "cell_line") |> filter(n_tested >= 0.5 * n_kd)
genes <- tested |> inner_join(rbp, by = c("cell_line", "gene_id")) |>
  mutate(hub_group = factor(GRP_LAB[as.integer(cut(n_rbps, HUB_BREAKS))], GRP_LAB))
grp <- genes |> group_by(cell_line, hub_group) |>
  summarise(n = n(), median_tpm = median(ntpm), q25_tpm = quantile(ntpm, 0.25), q75_tpm = quantile(ntpm, 0.75),
            median_span = median(gene_span), .groups = "drop")

# ---- f: change rate of unbound genes by hub group within deciles ----------------------------------------------------
rates <- bind_rows(lapply(kdl, function(x) {
  x$d |> filter(!bound) |> group_by(dec, hub) |> summarise(n = n(), rate = mean(changed), .groups = "drop") |>
    mutate(kd_experiment = x$info$kd_experiment, cell_line = x$info$cell_line)
})) |> mutate(hub_group = factor(GRP_LAB[hub], GRP_LAB))

# ---- g: attribution --------------------------------------------------------------------------------------------------
# expected counts for bound genes from unbound-gene change rates; codes: decile (1-10) and decile x hub (1-30)
attrib_stat <- function(dec, hub, bound, chg) {
  cell <- (dec - 1L) * 3L + hub
  u <- !bound
  r_dec <- tabulate(dec[u & chg], N_DEC) / tabulate(dec[u], N_DEC)
  n_cell <- tabulate(cell[u], 3L * N_DEC)
  r_cell <- tabulate(cell[u & chg], 3L * N_DEC) / n_cell
  use_cell <- n_cell[cell[bound]] >= MIN_CELL
  e_hub_g <- ifelse(use_cell, r_cell[cell[bound]], r_dec[dec[bound]])
  c(O = sum(chg[bound]), E_expr = sum(r_dec[dec[bound]], na.rm = TRUE), E_hub = sum(e_hub_g, na.rm = TRUE),
    n_fallback = sum(!use_cell))
}
attrib_one <- function(d, seed = 1) {
  s <- attrib_stat(d$dec, d$hub, d$bound, d$changed)
  set.seed(seed)
  idx <- split(seq_len(nrow(d)), paste(d$dec, d$bound))                 # resample within decile x bound/unbound
  boot <- replicate(B, { j <- unlist(lapply(idx, function(i) i[sample.int(length(i), replace = TRUE)]))
    b <- attrib_stat(d$dec[j], d$hub[j], d$bound[j], d$changed[j])
    c(hub = b[["E_hub"]] - b[["E_expr"]], binding = b[["O"]] - b[["E_hub"]], total = b[["O"]] - b[["E_expr"]]) })
  ci <- apply(boot, 1, quantile, c(0.025, 0.975), names = FALSE)
  data.frame(n_genes = nrow(d), n_bound = sum(d$bound), n_bound_hub = sum(d$bound & d$hub == 3L),
             O = s[["O"]], E_expr = s[["E_expr"]], E_hub = s[["E_hub"]], n_fallback = s[["n_fallback"]],
             extra_total = s[["O"]] - s[["E_expr"]], extra_hub = s[["E_hub"]] - s[["E_expr"]], extra_binding = s[["O"]] - s[["E_hub"]],
             total_lo = ci[1, "total"], total_hi = ci[2, "total"], hub_lo = ci[1, "hub"], hub_hi = ci[2, "hub"],
             binding_lo = ci[1, "binding"], binding_hi = ci[2, "binding"])
}
run_attrib <- function(null = FALSE) {
  bind_rows(lapply(seq_along(kdl), function(k) {
    x <- kdl[[k]]; d <- x$d
    if (null) { set.seed(2000 + k); d$bound <- seq_len(nrow(d)) %in% sample.int(nrow(d), sum(d$bound)) }
    cbind(info(x), attrib_one(d, seed = k))
  })) |> mutate(hub_share = ifelse(extra_total > 0, extra_hub / extra_total, NA))
}
message("Attribution and random-set null (", B, " bootstraps per knockdown) ...")
t0 <- Sys.time()
att <- run_attrib()
null_att <- run_attrib(null = TRUE)
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

set.seed(19)
boot_kd <- function(r, f) quantile(replicate(N_BOOT_KD, f(r[sample.int(nrow(r), replace = TRUE), ])), c(0.025, 0.975), names = FALSE)
summ_one <- function(r, label) {
  tot <- function(z) sum(z$extra_hub) / sum(z$extra_total)
  ci_share <- boot_kd(r, tot); ci_tot <- boot_kd(r, function(z) sum(z$extra_total))
  data.frame(set = label, n_kd = nrow(r), O = sum(r$O), E_expr = sum(r$E_expr), E_hub = sum(r$E_hub),
             extra_total = sum(r$extra_total), extra_total_lo = ci_tot[1], extra_total_hi = ci_tot[2],
             extra_hub = sum(r$extra_hub), extra_binding = sum(r$extra_binding),
             hub_share_total = tot(r), hub_share_lo = ci_share[1], hub_share_hi = ci_share[2],
             median_hub_share = median(r$hub_share, na.rm = TRUE), n_extra_pos = sum(r$extra_total > 0),
             n_hub_majority = sum(r$hub_share > 0.5, na.rm = TRUE),
             p_hub = wilcox.test(r$extra_hub, exact = FALSE)$p.value, p_binding = wilcox.test(r$extra_binding, exact = FALSE)$p.value,
             n_binding_ci_pos = sum(r$binding_lo > 0), n_hub_ci_pos = sum(r$hub_lo > 0))
}
summ <- bind_rows(summ_one(att, "main"), summ_one(att[att$cell_line == "K562", ], "K562"),
                  summ_one(att[att$cell_line == "HepG2", ], "HepG2"),
                  summ_one(att[att$flag %in% c("", NA), ], "without FXR2 K562"),
                  summ_one(null_att, "null (random bound sets)"))

# ---- write ------------------------------------------------------------------------------------------------------------
r5 <- function(x) signif(x, 5)
w <- function(df, f) write.table(df |> mutate(across(where(is.double), r5)), file.path("analysis/results/tables", f),
                                 sep = "\t", quote = FALSE, row.names = FALSE, na = "")
w(genes |> dplyr::select(cell_line, gene_id, gene_name, ntpm, gene_span, n_rbps, frac_rbps, hub_group), "hub_groups.tsv")
w(att, "hub_attribution.tsv"); w(summ, "hub_attribution_summary.tsv")
# figure f summary for the website: median over knockdowns of the unbound change rate per cell line x decile x hub group
w(rates |> filter(n >= MIN_CELL) |> group_by(cell_line, dec, hub_group) |>
    summarise(median_rate = median(rate), n_kd = n(), .groups = "drop"), "hub_unbound_rates.tsv")
saveRDS(list(genes = genes, grp = grp, rates = rates, att = att, summ = summ,
             params = list(GRP_LAB = GRP_LAB, MIN_CELL = MIN_CELL, B = B, POS = POS)),
        "data/analysis/hub_attribution_plotdata.rds")

# ---- checks and summary -----------------------------------------------------------------------------------------------
cat(sprintf("\nGenes with an RBP count: median %.1f%% of each knockdown's genes; fallback (sparse hub cells) for %d bound genes in total\n",
            100 * median(att$n_genes / n_setup[att$kd_experiment]), sum(att$n_fallback)))
cat("\ne. Hub groups:\n"); print(as.data.frame(grp), digits = 3)
cat("\nf. Unbound change rate by hub group (median over knockdowns, all deciles pooled by mean):\n")
print(rates |> group_by(cell_line, kd_experiment, hub_group) |> summarise(rate = weighted.mean(rate, n), .groups = "drop") |>
        group_by(cell_line, hub_group) |> summarise(median_rate = median(rate), .groups = "drop") |> as.data.frame(), digits = 3)
cat("\ng. Attribution:\n"); print(summ, digits = 3, row.names = FALSE)
cat("\nPositive controls:\n")
print(att |> filter(paste(target, cell_line) %in% POS) |>
        dplyr::select(target, cell_line, n_bound, O, E_expr, E_hub, extra_total, extra_hub, extra_binding, hub_share) |>
        as.data.frame(), digits = 3)
message("Wrote hub_groups.tsv, hub_attribution.tsv, hub_attribution_summary.tsv, hub_attribution_plotdata.rds")
