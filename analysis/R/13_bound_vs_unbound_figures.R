# Binding vs expression: FIGURES (no computation), one sequential set.
#   a  baseline expression of bound vs unbound genes (13_bound_vs_unbound.R)
#   b  change rate, bound vs unbound (18_change_rate.R)          c  expression-adjusted odds ratio of changing (18)
#   d  magnitude: all bound vs unbound genes (13)
#   hub block: e hub groups, f change rate of unbound genes by hub group, g extra changing targets explained by hub
#   status vs binding (19_hub_attribution.R); h what hub genes are, GO (17_hub_go.R)
# a-d share one gene set and one set of expression deciles (lib_bvu.R prep_kd_genes()).
# Shared colours / theme / saving: analysis/R/lib_plot.R. Per-figure sizes: FIG block below.
#
# Usage: Rscript analysis/R/13_bound_vs_unbound_figures.R   (after 13, 17, 18 and 19)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(patchwork) })
source("analysis/R/lib_plot.R")

FIG <- list(
  baseline = list(w = "page", h = 4.4),
  rates    = list(w = "page", h = 4.6),
  odds     = list(w = "page", h = 5.0),
  magnitude = list(w = "page", h = 5.4),
  hubgrp   = list(w = "page", h = 4.2),
  hubrates = list(w = "page", h = 4.4),
  attrib   = list(w = "page", h = 5.4),
  hubgo    = list(w = "page", h = 8.6)
)
cell_f <- function(x) factor(x, names(CELL_COL))
fmt_p <- function(p) ifelse(p < 0.001, sprintf("%.0e", p), sprintf("%.3f", p))

# ---- a: bound genes are more highly expressed -----------------------------------------------------------------------
fig_baseline <- function(d) {
  b <- d$base |> mutate(cell_line = cell_f(cell_line)) |> arrange(tpm_ratio) |> mutate(rank = row_number())
  ggplot(b, aes(rank)) +
    geom_segment(aes(y = tpm_unbound, yend = tpm_bound), colour = GRID, linewidth = 0.7) +
    geom_point(aes(y = tpm_unbound, shape = "Unbound genes"), colour = MUTED, size = 1.5) +
    geom_point(aes(y = tpm_bound, colour = cell_line, shape = "Bound genes"), size = 1.8) +
    scale_colour_manual(values = CELL_COL, name = NULL) +
    scale_shape_manual(values = c(`Bound genes` = 16, `Unbound genes` = 1), name = NULL) +
    scale_y_log10(labels = scales::label_comma()) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2, override.aes = list(colour = INK2))) +
    labs(title = sprintf("Bound genes are more highly expressed than unbound genes in %d of %d knockdowns",
                         sum(b$tpm_ratio > 1), nrow(b)),
         subtitle = sprintf("Median expression in the control samples (nuclear TPM, log scale); median %.0f vs %.0f TPM, %.1f-fold. Ranked by the ratio",
                            median(b$tpm_bound), median(b$tpm_unbound), median(b$tpm_ratio)),
         caption = "Highly expressed genes give larger, more easily detected changes, so the comparisons that follow are made within expression deciles.",
         x = "Knockdowns", y = "Median control TPM") +
    theme_viz() + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
}

# ---- b: change rate of bound and unbound genes (from 18_change_rate.R) -----------------------------------------------------
fig_rates <- function(cr) {
  r <- cr$res |> filter(version == "main") |> mutate(cell_line = cell_f(cell_line)) |>
    arrange(rate_unbound) |> mutate(rank = row_number())
  s <- cr$summ |> filter(set == "main")
  ggplot(r, aes(rank)) +
    geom_segment(aes(y = rate_unbound, yend = rate_bound), colour = GRID, linewidth = 0.7) +
    geom_point(aes(y = rate_unbound, shape = "Unbound genes"), colour = MUTED, size = 1.5) +
    geom_point(aes(y = rate_bound, colour = cell_line, shape = "Bound genes"), size = 1.8) +
    scale_colour_manual(values = CELL_COL, name = NULL) +
    scale_shape_manual(values = c(`Bound genes` = 16, `Unbound genes` = 1), name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2, override.aes = list(colour = INK2))) +
    labs(title = sprintf("Bound genes change at a higher rate than unbound genes in %d of %d knockdowns (before adjustment)",
                         sum(r$rate_bound > r$rate_unbound), nrow(r)),
         subtitle = sprintf("Share of genes significantly changed after knockdown (DESeq2 padj < 0.05, either direction). Median %s bound vs %s unbound.\nKnockdowns ranked by the rate in unbound genes",
                            scales::percent(s$median_rate_bound, 0.1), scales::percent(s$median_rate_unbound, 0.1)),
         caption = "Unadjusted rates: bound genes are more highly expressed, and changes are easier to detect at high expression (adjusted for in the next figure).",
         x = "Knockdowns", y = "Genes changed") +
    theme_viz() + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
}

# ---- c: expression-adjusted odds ratio of changing ---------------------------------------------------------------------------
fig_odds <- function(cr) {
  r <- cr$res |> filter(version == "main") |>
    mutate(cell_line = cell_f(cell_line), label = paste(target, cell_line), sig = q_mh < 0.05) |>
    arrange(or_mh) |> mutate(rank = row_number())
  s <- cr$summ |> filter(set == "main")
  pos <- cr$params$POS
  lab <- r |> filter(label %in% pos)
  ggplot(r, aes(rank, or_mh)) +
    geom_hline(yintercept = 1, colour = INK2, linewidth = 0.4, linetype = "22") +
    geom_linerange(aes(ymin = or_lo, ymax = or_hi), colour = GRID, linewidth = 1.2) +
    geom_point(aes(colour = cell_line, shape = sig), size = 1.8, stroke = 0.6) +
    geom_text_repel(data = lab, aes(label = label), size = pt_to_size(FS$annot), colour = INK2, seed = 1,
                    min.segment.length = 0, segment.colour = AXIS, nudge_x = -22, nudge_y = 0.35, box.padding = 0.4) +
    scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 5)) +
    scale_colour_manual(values = CELL_COL, name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."),
                       breaks = c("TRUE", "FALSE"), name = NULL) +
    labs(title = sprintf("At equal expression, bound genes have higher odds of changing in %d of %d knockdowns (lower in %d)",
                         s$n_sig_up, s$n_kd, s$n_sig_down),
         subtitle = sprintf("Odds ratio of changing, bound vs unbound genes, within expression deciles (Mantel–Haenszel; grey bar: 95%% CI).\nMedian %.2f (95%% CI %.2f–%.2f; signed-rank p = %s). 1 = no difference",
                            s$median_or, s$or_ci_lo, s$or_ci_hi, fmt_p(s$signed_rank_p)),
         x = "Knockdowns ranked by odds ratio", y = "Odds ratio (log scale)") +
    theme_viz() + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
}

# ---- d: magnitude, all bound vs unbound genes ------------------------------------------------------------------------------------
fig_magnitude <- function(d) {
  s <- d$summ |> filter(set == "main")
  r <- d$mag |>
    mutate(cell_line = cell_f(cell_line), label = paste(target, cell_line), sig = q < 0.05) |>
    arrange(diff_median) |> mutate(rank = row_number())
  pos <- d$params$POS
  lab <- r |> filter(label %in% pos)
  ylim <- range(c(r$diff_median, quantile(r$ci_lo, 0.01), quantile(r$ci_hi, 0.99)))   # every point visible; extreme CIs cut
  ttl <- if (s$signed_rank_p < 0.05 && s$median_diff > 0)
    sprintf("Bound genes change slightly more than unbound genes of the same expression (median %+.3f log2 FC)", s$median_diff)
  else sprintf("Bound genes do not change more than unbound genes of the same expression (median %+.3f log2 FC)", s$median_diff)
  ggplot(r, aes(rank, diff_median)) +
    geom_hline(yintercept = 0, colour = INK2, linewidth = 0.4, linetype = "22") +
    geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), colour = GRID, linewidth = 1.1) +
    geom_point(aes(colour = cell_line, shape = sig), size = 1.8, stroke = 0.6) +
    geom_text_repel(data = lab, aes(label = label), size = pt_to_size(FS$annot), colour = INK2, seed = 1,
                    min.segment.length = 0, segment.colour = AXIS, nudge_x = -25, nudge_y = 0.06, box.padding = 0.5,
                    max.overlaps = Inf) +
    coord_cartesian(ylim = ylim) +
    scale_colour_manual(values = CELL_COL, name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."),
                       breaks = c("TRUE", "FALSE"), name = NULL) +
    labs(title = ttl,
         subtitle = sprintf("Difference in median |log2 fold change| of all bound − unbound genes (changed or not), within expression deciles,\nminus the value expected with no effect (permutation mean). Grey bar: 95%% CI (long bars cut at the plot edge).\nMedian %+.4f (95%% CI %+.4f to %+.4f), p = %s; larger in %d, smaller in %d",
                            s$median_diff, s$ci_lo, s$ci_hi, fmt_p(s$signed_rank_p), s$n_sig_up, s$n_sig_down),
         x = "Knockdowns ranked by difference", y = "Difference in median |log2 FC|") +
    theme_viz() + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
}

# ---- e: hub gene groups -------------------------------------------------------------------------------------------------
fig_hubgrp <- function(h) {
  g <- h$genes |> mutate(cell_line = cell_f(cell_line))
  s <- h$grp |> mutate(cell_line = cell_f(cell_line), hub_group = factor(hub_group, h$params$GRP_LAB))
  p1 <- ggplot(s, aes(hub_group, n, fill = cell_line)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    geom_text(aes(label = scales::comma(n)), position = position_dodge(0.75), vjust = -0.4, size = pt_to_size(FS$annot), colour = INK2) +
    scale_fill_manual(values = CELL_COL, name = NULL) +
    scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.12))) +
    labs(subtitle = "Genes per group", x = "RBPs binding the gene", y = "Genes") +
    theme_viz() + theme(panel.grid.major.x = element_blank())
  p2 <- ggplot(g, aes(hub_group, ntpm, fill = cell_line)) +
    geom_boxplot(position = position_dodge(0.75), width = 0.6, outlier.shape = NA, colour = INK2, linewidth = 0.3, alpha = 0.85) +
    scale_fill_manual(values = CELL_COL, guide = "none") +
    scale_y_log10(labels = scales::label_comma()) + coord_cartesian(ylim = c(1, 1000)) +
    labs(subtitle = "Expression in control samples (box: quartiles)", x = "RBPs binding the gene", y = "Nuclear TPM (log scale)") +
    theme_viz() + theme(panel.grid.major.x = element_blank())
  (p1 | p2) + plot_annotation(
    title = "Genes bound by many RBPs are few and highly expressed",
    subtitle = sprintf("Expressed genes tested in most knockdowns, grouped by the number of profiled RBPs binding them. Median TPM: %s",
                       paste(sprintf("%s %s", s$hub_group[s$cell_line == "K562"], round(s$median_tpm[s$cell_line == "K562"])), collapse = ", ") |>
                         (\(x) paste0(x, " (K562)"))()),
    theme = theme_viz())
}

# ---- f: does hub status alone raise the chance of changing? (unbound genes) --------------------------------------------------
HUB_COL <- setNames(SEQ_BLUE[2:4], c("< 5", "5–19", "≥ 20"))
fig_hub_rates <- function(h) {
  r <- h$rates |> filter(n >= h$params$MIN_CELL) |>
    group_by(cell_line, dec, hub_group) |> summarise(rate = median(rate), n_kd = n(), .groups = "drop") |>
    filter(n_kd >= 20) |>                                             # cells seen in >= 20 knockdowns only
    mutate(cell_line = cell_f(cell_line), hub_group = factor(hub_group, h$params$GRP_LAB))
  w <- r |> dplyr::select(cell_line, dec, hub_group, rate) |> tidyr::pivot_wider(names_from = hub_group, values_from = rate)
  n_more <- sum(w$`≥ 20` > w$`< 5`, na.rm = TRUE); n_cmp <- sum(!is.na(w$`≥ 20`) & !is.na(w$`< 5`))
  ggplot(r, aes(dec, rate, colour = hub_group)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
    facet_wrap(~ cell_line) +
    scale_colour_manual(values = HUB_COL, name = "Other RBPs\nbinding the gene") +
    scale_x_continuous(breaks = 1:10) + scale_y_continuous(labels = scales::percent) +
    labs(title = sprintf("Among unbound genes, hub genes change more often than other genes of the same expression (%d of %d deciles)",
                         n_more, n_cmp),
         subtitle = "Share of genes NOT bound by the knocked-down RBP that change after knockdown (padj < 0.05), by expression decile and hub group;\nmedian over knockdowns (cells with ≥ 10 genes, shown when present in ≥ 20 knockdowns).\nHub group = number of other profiled RBPs binding the gene",
         x = "Expression decile within the knockdown (1 = lowest)", y = "Unbound genes changed") +
    theme_viz() + theme(legend.position = "right", panel.grid.major.x = element_blank())
}

# ---- g: how many extra changing targets are explained by hub status? -----------------------------------------------------------
ATT_COL <- c(`Explained by hub status` = "#eda100", `Attributable to binding` = "#6250d6")
fig_hub_attribution <- function(h) {
  s <- h$summ |> filter(set == "main")
  a <- h$att |> mutate(label = paste(target, cell_line)) |> arrange(extra_total) |> mutate(rank = row_number())
  l <- a |> dplyr::select(rank, label, `Explained by hub status` = extra_hub, `Attributable to binding` = extra_binding) |>
    tidyr::pivot_longer(-c(rank, label), names_to = "part", values_to = "n") |> mutate(part = factor(part, names(ATT_COL)))
  pos <- h$params$POS
  lab <- a |> filter(label %in% pos)
  ggplot() +
    geom_hline(yintercept = 0, colour = INK2, linewidth = 0.4) +
    geom_col(data = l, aes(rank, n, fill = part), width = 0.85) +
    geom_point(data = a, aes(rank, extra_total, shape = "Total extra"), colour = INK, size = 1.1) +
    geom_text_repel(data = lab, aes(rank, extra_total, label = label), size = pt_to_size(FS$annot), colour = INK2, seed = 1,
                    min.segment.length = 0, segment.colour = AXIS, nudge_x = -20, nudge_y = 40, box.padding = 0.4, max.overlaps = Inf) +
    scale_fill_manual(values = ATT_COL, name = NULL) +
    scale_shape_manual(values = c(`Total extra` = 16), name = NULL) +
    labs(title = sprintf("Of %s extra changing targets, %.0f%% are explained by hub status and %.0f%% by binding",
                         scales::comma(round(s$extra_total)), 100 * s$hub_share_total, 100 * (1 - s$hub_share_total)),
         subtitle = sprintf("Per knockdown: bound genes that change, minus the number expected from unbound genes of the same expression (total extra),\nsplit into the part expected from unbound genes of the same expression and hub group (hub) and the rest (binding).\nHub share over all knockdowns %.0f%% (95%% CI %.0f–%.0f%%); hub status explains more than half of the extra in %d of %d knockdowns",
                            100 * s$hub_share_total, 100 * s$hub_share_lo, 100 * s$hub_share_hi, s$n_hub_majority, s$n_extra_pos),
         caption = sprintf("Genes with a count of binding RBPs only (most small RNAs have none). Hub group: < 5 / 5–19 / ≥ 20 other RBPs binding the gene.
Null check, random gene sets of the same size: total extra %+.0f, of which binding %+.0f (a small positive bias).",
                           h$summ$extra_total[h$summ$set == "null (random bound sets)"], h$summ$extra_binding[h$summ$set == "null (random bound sets)"]),
         x = "Knockdowns ranked by total extra changing targets", y = "Extra changing bound genes") +
    theme_viz() + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
}

# ---- h: what are the hub genes? (GO, expression-matched background) --------------------------------------------------------
fig_hub_go <- function(go, n_top = 8) {
  g <- go$go |> filter(background == "matched")
  pick <- g |> group_by(ontology, cell_line) |> slice_min(p.adjust, n = n_top, with_ties = FALSE) |> ungroup() |>
    distinct(ontology, ID)
  x <- g |> semi_join(pick, by = c("ontology", "ID")) |>
    mutate(cell_line = cell_f(cell_line), ont = factor(c(BP = "Biological process", CC = "Cellular component")[ontology],
                                                         c("Biological process", "Cellular component")),
           term = stringr::str_trunc(Description, 48))
  ord <- x |> group_by(ont, term) |> summarise(f = max(fold_enrichment), .groups = "drop") |> arrange(ont, f)
  x$term <- factor(x$term, unique(ord$term))
  p1 <- ggplot(x, aes(fold_enrichment, term)) +
    geom_point(aes(size = k, colour = cell_line), alpha = 0.85) +
    facet_grid(ont ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = CELL_COL, name = NULL) +
    scale_size_area(max_size = 4.5, name = "Hub genes\nin term") +
    labs(subtitle = sprintf("GO terms enriched in hub genes vs expression-matched non-hub genes (top %d per cell line; q < 0.05)", n_top),
         x = "Fold enrichment", y = NULL) +
    theme_viz() + theme(legend.position = "right", axis.text.y = element_text(size = FS$annot, colour = INK2),
                        strip.text.y = element_text(angle = 0, hjust = 0))
  pr <- go$prof |> filter(set != "all non-hub") |>
    tidyr::pivot_longer(c(median_tpm, median_span, median_utr3), names_to = "m", values_to = "v") |>
    mutate(m = factor(c(median_tpm = "Expression (TPM)", median_span = "Gene length (kb)", median_utr3 = "3′UTR length (kb)")[m],
                      c("Expression (TPM)", "Gene length (kb)", "3′UTR length (kb)")),
           v = ifelse(m == "Expression (TPM)", v, v / 1000), cell_line = cell_f(cell_line),
           set = factor(set, c("expression-matched non-hub", "hub (≥ 20 RBPs)")))
  p2 <- ggplot(pr, aes(v, set, colour = cell_line)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.5)) +
    facet_wrap(~ m, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = CELL_COL, guide = "none") +
    expand_limits(x = 0) +
    labs(subtitle = "Medians: hub genes are longer, with longer 3′UTRs (matching on expression is incomplete at the top end)", x = NULL, y = NULL) +
    theme_viz() + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(1.6, "lines"))
  (p1 / p2) + plot_layout(heights = c(3.2, 1)) +
    plot_annotation(title = "Hub genes are long regulatory genes: development, adhesion, signalling and RNA-granule components",
                    theme = theme_viz())
}

# ---- draw everything (a-d: binding vs change; e-h: hub genes) --------------------------------------------------------------
d  <- readRDS("data/analysis/bound_vs_unbound_plotdata.rds")      # 13_bound_vs_unbound.R
cr <- readRDS("data/analysis/change_rate_plotdata.rds")           # 18_change_rate.R
h  <- readRDS("data/analysis/hub_attribution_plotdata.rds")       # 19_hub_attribution.R
go <- readRDS("data/analysis/hub_go_plotdata.rds")                # 17_hub_go.R
save_fig(fig_baseline(d), "bvu_a_baseline", FIG$baseline$w, FIG$baseline$h)
save_fig(fig_rates(cr), "bvu_b_rates", FIG$rates$w, FIG$rates$h)
save_fig(fig_odds(cr), "bvu_c_odds", FIG$odds$w, FIG$odds$h)
save_fig(fig_magnitude(d), "bvu_d_magnitude", FIG$magnitude$w, FIG$magnitude$h)
save_fig(fig_hubgrp(h), "bvu_e_hub_groups", FIG$hubgrp$w, FIG$hubgrp$h)
save_fig(fig_hub_rates(h), "bvu_f_hub_rates", FIG$hubrates$w, FIG$hubrates$h)
save_fig(fig_hub_attribution(h), "bvu_g_hub_attribution", FIG$attrib$w, FIG$attrib$h)
save_fig(fig_hub_go(go), "bvu_h_hub_go", FIG$hubgo$w, FIG$hubgo$h)
message("Figures written to ", FIG_DIR)
