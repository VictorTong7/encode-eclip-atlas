# Aim 1, part A: dataset overview FIGURES (no computation; reads data/analysis/aim1_overview_plotdata.rds).
# Shared colours / theme / saving: analysis/R/lib_plot.R. Per-figure sizes and choices: FIG block below.
#
# Usage: Rscript analysis/R/05_dataset_overview_figures.R   (run 05_dataset_overview.R first)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(ggrepel); library(ComplexHeatmap) })
source("analysis/R/lib_plot.R")

# ---- per-figure settings --------------------------------------------------------------------------------
# w = width class on the website ("text" = text column, "page" = .column-page; inches set in lib_plot.R so that
# every figure is drawn at its display size). h = height in inches: the only size to tune per figure.
FIG <- list(
  profiled     = list(w = "text", h = 4.8),
  peaks        = list(w = "text", h = 4.6),
  properties   = list(w = "page", h = 3.9, bins = 60),
  batch        = list(w = "page", h = 4.7),
  genes        = list(w = "page", h = 4.7),
  shared       = list(w = "page", h = 8.2, n_hubs = 15, bins = 45)
)
SET_LAB <- c(eCLIP_K562 = "eCLIP · K562", eCLIP_HepG2 = "eCLIP · HepG2", knockdown_K562 = "Knockdown · K562",
             knockdown_HepG2 = "Knockdown · HepG2", RBNS = "RNA Bind-n-Seq")
lab_cell <- function(p) p + scale_colour_manual(values = CELL_COL, name = NULL)
fmt_p <- function(p) if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)

# ---- a: what was profiled? ------------------------------------------------------------------------------------
fig_profiled <- function(d) {
  m <- as.matrix(d$avail[, names(SET_LAB)]); rownames(m) <- d$avail$target; colnames(m) <- SET_LAB[colnames(m)]
  cm <- make_comb_mat(m)
  cm <- cm[comb_size(cm) >= 2]                                   # drop combinations of a single RBP
  e <- d$exps
  a <- d$avail
  n_kd <- sum((a$eCLIP_K562 & a$knockdown_K562) | (a$eCLIP_HepG2 & a$knockdown_HepG2))
  ttl <- sprintf("%d RBPs profiled by eCLIP: %d with a matching knockdown, %d with cell-free data",
                 n_distinct(e$target), n_kd, sum(a$RBNS))
  cap <- sprintf("Bars: number of RBPs with exactly that combination of datasets (combinations with ≥ 2 RBPs).
%d eCLIP sets lifted over from GRCh38; %d RBP × cell line pairs have two experiments.",
                 sum(e$lifted), n_distinct(paste(e$target, e$cell_line)[e$n_same_rbp_cell > 1]))
  ht <- UpSet(cm, comb_order = order(-comb_size(cm)), set_order = names(SET_LAB) |> (\(x) SET_LAB[x])(),
              pt_size = unit(2.8, "mm"), lwd = 1.6, comb_col = INK2, bg_col = c(MID, SURFACE), bg_pt_col = GRID,
              top_annotation = upset_top_annotation(cm, add_numbers = TRUE, numbers_rot = 0, gp = gpar(fill = HIGHLIGHT, col = NA),
                                                    numbers_gp = gpar(fontsize = FS$annot, col = INK2),
                                                    annotation_name_gp = gpar(fontsize = FS$annot, col = INK2),
                                                    axis_param = list(gp = gpar(fontsize = FS$annot, col = MUTED))),
              right_annotation = upset_right_annotation(cm, add_numbers = TRUE, gp = gpar(fill = MUTED, col = NA),
                                                        numbers_gp = gpar(fontsize = FS$annot, col = INK2),
                                                        annotation_name_gp = gpar(fontsize = FS$annot, col = INK2),
                                                        axis_param = list(gp = gpar(fontsize = FS$annot, col = MUTED))),
              row_names_gp = gpar(fontsize = FS$axis_text, col = INK2),
              column_title = ttl, column_title_gp = gpar(fontsize = FS$title, fontface = "bold", col = INK), column_title_side = "top")
  function() {
    draw(ht, background = SURFACE, padding = unit(c(16, 4, 4, 4), "mm"))
    grid::grid.text(cap, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"),
                    gp = gpar(fontsize = FS$annot, col = MUTED))
  }
}

# ---- b: how many peaks per experiment? ---------------------------------------------------------------------------
fig_peaks <- function(d) {
  e <- d$exps |> arrange(desc(n_peaks)) |> mutate(rank = row_number())
  th <- d$thresholds
  lines <- tibble(y = th, lab = c(sprintf("%s peaks · %d experiments below (excluded from binding profiles)", th[1],
                                          sum(e$n_peaks < th[1])),
                                  sprintf("%s peaks · %d below (excluded from sequence models)",
                                          scales::comma(th[2]), sum(e$n_peaks < th[2]))))
  p <- ggplot(e, aes(rank, n_peaks)) +
    geom_hline(data = lines, aes(yintercept = y), colour = AXIS, linewidth = 0.4, linetype = "22") +
    geom_text(data = lines, aes(x = 1, y = y, label = lab), hjust = 0, vjust = 1.5, size = pt_to_size(FS$annot), colour = INK2) +
    geom_point(aes(colour = cell_line, shape = lifted), size = 1.8, stroke = 0.7) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = "hg19 native", `TRUE` = "lifted from GRCh38"),
                       name = NULL) +
    scale_y_log10(labels = scales::label_comma(), breaks = c(10, 100, 1000, 10000)) +
    labs(title = sprintf("Experiment depth varies %s-fold", scales::comma(max(e$n_peaks) / min(e$n_peaks), accuracy = 100)),
         subtitle = sprintf("Peaks per eCLIP experiment (IDR peaks, log scale), ranked; median %s",
                            scales::comma(median(e$n_peaks))),
         x = "Experiments ranked by peak count", y = "Peaks") +
    theme_viz() + theme(panel.grid.major.x = element_blank())
  lab_cell(p) + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
}

# ---- c: what do the peaks look like? --------------------------------------------------------------------------------
fig_properties <- function(d) {
  v <- d$peak_vals
  one <- function(col, xlab, ttl, log = FALSE, cut = NULL, cut_lab = NULL) {
    p <- ggplot(v, aes(.data[[col]], after_stat(density), colour = cell_line)) +
      geom_freqpoly(bins = FIG$properties$bins, linewidth = 0.6) +
      labs(subtitle = ttl, x = xlab, y = if (col == "width") "Density" else NULL) + theme_viz() +
      theme(legend.position = if (col == "width") "top" else "none", axis.text.y = element_blank(),
            panel.grid.major.y = element_blank())
    if (log) p <- p + scale_x_log10(labels = scales::label_comma())
    if (!is.null(cut)) p <- p + geom_vline(xintercept = cut, colour = INK2, linewidth = 0.4, linetype = "22") +
      annotate("text", x = cut, y = Inf, label = cut_lab, hjust = -0.05, vjust = 1.4, size = pt_to_size(FS$annot), colour = INK2)
    lab_cell(p)
  }
  pw <- one("width", "Peak width (nt)", sprintf("Width: median %d nt", as.integer(median(v$width))), log = TRUE)
  pf <- one("log2fc", "log2 enrichment over input", "Enrichment: hard cut at 3", cut = 3, cut_lab = "ENCODE cut-off")
  pp <- one("neglog10p", "−log10 p", sprintf("Significance: cut at 3, cap at 400 (%.1f%%)",
                                              100 * mean(v$neglog10p >= 400)), log = TRUE, cut = 3, cut_lab = "cut-off")
  (pw | pf | pp) + plot_annotation(
    title = "Peaks are short and sit just above ENCODE's thresholds",
    subtitle = sprintf("Distributions over all %s peaks, by cell line", scales::comma(nrow(v))),
    theme = theme_viz())
}

# ---- d: is depth driven by batch? ----------------------------------------------------------------------------------
fig_batch <- function(d) {
  e <- d$exps |> mutate(year = factor(release_year))
  b <- d$batch
  pv <- function(cl, m) fmt_p(b$wilcoxon_p[b$cell_line == cl & b$metric == m])
  base <- function(y, ylab, sub, logy = FALSE) {
    p <- ggplot(e, aes(year, .data[[y]], colour = cell_line)) +
      geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.6, seed = 1), size = 1.5, alpha = 0.8) +
      stat_summary(aes(group = cell_line), fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3,
                   position = position_dodge(0.6), colour = INK) +
      labs(subtitle = sub, x = "Release year", y = ylab) +
      theme_viz() + theme(panel.grid.major.x = element_blank())
    if (logy) p <- p + scale_y_log10(labels = scales::label_comma())
    lab_cell(p)
  }
  p1 <- base("n_peaks", "Peaks per experiment", sprintf("Depth · ENCODE4 vs ENCODE3: K562 p %s, HepG2 p %s",
                                                       pv("K562", "n_peaks"), pv("HepG2", "n_peaks")), logy = TRUE)
  p2 <- base("median_log2fc", "Median peak log2 enrichment",
             sprintf("Peak strength · K562 p %s, HepG2 p %s", pv("K562", "median_log2fc"), pv("HepG2", "median_log2fc"))) +
    theme(legend.position = "none")
  (p1 | p2) + plot_annotation(title = "Later (ENCODE4) K562 experiments are shallower; peak strength is stable",
                              subtitle = "One point per experiment; bars = median. ENCODE3 = 2015–2018, ENCODE4 = 2018–2023; Wilcoxon tests by phase within each cell line",
                              theme = theme_viz())
}

# ---- e: how many genes does each RBP bind? -------------------------------------------------------------------------
fig_genes <- function(d) {
  e <- d$exps
  p1 <- ggplot(e, aes(n_peaks, n_genes, colour = cell_line)) +
    geom_abline(slope = 1, intercept = 0, colour = AXIS, linewidth = 0.4, linetype = "22") +
    annotate("text", x = 20, y = 20, label = "1 peak per gene", angle = 38, vjust = -0.6, size = pt_to_size(FS$annot), colour = MUTED) +
    geom_point(size = 1.6, alpha = 0.85) +
    scale_x_log10(labels = scales::label_comma()) + scale_y_log10(labels = scales::label_comma()) +
    labs(subtitle = sprintf("Bound genes vs peaks; median %.1f peaks per bound gene", median(e$n_peaks / e$n_genes)),
         x = "Peaks", y = "Bound genes") + theme_viz()
  qc <- e |> select(experiment, cell_line, `On expressed genes` = pct_expressed_gene, Intergenic = pct_intergenic) |>
    pivot_longer(c(`On expressed genes`, Intergenic)) |> mutate(name = factor(name, c("On expressed genes", "Intergenic")))
  p2 <- ggplot(qc, aes(cell_line, value, colour = cell_line)) +
    geom_point(position = position_jitter(width = 0.25, seed = 1), size = 1.2, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.6, linewidth = 0.3, colour = INK) +
    facet_wrap(~ name, scales = "free_y") +
    labs(subtitle = "QC: % of each experiment's peaks", x = NULL, y = "% of peaks") +
    theme_viz() + theme(legend.position = "none", panel.grid.major.x = element_blank())
  (lab_cell(p1) | lab_cell(p2)) + plot_layout(widths = c(1.2, 1)) +
    plot_annotation(title = "Most RBPs bind thousands of genes, almost all of them expressed",
                    subtitle = sprintf("Median per experiment: %s bound genes, %.0f%% of peaks on expressed genes, %.1f%% intergenic",
                                       scales::comma(median(e$n_genes)), median(e$pct_expressed_gene), median(e$pct_intergenic)),
                    theme = theme_viz())
}

# ---- f: do many RBPs bind the same genes? -----------------------------------------------------------------------------
fig_shared <- function(d) {
  g <- d$genes; r <- d$gene_rho; n <- FIG$shared$n_hubs
  p1 <- ggplot(g, aes(frac_rbps, after_stat(density), colour = cell_line)) +
    geom_freqpoly(binwidth = 0.025, boundary = 0, linewidth = 0.6) +
    scale_x_continuous(labels = scales::percent) +
    labs(subtitle = sprintf("Share of profiled RBPs binding each expressed gene
(%.0f%% bound by none)", mean(r$pct_unbound)),
         x = "Share of the cell line's RBPs", y = "Density") + theme_viz() + theme(axis.text.y = element_blank())
  hubs <- g |> group_by(gene_name) |> filter(n() == 2) |> summarise(frac = mean(frac_rbps), .groups = "drop") |>
    slice_max(frac, n = n, with_ties = FALSE) |> pull(gene_name)
  hd <- g |> filter(gene_name %in% hubs) |> mutate(gene_name = factor(gene_name, rev(hubs)))
  p2 <- ggplot(hd, aes(frac_rbps, gene_name, colour = cell_line)) +
    geom_line(aes(group = gene_name), colour = GRID, linewidth = 1) + geom_point(size = 2) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(subtitle = sprintf("Top %d 'hub' genes
(bound by the most RBPs in both lines)", n), x = "Share of RBPs", y = NULL) +
    theme_viz() + theme(legend.position = "none", panel.grid.major.y = element_blank())
  hex <- function(x, xlab, sub) {
    ggplot(g, aes(.data[[x]], n_rbps)) + geom_bin_2d(bins = FIG$shared$bins) +
      scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                    breaks = scales::breaks_log(n = 4)) +
      scale_fill_gradientn(colours = SEQ_BLUE, trans = "log10", name = "Genes") +
      facet_wrap(~ cell_line) + labs(subtitle = sub, x = xlab, y = "RBPs binding the gene") +
      theme_viz() + theme(legend.position = "right")
  }
  p3 <- hex("median_tpm", "Expression in control RNA-seq (TPM)",
            sprintf("Expression: Spearman ρ %.2f (K562), %.2f (HepG2)
 ",
                    r$spearman_rbps_vs_tpm[r$cell_line == "K562"], r$spearman_rbps_vs_tpm[r$cell_line == "HepG2"]))
  p4 <- hex("gene_span", "Gene length (nt)",
            sprintf("Length: ρ %.2f (K562), %.2f (HepG2)
Joint model β (K562): expression %.2f, length %.2f",
                    r$spearman_rbps_vs_span[r$cell_line == "K562"], r$spearman_rbps_vs_span[r$cell_line == "HepG2"],
                    r$beta_log_tpm[r$cell_line == "K562"], r$beta_log_span[r$cell_line == "K562"]))
  ((lab_cell(p1) | lab_cell(p2)) / (p3 | p4)) + plot_layout(heights = c(1, 1)) +
    plot_annotation(title = "Highly expressed genes are bound by many RBPs: eCLIP detection follows expression",
                    subtitle = "Expressed genes only (median TPM ≥ 1); small ncRNAs excluded. 'Hub' genes are mainly abundant transcripts, partly long ones.",
                    theme = theme_viz())
}

# ---- draw everything ----------------------------------------------------------------------------------------------------
d <- readRDS("data/analysis/aim1_overview_plotdata.rds")
for (nm in c("exps", "genes", "peak_vals")) d[[nm]]$cell_line <- factor(d[[nm]]$cell_line, names(CELL_COL))   # K562 first
save_drawing(fig_profiled(d), "aim1_overview_a_profiled", FIG$profiled$w, FIG$profiled$h)
save_fig(fig_peaks(d), "aim1_overview_b_peaks_per_experiment", FIG$peaks$w, FIG$peaks$h)
save_fig(fig_properties(d), "aim1_overview_c_peak_properties", FIG$properties$w, FIG$properties$h)
save_fig(fig_batch(d), "aim1_overview_d_batch", FIG$batch$w, FIG$batch$h)
save_fig(fig_genes(d), "aim1_overview_e_genes_per_experiment", FIG$genes$w, FIG$genes$h)
save_fig(fig_shared(d), "aim1_overview_f_shared_genes", FIG$shared$w, FIG$shared$h)
message("Overview figures written to ", FIG_DIR)
