# Aim 1, part B: binding landscape FIGURES (no computation; reads data/analysis/aim1_landscape_plotdata.rds).
# Shared colours / font scale (FS) / theme / saving: analysis/R/lib_plot.R. Per-figure sizes and choices: FIG block below.
#
# Usage: Rscript analysis/R/05_binding_landscape_figures.R   (run 05_binding_landscape.R first)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(ggrepel); library(ComplexHeatmap); library(circlize) })
source("analysis/R/lib_plot.R")

# ---- per-figure settings --------------------------------------------------------------------------------------
# w = width class on the website ("text" = text column, "page" = .column-page; inches set in lib_plot.R so that
# every figure is drawn at its display size). h = height in inches: the only size to tune per figure.
FIG <- list(
  composition = list(w = "text", h = 5.0, min_label = 0.08),
  heatmap     = list(w = "text", h = 10.5, z_cap = 3),
  pca         = list(w = "page", h = 5.3),
  similarity  = list(w = "page", h = 7.6),
  concordance = list(w = "page", h = 5.2)
)
KNOWN_RBPS <- c("UPF1", "PUM2", "RBFOX2", "PTBP1", "HNRNPC", "TARDBP", "U2AF2", "DKC1", "SMNDC1", "EIF3G",
                "LARP4", "SRSF1", "QKI", "IGF2BP1")
cat_lab_n <- function(meta) {                      # "CDS (n = 44)" labels in category order
  n <- table(meta$main_category)
  setNames(sprintf("%s (n = %d)", CATEGORY_LAB[names(n)], as.integer(n)), names(n))
}
hm_legend <- list(title_gp = gpar(fontsize = FS$legend, col = INK2), labels_gp = gpar(fontsize = FS$annot, col = MUTED))
annot_gp <- function(col = INK2) gpar(fontsize = FS$annot, col = col)

# ---- a: region composition per main category (tall segmented bars; legend kept) -------------------------------------
fig_composition <- function(d) {
  lab <- cat_lab_n(d$meta)
  pf <- d$peak_frac |>
    mutate(region = factor(region, rev(REGIONS)),
           ylab = factor(lab[as.character(main_category)], rev(lab[CATEGORIES])),
           txt = ifelse(frac >= FIG$composition$min_label, sprintf("%.0f%%", 100 * frac), ""),
           txt_col = ifelse(region %in% c("small_ncRNA", "5UTR"), INK, "#ffffff"))
  ggplot(pf, aes(frac, ylab, fill = region)) +
    geom_col(width = 0.66, colour = SURFACE, linewidth = 0.6) +
    geom_text(aes(label = txt, colour = txt_col, group = region), position = position_stack(vjust = 0.5),
              size = pt_to_size(FS$annot)) +
    scale_colour_identity() +
    scale_fill_manual(values = REGION_COL, labels = REGION_LAB[REGIONS], breaks = REGIONS, name = NULL) +
    scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
    labs(title = "Where RBPs bind: each RBP has a main binding region",
         subtitle = sprintf("Mean share of peaks per region, by main binding category: the region with the highest\npeak share relative to the average RBP (a lead of < %.2f log2 = no clear preference)", d$margin),
         x = "Share of peaks on expressed genes", y = NULL) +
    theme_viz() + theme(panel.grid.major.y = element_blank()) + guides(fill = guide_legend(nrow = 2))
}

# ---- b: z-score heatmap split by category --------------------------------------------------------------------------------
fig_heatmap <- function(d) {
  m <- d$meta; cap <- FIG$heatmap$z_cap
  o <- order(m$main_category, -m$category_margin)
  Zp <- pmin(pmax(d$Z[o, ], -cap), cap); colnames(Zp) <- REGION_LAB[colnames(Zp)]; mo <- m[o, ]
  show <- which(mo$target %in% KNOWN_RBPS)
  ra <- rowAnnotation(
    `Cell line` = mo$cell_line,
    Peaks = anno_barplot(log10(mo$n_peaks), gp = gpar(fill = MUTED, col = NA), border = FALSE, width = unit(1.4, "cm"),
                         axis_param = list(at = 2:4, labels = c("100", "1K", "10K"), gp = annot_gp(MUTED))),
    RBP = anno_mark(at = show, labels = mo$label[show], labels_gp = annot_gp()),
    col = list(`Cell line` = CELL_COL), annotation_legend_param = list(`Cell line` = hm_legend),
    annotation_name_gp = annot_gp(), annotation_name_side = "bottom", annotation_name_rot = 45)
  la <- rowAnnotation(Category = mo$main_category, col = list(Category = CATEGORY_COL), show_legend = FALSE,
                      show_annotation_name = FALSE, simple_anno_size = unit(2.5, "mm"))
  ht <- Heatmap(Zp, name = "z-score", col = colorRamp2(c(-cap, 0, cap), c(DIV_BLUE, MID, DIV_RED)),
                cluster_columns = FALSE, cluster_rows = FALSE,
                row_split = factor(CATEGORY_LAB[as.character(mo$main_category)], CATEGORY_LAB[CATEGORIES]),
                row_title_rot = 0, row_title_gp = annot_gp(),
                show_row_names = FALSE, row_gap = unit(1.5, "pt"),
                column_names_gp = gpar(fontsize = FS$axis_text, col = INK2), column_names_rot = 45,
                left_annotation = la, right_annotation = ra,
                heatmap_legend_param = c(hm_legend, list(at = c(-cap, 0, cap), labels = c(paste0("≤−", cap), "0", paste0("≥", cap)))),
                column_title = sprintf("Region enrichment of all %d experiments (z-scored)\nRows grouped by main category, strongest preference first", nrow(Zp)),
                column_title_gp = gpar(fontsize = FS$title, fontface = "bold", col = INK))
  function() draw(ht, background = SURFACE, merge_legend = TRUE)
}

# ---- c: PCA map, one panel per main binding category (replaces UMAP) ----------------------------------------------------
# Axis titles carry the variance explained and what each end of the axis means (from the loadings; see pca_loadings.tsv).
fig_pca <- function(d) {
  m <- d$meta; v <- d$pca$var
  xl <- sprintf("← coding sequence          PC1 (%.0f%% of variance)          introns →", v$pct_variance[1])
  yl <- sprintf("← 3′UTR, deep intron     PC2 (%.0f%%)     5′UTR, small ncRNA →", v$pct_variance[2])
  facet_df <- bind_rows(lapply(CATEGORIES, function(g) mutate(m, panel = g, hit = main_category == g))) |>
    mutate(panel = factor(panel, CATEGORIES, labels = cat_lab_n(m)[CATEGORIES]))
  ggplot(facet_df, aes(pc1, pc2)) +
    geom_hline(yintercept = 0, colour = GRID, linewidth = 0.3) + geom_vline(xintercept = 0, colour = GRID, linewidth = 0.3) +
    geom_point(data = ~ filter(.x, !hit), colour = GRID, size = 1.1) +
    geom_point(data = ~ filter(.x, hit), aes(fill = as.character(main_category)), shape = 21, colour = SURFACE,
               size = 2.2, stroke = 0.4) +
    scale_fill_manual(values = CATEGORY_COL, guide = "none") +
    facet_wrap(~ panel, nrow = 2) +
    labs(title = sprintf("Two axes explain %.0f%% of binding variation: introns vs coding sequence, and 5′UTR/small RNA vs 3′UTR",
                         sum(v$pct_variance[1:2])),
         subtitle = "PCA of region-enrichment z-scores (231 experiments, 7 regions). Each panel highlights one main binding category (grey = all others).",
         x = xl, y = yl) +
    theme_viz() + theme(axis.text = element_blank(), panel.grid.major = element_blank()) + coord_equal()
}

# ---- d: experiment x experiment similarity -----------------------------------------------------------------------------------
fig_similarity <- function(d) {
  m <- d$meta; o <- d$within_order; S <- d$S[o, o]; mo <- m[o, ]
  col_fun <- colorRamp2(c(-1, 0, 1), c(DIV_BLUE, MID, DIV_RED))
  split <- factor(CATEGORY_LAB[as.character(mo$main_category)], CATEGORY_LAB[CATEGORIES])
  ann <- function(which) HeatmapAnnotation(
    which = which, Category = mo$main_category, `Cell line` = mo$cell_line,
    col = list(Category = CATEGORY_COL, `Cell line` = CELL_COL), simple_anno_size = unit(2.5, "mm"),
    show_annotation_name = which == "column", annotation_name_gp = annot_gp(),
    show_legend = which == "column", annotation_legend_param = list(Category = c(hm_legend, list(labels = CATEGORY_LAB[CATEGORIES], at = CATEGORIES)),
                                                                    `Cell line` = hm_legend))
  h1 <- Heatmap(S, name = "Correlation", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                row_split = split, column_split = split, row_gap = unit(1, "pt"), column_gap = unit(1, "pt"),
                row_title_rot = 0, row_title_gp = annot_gp(), column_title = NULL,
                show_row_names = FALSE, show_column_names = FALSE, top_annotation = ann("column"), left_annotation = ann("row"),
                heatmap_legend_param = c(hm_legend, list(at = c(-1, 0, 1))), use_raster = TRUE)
  c <- d$conc; ko <- c$kh_order; Skh <- c$S_kh[ko, ko]; rbps <- c$both[ko]
  kcat <- m$main_category[which(m$cell_line == "K562")[match(rbps, m$target[m$cell_line == "K562"])]]
  mark <- which(rbps %in% KNOWN_RBPS)                       # label a few known RBPs instead of all 70
  h2 <- Heatmap(Skh, name = "Correlation ", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE,
                left_annotation = rowAnnotation(K562 = kcat, col = list(K562 = CATEGORY_COL), show_legend = FALSE,
                                                show_annotation_name = FALSE, simple_anno_size = unit(2.5, "mm")),
                right_annotation = rowAnnotation(RBP = anno_mark(at = mark, labels = rbps[mark], labels_gp = annot_gp())),
                row_title = "K562 experiment", column_title = "HepG2 experiment (same RBP order)",
                row_title_gp = annot_gp(), column_title_gp = annot_gp(), column_title_side = "bottom",
                layer_fun = function(j, i, x, y, w, h, fill) {
                  on <- i == j
                  if (any(on)) grid.rect(x[on], y[on], w[on], h[on], gp = gpar(col = INK, fill = NA, lwd = 0.6))
                })
  same <- diag(c$S_kh); other <- c$S_kh[row(c$S_kh) != col(c$S_kh)]
  function() {
    grid.newpage(); grid.rect(gp = gpar(fill = SURFACE, col = NA))
    grid.text("The same RBP has very similar binding profiles in both cell lines; different RBPs share little",
              x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"), just = c("left", "top"),
              gp = gpar(fontsize = FS$title, fontface = "bold", col = INK))
    grid.text(sprintf("Pearson correlation of region-enrichment profiles. Left: all %d experiments, grouped by main category.\nRight: the %d RBPs profiled in both lines; outlined diagonal = same RBP (median r %.2f vs %.2f for different RBPs).",
                      nrow(S), length(c$both), median(same), median(other)),
              x = unit(4, "mm"), y = unit(1, "npc") - unit(11, "mm"), just = c("left", "top"),
              gp = gpar(fontsize = FS$subtitle, col = INK2, lineheight = 1.1))
    pushViewport(viewport(x = 0, y = 0, width = 0.57, height = 0.86, just = c("left", "bottom")))
    draw(h1, newpage = FALSE, background = SURFACE, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    popViewport()
    pushViewport(viewport(x = 0.58, y = 0.04, width = 0.41, height = 0.82, just = c("left", "bottom")))
    draw(h2, newpage = FALSE, background = SURFACE)
    popViewport()
  }
}

# ---- e: cross-cell concordance (profile distance + category agreement) -----------------------------------------------------------
fig_concordance <- function(d) {
  c <- d$conc
  pa <- ggplot() +
    geom_histogram(data = data.frame(x = c$diff), aes(x, after_stat(density)), bins = 40, fill = AXIS, colour = SURFACE, linewidth = 0.3) +
    geom_rug(data = data.frame(x = c$same), aes(x), colour = HIGHLIGHT, linewidth = 0.6, length = unit(0.06, "npc")) +
    geom_vline(xintercept = c$obs, colour = HIGHLIGHT, linewidth = 0.6) +
    geom_vline(xintercept = mean(c$diff), colour = INK2, linewidth = 0.4) +
    # labels sit in the empty upper-right corner (they would cross the lines at narrow widths)
    annotate("text", x = Inf, y = Inf, vjust = 1.3, hjust = 1.02, size = pt_to_size(FS$annot), colour = INK,
             label = sprintf("Blue line: same RBP, K562 vs HepG2\nmean %.2f (n = %d; ticks = each RBP)", c$obs, length(c$both))) +
    annotate("text", x = Inf, y = Inf, vjust = 4.6, hjust = 1.02, size = pt_to_size(FS$annot), colour = INK2,
             label = sprintf("Grey line: different RBPs\nmean %.2f (bars)", mean(c$diff))) +
    labs(subtitle = sprintf("Profile distance between cell lines\npermutation p %s (10,000 shuffles)",
                            if (c$p_perm < 1e-3) "< 0.001" else sprintf("= %.3f", c$p_perm)),
         x = "Euclidean distance between region profiles (z-scores)", y = "Density") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +          # headroom for the labels
    theme_viz() + theme(panel.grid.major.x = element_blank())
  tab <- c$category |> count(category_K562, category_HepG2) |>
    complete(category_K562 = factor(CATEGORIES, CATEGORIES), category_HepG2 = factor(CATEGORIES, CATEGORIES), fill = list(n = 0)) |>
    mutate(diag = as.character(category_K562) == as.character(category_HepG2))
  agree <- sum(tab$n[tab$diag]); mixed <- sum(tab$n[!tab$diag & (tab$category_K562 == "mixed" | tab$category_HepG2 == "mixed")])
  tab$category_K562 <- factor(tab$category_K562, rev(CATEGORIES))
  pb <- ggplot(tab, aes(category_HepG2, category_K562)) +
    geom_tile(aes(fill = ifelse(n == 0, NA, n)), colour = SURFACE, linewidth = 0.8) +
    geom_tile(data = ~ filter(.x, diag), fill = NA, colour = INK2, linewidth = 0.4) +
    geom_text(aes(label = ifelse(n > 0, n, "")), size = pt_to_size(FS$annot), colour = INK) +
    scale_fill_gradientn(colours = SEQ_BLUE[1:3], na.value = MID, guide = "none") +
    scale_x_discrete(labels = CATEGORY_LAB) + scale_y_discrete(labels = CATEGORY_LAB) +
    labs(subtitle = sprintf("Main category: same in %d/%d RBPs\n(%d of %d differences involve 'no clear preference')",
                            agree, sum(tab$n), mixed, sum(tab$n) - agree), x = "HepG2", y = "K562") +
    theme_viz() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
  (pa | pb) + plot_layout(widths = c(1.15, 1)) +
    plot_annotation(title = "Binding profiles are reproducible across cell lines",
                    subtitle = "The same RBP's K562 and HepG2 profiles are far closer than different RBPs'; the coarse category label agrees less often,\nmostly at the 'no clear preference' boundary",
                    theme = theme_viz())
}

# ---- draw everything ----------------------------------------------------------------------------------------------------------------
d <- readRDS("data/analysis/aim1_landscape_plotdata.rds")
d$meta$cell_line <- factor(d$meta$cell_line, names(CELL_COL))
save_fig(fig_composition(d), "aim1_landscape_a_region_composition", FIG$composition$w, FIG$composition$h)
save_drawing(fig_heatmap(d), "aim1_landscape_b_enrichment_heatmap", FIG$heatmap$w, FIG$heatmap$h)
save_fig(fig_pca(d), "aim1_landscape_c_pca", FIG$pca$w, FIG$pca$h)
save_drawing(fig_similarity(d), "aim1_landscape_d_similarity", FIG$similarity$w, FIG$similarity$h)
save_fig(fig_concordance(d), "aim1_landscape_e_concordance", FIG$concordance$w, FIG$concordance$h)
message("Landscape figures written to ", FIG_DIR)
