# Binding vs expression, step 1: knockdown QC FIGURES (no computation; reads data/analysis/kd_qc_plotdata.rds).
# Shared colours / theme / saving: analysis/R/lib_plot.R. Per-figure sizes and choices: FIG block below.
#
# Usage: Rscript analysis/R/12_kd_qc_figures.R   (run 12_kd_qc.R first)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(grid) })
source("analysis/R/lib_plot.R")

# ---- per-figure settings --------------------------------------------------------------------------------
FIG <- list(
  samples    = list(w = "page", h = 6.6),
  efficiency = list(w = "page", h = 4.6, label = c("NOL12 HepG2", "FXR2 K562")),
  response   = list(w = "text", h = 5.4, n_label = 2),
  eclip      = list(w = "page", h = 5.4, label = c("HLTF HepG2", "SRSF9 HepG2")),
  flow       = list(w = "text", h = 5.8)
)
GROUP_COL <- c(CATEGORY_COL, uncategorised = INK2)
GROUP_LAB <- c(CATEGORY_LAB, uncategorised = "Uncategorised (< 100 peaks)")
STATUS <- c(none = "Analysed", `weak knockdown` = "Excluded: weak knockdown", `ncRNA binder` = "Excluded: ncRNA binder")
STATUS_SHAPE <- c(16, 1, 4)
pct <- function(x) sprintf("%.0f%%", 100 * x)

# ---- a: RNA-seq sample QC ----------------------------------------------------------------------------------
fig_samples <- function(d) {
  s <- d$samples |> filter(shrna) |> mutate(cell_line = factor(cell_line, names(CELL_COL))) |>
    arrange(cell_line, control_batch, desc(kind), experiment, replicate) |>
    group_by(cell_line) |> mutate(x = row_number()) |> ungroup()
  mets <- c(counts_M = "Reads assigned to genes (M)", genes_c10 = "Genes with ≥ 10 counts", chrM_share = "Mitochondrial share of reads")
  l <- s |> dplyr::select(cell_line, x, kind, experiment, target, qc_flag, all_of(names(mets))) |>
    pivot_longer(all_of(names(mets)), names_to = "metric") |> mutate(metric = factor(mets[metric], mets))
  odd <- s |> filter(qc_flag != "") |> group_by(cell_line, control_batch) |> summarise(x0 = min(x) - 0.5, x1 = max(x) + 0.5, .groups = "drop")
  fx <- l |> filter(experiment == "ENCSR139BIJ")
  fxl <- fx |> group_by(cell_line, metric) |> slice_head(n = 1) |> ungroup() |> filter(metric == mets[["counts_M"]])
  ggplot(l, aes(x, value)) +
    geom_rect(data = odd, aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = MID) +
    geom_point(aes(colour = cell_line, shape = kind), size = 1.3, stroke = 0.5) +
    geom_point(data = fx, shape = 21, size = 2.6, colour = INK, stroke = 0.6) +
    geom_text_repel(data = fxl, aes(label = "FXR2 K562 (kept, flagged)"), size = pt_to_size(FS$annot), colour = INK2,
                    nudge_y = 42, nudge_x = 12, hjust = 0, min.segment.length = 0, segment.colour = AXIS, seed = 1) +
    facet_grid(metric ~ cell_line, scales = "free", space = "free_x", switch = "y") +
    scale_colour_manual(values = CELL_COL, guide = "none") +
    scale_shape_manual(values = c(knockdown = 16, control = 2), labels = c(knockdown = "Knockdown replicate", control = "Control replicate"), name = NULL) +
    scale_y_continuous(labels = function(v) ifelse(v < 1, scales::percent(v, accuracy = 1), scales::comma(v))) +
    labs(title = "No sample is too shallow: every replicate detects over 12,000 genes",
         subtitle = sprintf("%d shRNA RNA-seq replicates, grouped by control batch. Shaded: %d small-RNA-rich K562 control batches (≥ %s of TPM in genes < %d nt)",
                            nrow(s), length(d$odd_batches), pct(0.4), d$params$SHORT_NT),
         caption = "Only FXR2 K562 from the shaded batches is analysed (the others are weaker duplicate knockdowns).
Mitochondrial and small RNAs distort TPM, so TPM-based values use nuclear TPM (chrM genes removed, renormalised).",
         x = NULL, y = NULL) +
    theme_viz() + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), strip.placement = "outside",
                        strip.text.y.left = element_text(angle = 0, hjust = 1, colour = INK, face = "bold", size = FS$strip))
}

# ---- b: how strong is each knockdown? ----------------------------------------------------------------------
fig_efficiency <- function(d) {
  k <- d$kd |> mutate(cell_line = factor(cell_line, names(CELL_COL)), label = paste(target, cell_line)) |>
    group_by(cell_line) |> arrange(kd_efficiency, .by_group = TRUE) |> mutate(rank = row_number()) |>
    mutate(strip = sprintf("%s · %d kept, %d excluded (< %s)", cell_line, sum(!weak_kd), sum(weak_kd), pct(d$params$MIN_KD))) |>
    ungroup()
  k$strip <- factor(k$strip, unique(k$strip[order(k$cell_line)]))
  lab <- k |> group_by(cell_line) |> filter(rank %in% c(1, max(rank)) | label %in% FIG$efficiency$label) |> ungroup()
  ggplot(k, aes(rank, kd_efficiency)) +
    geom_hline(yintercept = d$params$MIN_KD, colour = INK2, linewidth = 0.4, linetype = "22") +
    geom_point(data = filter(k, weak_kd), colour = AXIS, size = 1.8) +
    geom_point(data = filter(k, !weak_kd), aes(colour = cell_line), size = 1.8) +
    geom_text_repel(data = lab, aes(label = label), size = pt_to_size(FS$annot), colour = INK2, seed = 1,
                    min.segment.length = 0, segment.colour = AXIS, box.padding = 0.5,
                    nudge_y = ifelse(lab$kd_efficiency > 0.8, 0, 0.08), nudge_x = ifelse(lab$kd_efficiency > 0.8, -30, 0),
                    ylim = c(0, 1)) +
    facet_wrap(~ strip, scales = "free_x") +
    scale_colour_manual(values = CELL_COL, guide = "none") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = sprintf("%d of %d knockdowns remove at least half of the target RNA",
                         sum(!k$weak_kd), nrow(k)),
         subtitle = sprintf("Knockdown efficiency = 1 − 2^(DESeq2 log2 fold change of the target gene); median %s. Grey: excluded (dashed line = %s)",
                            pct(median(k$kd_efficiency)), pct(d$params$MIN_KD)),
         x = "Knockdowns ranked by efficiency", y = "Target RNA removed") +
    theme_viz() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_blank())
}

# ---- c: does a stronger knockdown change more genes? ---------------------------------------------------------
fig_response <- function(d) {
  k <- d$kd |> mutate(cell_line = factor(cell_line, names(CELL_COL)), label = paste(target, cell_line))
  rho <- k |> group_by(cell_line) |> summarise(r = cor(kd_efficiency, n_de, method = "spearman"),
                                               p = cor.test(kd_efficiency, n_de, method = "spearman", exact = FALSE)$p.value)
  rs <- paste(sprintf("%s ρ = %.2f (p %s)", rho$cell_line, rho$r, ifelse(rho$p < 0.001, "< 0.001", sprintf("= %.2f", rho$p))), collapse = "; ")
  lab <- k |> group_by(cell_line) |> slice_max(n_de, n = FIG$response$n_label) |> ungroup()
  ggplot(k, aes(kd_efficiency, n_de)) +
    geom_vline(xintercept = d$params$MIN_KD, colour = INK2, linewidth = 0.4, linetype = "22") +
    geom_point(aes(colour = cell_line, shape = weak_kd, size = n_kd_same_control), stroke = 0.7, alpha = 0.85) +
    geom_text_repel(data = lab, aes(label = label), size = pt_to_size(FS$annot), colour = INK2, seed = 1,
                    min.segment.length = 0, segment.colour = AXIS, nudge_y = 0.35, direction = "x") +
    scale_colour_manual(values = CELL_COL, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = "Kept", `TRUE` = "Excluded (< 50%)"), name = NULL) +
    scale_size_area(max_size = 4, breaks = c(1, 4, 8), name = "Knockdowns\nsharing control") +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_log10(labels = scales::label_comma(), expand = expansion(mult = c(0.04, 0.12))) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2.5)), shape = guide_legend(order = 2, override.aes = list(size = 2.5))) +
    labs(title = if (all(rho$p > 0.05)) "Knockdown strength does not predict how many genes change"
                 else "Stronger knockdowns change more genes",
         subtitle = sprintf("Differentially expressed genes (DESeq2 padj < 0.05) per knockdown, all %d knockdowns
Spearman: %s",
                            nrow(k), rs),
         x = "Target RNA removed", y = "DE genes (log scale)") +
    theme_viz() + theme(legend.box = "horizontal")
}

# ---- d: how many genes does each knockdown's eCLIP bind? ----------------------------------------------------
fig_eclip <- function(d) {
  k <- d$kd |> filter(in_binding_set) |> arrange(n_bound) |> mutate(rank = row_number(), label = paste(target, cell_line),
         status = factor(STATUS[excluded_reason], STATUS), group = factor(group, names(GROUP_COL)))
  n_st <- table(k$status)
  st_lab <- sprintf("%s (%d)", levels(k$status), as.integer(n_st))
  tl <- d$params$tier_lab; th <- d$params$TIERS
  tiers <- tibble(y = c(th[1], th[1], th[2]), vj = c(1.6, -0.6, -0.6), lab = tl,
                  n = as.integer(table(factor(k$tier[k$analysed], tl))))
  re <- k |> filter(reassigned) |> mutate(txt = paste0(label, ": reassigned to ", CATEGORY_LAB[as.character(group)]))
  ggplot(k, aes(rank, n_bound)) +
    geom_hline(yintercept = th, colour = INK2, linewidth = 0.4, linetype = "22") +
    geom_text(data = tiers, aes(x = nrow(k), y = y, label = sprintf("%s genes · %d analysed", lab, n), vjust = vj), hjust = 1,
              size = pt_to_size(FS$annot), colour = INK2) +
    geom_point(aes(colour = group, shape = status), size = 2, stroke = 0.8) +
    geom_text_repel(data = re, aes(label = txt), size = pt_to_size(FS$annot), colour = INK2, seed = 1,
                    min.segment.length = 0, segment.colour = AXIS, nudge_y = -0.85, nudge_x = 25, hjust = 0) +
    scale_colour_manual(values = GROUP_COL, labels = GROUP_LAB, name = NULL, drop = TRUE) +
    scale_shape_manual(values = setNames(STATUS_SHAPE, STATUS), labels = st_lab, name = NULL) +
    scale_y_log10(labels = scales::label_comma(), breaks = c(30, 100, 300, 1000, 3000)) +
    guides(shape = guide_legend(order = 1, ncol = 1), colour = guide_legend(order = 2, ncol = 2, override.aes = list(shape = 16))) +
    labs(title = sprintf("%d knockdowns enter the analysis; the ncRNA binders sit at the low end", sum(k$analysed)),
         subtitle = sprintf("Bound genes = genes with ≥ 1 eCLIP peak among the genes tested after knockdown (expressed, nuclear TPM ≥ 1); median %s",
                            scales::comma(median(k$n_bound))),
         caption = sprintf("%d knockdowns with < 30 bound genes are not shown (excluded in Aim 2). ncRNA binders: Aim 1 category 'ncRNA'; an ncRNA-category RBP with ≥ %d bound mRNA genes\nand < %d%% of peaks on ncRNA is moved to its second category.",
                           sum(!d$kd$in_binding_set), d$params$REASSIGN_GENES, d$params$REASSIGN_NC_PCT),
         x = "Knockdowns ranked by bound genes", y = "Bound genes (log scale)") +
    theme_viz() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_blank(),
                        legend.position = "right", legend.justification = "top")
}

# ---- e: analysis-set flow chart ----------------------------------------------------------------------------
fig_flow <- function(d) {
  f <- d$flow; k <- d$kd |> filter(analysed); tl <- d$params$tier_lab
  final <- paste(vapply(tl, function(t) sprintf("%s bound genes: %d (HepG2 %d, K562 %d)", t, sum(k$tier == t),
                                                sum(k$tier == t & k$cell_line == "HepG2"), sum(k$tier == t & k$cell_line == "K562")), ""),
                 collapse = "\n")
  function() {
    grid.newpage(); grid.rect(gp = gpar(fill = SURFACE, col = NA))
    pushViewport(viewport(x = 0.5, y = 0.5, width = unit(1, "npc") - unit(8, "mm"), height = unit(1, "npc") - unit(8, "mm")))
    grid.text(sprintf("From %d knockdowns to %d analysed", f$n[1], f$n[nrow(f)]), x = 0, y = 1, just = c("left", "top"),
              gp = gpar(fontsize = FS$title, fontface = "bold", col = INK))
    grid.text("Each box: knockdowns remaining. Right: knockdowns removed at that step and why.", x = 0, y = unit(1, "npc") - unit(7, "mm"),
              just = c("left", "top"), gp = gpar(fontsize = FS$subtitle, col = INK2))
    n <- nrow(f); top <- 0.84; bot <- 0.20; ys <- seq(top, bot, length.out = n); bw <- 0.42; bh <- 0.075; bx <- 0.24
    for (i in seq_len(n)) {
      last <- i == n
      grid.roundrect(x = bx, y = ys[i], width = bw, height = bh, r = unit(1.5, "mm"),
                     gp = gpar(fill = if (last) HIGHLIGHT else MID, col = NA))
      grid.text(sprintf("%d  %s", f$n[i], if (last) "analysed" else f$step[i]), x = bx, y = ys[i],
                gp = gpar(fontsize = FS$axis_text, col = if (last) SURFACE else INK, fontface = if (last) "bold" else "plain"))
      if (i > 1) {
        grid.lines(x = c(bx, bx), y = c(ys[i - 1] - bh / 2, ys[i] + bh / 2), arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
                   gp = gpar(col = AXIS, fill = AXIS, lwd = 1))
        ym <- (ys[i - 1] + ys[i]) / 2
        grid.lines(x = c(bx, bx + bw / 2 + 0.03), y = c(ym, ym), gp = gpar(col = AXIS, lwd = 0.8, lty = "22"))
        grid.text(sprintf("−%d  %s", f$dropped[i], f$dropped_reason[i]), x = bx + bw / 2 + 0.04, y = ym, just = "left",
                  gp = gpar(fontsize = FS$annot, col = INK2))
      }
    }
    grid.text(final, x = bx - bw / 2, y = ys[n] - bh / 2 - 0.02, just = c("left", "top"),
              gp = gpar(fontsize = FS$annot, col = INK2, lineheight = 1.15))
    popViewport()
  }
}

# ---- draw everything (order: all samples -> knockdown strength -> response -> filtered eCLIP set -> flow) ----------
d <- readRDS("data/analysis/kd_qc_plotdata.rds")
save_fig(fig_samples(d), "kd_qc_a_samples", FIG$samples$w, FIG$samples$h)
save_fig(fig_efficiency(d), "kd_qc_b_efficiency", FIG$efficiency$w, FIG$efficiency$h)
save_fig(fig_response(d), "kd_qc_c_response", FIG$response$w, FIG$response$h)
save_fig(fig_eclip(d), "kd_qc_d_eclip_genes", FIG$eclip$w, FIG$eclip$h)
save_drawing(fig_flow(d), "kd_qc_e_flow", FIG$flow$w, FIG$flow$h)
message("Knockdown QC figures written to ", FIG_DIR)
