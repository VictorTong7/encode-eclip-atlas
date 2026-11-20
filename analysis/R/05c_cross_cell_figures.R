# Figure 2: is RBP binding intrinsic? Draws figures from data/analysis/cross_cell.rds (R/05b).
#
# Outputs: analysis/results/figures/fig2a_sharing.{pdf,png}, fig2b_difference_decomposition,
#          fig2c_binding_maps, fig2d_group_features, fig2e_example_targets
# Usage: Rscript analysis/R/05c_cross_cell_figures.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })
source("preprocess/R/lib_regions.R")
cc <- readRDS("data/analysis/cross_cell.rds")

# palette: dataviz reference instance, light mode; categorical slots in fixed order
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
SURFACE <- "#fcfcfb"; S <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
theme_viz <- function(base = 9) {
  theme_minimal(base_size = base, base_family = "sans") +
    theme(plot.background = element_rect(fill = SURFACE, colour = NA),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = AXIS, linewidth = 0.3),
          text = element_text(colour = INK2), axis.text = element_text(colour = MUTED),
          strip.text = element_text(colour = INK, face = "bold", hjust = 0),
          plot.title = element_text(colour = INK, face = "bold", size = base + 2),
          plot.subtitle = element_text(colour = INK2), legend.position = "top")
}
save_fig <- function(p, name, w, h) for (ext in c("pdf", "png"))
  ggsave(file.path("analysis/results/figures", paste0(name, ".", ext)), p, width = w, height = h, dpi = 200, bg = SURFACE)

# ---- 2a: sharing at three levels, test vs ceiling vs nulls ------------------------------------------
TYPE_LAB <- c(null1 = "Different RBPs", null2 = "Different RBPs,\nsame binding class",
              test = "Same RBP,\nK562 vs HepG2", ceiling = "Same RBP, two K562\nexperiments (ceiling)")
MET_LAB <- c(site_matched = "Binding sites shared\n(depth-matched)", gene_frac = "Target genes shared",
             map_r = "Positional map correlation\n(splice sites, start/stop codons)")
long <- cc$lv12 |> select(type, target, all_of(names(MET_LAB))) |>
  pivot_longer(all_of(names(MET_LAB)), names_to = "metric", values_to = "v") |> filter(!is.na(v)) |>
  mutate(type = factor(type, names(TYPE_LAB), TYPE_LAB), metric = factor(metric, names(MET_LAB), MET_LAB))
nulls <- filter(long, grepl("Different", type)); pts <- filter(long, !grepl("Different", type))
p2a <- ggplot(long, aes(v, type)) +
  geom_boxplot(data = nulls, fill = GRID, colour = MUTED, width = 0.5, outlier.size = 0.4,
               outlier.colour = AXIS, linewidth = 0.3) +
  geom_jitter(data = pts, aes(fill = type), shape = 21, colour = SURFACE, stroke = 0.4, size = 2.2,
              position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  stat_summary(data = pts, fun = median, geom = "crossbar", width = 0.6, colour = INK, linewidth = 0.3) +
  scale_fill_manual(values = setNames(S[1:2], TYPE_LAB[c("test", "ceiling")]), guide = "none") +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  labs(title = "The same RBP binds the same places in both cell lines",
       subtitle = "Agreement between experiments; bars mark medians, boxes summarise all pairs of different RBPs",
       x = NULL, y = NULL) +
  theme_viz() + theme(panel.grid.major.y = element_blank(), panel.spacing = unit(1.2, "lines"))
save_fig(p2a, "fig2a_sharing", 10, 4)

# ---- 2b: why are some peaks found in only one line? ---------------------------------------------------
CLS <- c("shared", "not expressed in other line", "expressed ≥2× lower", "weak peak", "unexplained")
CLS_LAB <- c("Found in both", "Gene not expressed in other", "Gene ≥2× lower in other", "Weak peak", "Unexplained")
dec <- cc$decomp |> mutate(cls = factor(as.character(cls), rev(CLS), rev(CLS_LAB)))
pooled <- dec |> group_by(type, cls) |> summarise(n = sum(n), .groups = "drop") |>
  group_by(type) |> mutate(frac = n / sum(n)) |> ungroup() |>
  mutate(row = factor(ifelse(type == "test", "K562 vs HepG2\n(70 RBPs)", "Two K562 experiments\n(6 RBPs, ceiling)"),
                      c("Two K562 experiments\n(6 RBPs, ceiling)", "K562 vs HepG2\n(70 RBPs)")),
         lab = ifelse(frac >= 0.06, sprintf("%.0f%%", 100 * frac), ""),
         lab_col = ifelse(cls == "Weak peak", INK, "#ffffff"))
fills <- setNames(S, CLS_LAB)
p2b1 <- ggplot(pooled, aes(frac, row, fill = cls)) +
  geom_col(width = 0.6, colour = SURFACE, linewidth = 0.6) +
  geom_text(aes(label = lab, colour = lab_col, group = cls), position = position_stack(vjust = 0.5), size = 2.8) +
  scale_colour_identity() + scale_fill_manual(values = fills, breaks = CLS_LAB, name = NULL) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(title = "Most cell-line differences are explained by expression or are within experimental noise",
       subtitle = "Fate of each peak in the other experiment (pooled over both directions)", x = NULL, y = NULL) +
  theme_viz() + theme(panel.grid.major.y = element_blank()) + guides(fill = guide_legend(nrow = 1))
per <- dec |> filter(type == "test") |> group_by(target) |> mutate(un = frac[cls == "Unexplained"]) |> ungroup() |>
  mutate(target = reorder(target, -un))
p2b2 <- ggplot(per, aes(frac, target, fill = cls)) +
  geom_col(width = 0.8, colour = SURFACE, linewidth = 0.25) +
  scale_fill_manual(values = fills, breaks = CLS_LAB, guide = "none") +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(subtitle = "Per RBP, ordered by unexplained share", x = "Share of peaks", y = NULL) +
  theme_viz(8) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 5.5))
save_fig(p2b1 / p2b2 + plot_layout(heights = c(1, 6)), "fig2b_difference_decomposition", 8, 11)

# ---- 2c: positional maps for example RBPs -------------------------------------------------------------
nb <- length(cc$bins) - 1L
centers <- head(cc$bins, -1) + diff(cc$bins) / 2
map_df <- function(target) {
  ex <- cc$pk[cc$pk$target == target, ]
  bind_rows(lapply(c("K562", "HepG2"), function(cell) {
    v <- cc$maps[, ex[[cell]]]
    data.frame(target = target, cell_line = cell, landmark = rep(cc$landmarks, each = nb),
               pos = rep(centers, length(cc$landmarks)), density = 100 * v)
  }))
}
ex_rbps <- c("RBFOX2", "U2AF2", "UPF1")
md <- bind_rows(lapply(ex_rbps, map_df)) |>
  mutate(landmark = factor(landmark, cc$landmarks), target = factor(target, ex_rbps),
         cell_line = factor(cell_line, c("K562", "HepG2")))
rlab <- cc$lv12 |> filter(type == "test", target %in% ex_rbps) |>
  transmute(target = factor(target, ex_rbps), lab = sprintf("r = %.2f", map_r),
            landmark = factor(cc$landmarks[1], cc$landmarks))
p2c <- ggplot(md, aes(pos, density, colour = cell_line)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_text(data = rlab, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = 1.1, vjust = 1.5,
            size = 2.6, colour = INK2) +
  scale_colour_manual(values = c(K562 = S[1], HepG2 = S[2]), name = NULL) +
  facet_grid(target ~ landmark, scales = "free_y") +
  scale_x_continuous(breaks = c(-300, 0, 300)) +
  labs(title = "Binding maps around transcript landmarks are reproduced across cell lines",
       subtitle = "Peak midpoints per 10-nt bin, % of the RBP's peaks; r = correlation of the full map (all four landmarks)",
       x = "Position relative to landmark (nt; positive = downstream)", y = "% of peaks") +
  theme_viz(8) + theme(strip.text.y = element_text(angle = 0))
save_fig(p2c, "fig2c_binding_maps", 9, 5.5)

# ---- 2d: what distinguishes conserved from not-conserved RBPs ------------------------------------------
FEAT_LAB <- c(log10_min_peaks = "Peaks in weaker experiment (log10)", abs_log2_peak_ratio = "Peak-count imbalance",
              lifted_2022_batch = "Includes a 2022-batch experiment", same_release_year = "Both released same year",
              frac_peaks_shared_expressed = "Peaks on genes expressed in both", rbp_log2_tpm = "RBP's own expression",
              rbp_abs_log2_ratio = "RBP expression difference", target_expr_divergence = "Target expression divergence",
              intronic_binder = "Intronic binder", utr3_or_cds_binder = "3′UTR / CDS binder",
              nucleus = "GO: nucleus", cytoplasm = "GO: cytoplasm", rna_splicing = "GO: RNA splicing",
              translation = "GO: translation", mrna_stability = "GO: mRNA stability",
              ribosome_biogenesis = "GO: ribosome biogenesis")
ef <- cc$eff |> mutate(label = FEAT_LAB[feature], sig = ifelse(p_bh < 0.05, "BH-adjusted p < 0.05", "not significant"),
                       label = reorder(label, estimate))
p2d <- ggplot(ef, aes(estimate, label)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", colour = MUTED, linewidth = 0.5) +
  geom_point(aes(fill = sig), shape = 21, colour = SURFACE, stroke = 0.5, size = 2.6) +
  scale_fill_manual(values = c(`BH-adjusted p < 0.05` = S[1], `not significant` = MUTED), guide = "none") +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title = "No feature cleanly separates conserved from less-conserved RBPs",
       subtitle = sprintf(paste0("Conserved (n = %d) vs not conserved (n = %d); Cliff's delta or difference in proportion, 95%% bootstrap CI.
",
                                 "None survives BH correction. Experiment depth is the strongest trend
",
                                 "and correlates with site sharing (Spearman rho = %.2f, BH p = %.3f)"),
                          ef$conserved_n[1], ef$not_conserved_n[1],
                          ef$spearman_with_site_sharing[ef$feature == "log10_min_peaks"],
                          ef$spearman_p_bh[ef$feature == "log10_min_peaks"]),
       x = "Effect size (positive = higher in conserved RBPs)", y = NULL) +
  theme_viz() + theme(panel.grid.major.y = element_blank())
save_fig(p2d, "fig2d_group_features", 8, 5)

# ---- 2e: per-gene binding for one conserved and one less-conserved RBP ------------------------------------
gene_scatter <- function(target) {
  ex <- cc$pk[cc$pk$target == target, ]
  sumg <- function(e) { d <- cc$gr_targets[[e]]; d <- d[d$gene_id %in% cc$universe, ]; tapply(d$signal, d$gene_id, sum) }
  a <- sumg(ex$K562); b <- sumg(ex$HepG2); u <- union(names(a), names(b))
  data.frame(target = target, K562 = ifelse(is.na(a[u]), 0, a[u]), HepG2 = ifelse(is.na(b[u]), 0, b[u]))
}
grp <- cc$groups
pick_c <- "QKI"; pick_n <- "HNRNPC"
sc <- bind_rows(gene_scatter(pick_c), gene_scatter(pick_n)) |>
  mutate(panel = factor(sprintf("%s (%s)", target, grp$group[match(target, grp$target)]),
                        sprintf("%s (%s)", c(pick_c, pick_n), grp$group[match(c(pick_c, pick_n), grp$target)])),
         status = case_when(K562 > 0 & HepG2 > 0 ~ "Bound in both", TRUE ~ "Bound in one line"))
p2e <- ggplot(sc, aes(log2(K562 + 1), log2(HepG2 + 1))) +
  geom_point(aes(fill = status), shape = 21, colour = SURFACE, stroke = 0.2, size = 1.6, alpha = 0.85) +
  scale_fill_manual(values = c(`Bound in both` = S[1], `Bound in one line` = AXIS), name = NULL) +
  facet_wrap(~ panel) + coord_equal() +
  labs(title = "Target genes: a conserved RBP vs a less-conserved one",
       subtitle = "Summed peak log2 enrichment per gene (genes expressed in both lines); 0 = not bound",
       x = "K562, log2(1 + summed log2FC)", y = "HepG2, log2(1 + summed log2FC)") +
  theme_viz()
save_fig(p2e, "fig2e_example_targets", 8, 4.6)
message("Figures written to analysis/results/figures/fig2*")
