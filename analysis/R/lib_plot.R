# Shared plot styling for the analysis section. EDIT AESTHETICS HERE.
# Sourced by the *_figures.R scripts (Aim 1 so far; other aims keep their own copies until overhauled).
#
# Palette: dataviz reference instance, light mode. Categorical slots are used in fixed order
# (validated: adjacent-pair CVD dE >= 9.1, normal-vision >= 22.9; aqua and yellow are < 3:1 contrast
# on the surface, so groups drawn in those colours always carry a visible text label).

suppressPackageStartupMessages({ library(ggplot2) })

# ---- text, surfaces, grid ----------------------------------------------------------------------
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
SURFACE <- "#fcfcfb"; MID <- "#f0efec"
DIV_BLUE <- "#1c5cab"; DIV_RED <- "#c23b3b"              # diverging poles (neutral midpoint = MID)
SEQ_BLUE <- c("#e8f0fb", "#9fc0ec", "#2a78d6", "#0f3f80")  # sequential ramp (light -> dark)
HIGHLIGHT <- "#2a78d6"                                    # single-series accent
# ---- figure sizing: DRAW EVERY FIGURE AT ITS DISPLAY SIZE ---------------------------------------------
# On the website a figure is shown at a fixed CSS width. Drawing it at (display px / 96) inches makes
# 1 pt = 1.33 screen px in every figure, so the font sizes below look identical everywhere.
# Never rescale text by hand to fit a figure: change its height or simplify it instead.
PX_PER_IN <- 96
DISPLAY_PX <- c(text = 799, page = 894)    # measured 2026-11-19 at a 1400 px viewport (text column; .column-page-left, which leaves the right margin to the "On this page" menu)
FIG_W <- DISPLAY_PX / PX_PER_IN            # inches: text ~8.3, page ~11.8
# font scale in pt (x 1.33 = px on screen)
FS <- list(title = 12,       # ~16 px
           subtitle = 10,    # ~13 px
           axis_title = 10, axis_text = 10, legend = 10, strip = 10,
           annot = 8,        # data labels, RBP names, counts (~11 px)
           min = 7.5)        # dense labels only (~10 px)
pt_to_size <- function(pt) pt / ggplot2::.pt   # geom_text/geom_label size (mm) for a font size in pt
BASE_SIZE <- FS$axis_text
FONT <- "sans"

# ---- transcript regions (fixed order = annotation priority) ------------------------------------------
REGIONS <- c("CDS", "3UTR", "5UTR", "small_ncRNA", "noncoding_exon", "proximal_intron", "distal_intron")
REGION_COL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"), REGIONS)
REGION_LAB <- c(CDS = "CDS", `3UTR` = "3′UTR", `5UTR` = "5′UTR", small_ncRNA = "Small ncRNA",
                noncoding_exon = "Noncoding exon", proximal_intron = "Proximal intron",
                distal_intron = "Distal intron", intergenic = "Intergenic")

# ---- main binding category (Aim 1 grouping; each group shares its region's colour) -------------------
CATEGORIES <- c("CDS", "3UTR", "5UTR", "ncRNA", "splice_site", "deep_intron", "mixed")
CATEGORY_COL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#008300", "#6250d6", "#b4b3ad"),
                         CATEGORIES)
CATEGORY_LAB <- c(CDS = "CDS", `3UTR` = "3′UTR", `5UTR` = "5′UTR", ncRNA = "ncRNA",
                  splice_site = "Splice-site intron", deep_intron = "Deep intron", mixed = "No clear preference")

# ---- cell lines -------------------------------------------------------------------------------------------
CELL_COL <- c(K562 = "#2a78d6", HepG2 = "#eb6834")

# ---- theme and saving ------------------------------------------------------------------------------------
theme_viz <- function(base = BASE_SIZE) {
  theme_minimal(base_size = base, base_family = FONT) +
    theme(plot.background = element_rect(fill = SURFACE, colour = NA),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = AXIS, linewidth = 0.3),
          text = element_text(colour = INK2, size = base),
          axis.text = element_text(colour = MUTED, size = FS$axis_text),
          axis.title = element_text(size = FS$axis_title),
          plot.title = element_text(colour = INK, face = "bold", size = FS$title),
          plot.subtitle = element_text(colour = INK2, size = FS$subtitle),
          plot.caption = element_text(colour = MUTED, hjust = 0, size = FS$annot),
          strip.text = element_text(colour = INK, face = "bold", hjust = 0, size = FS$strip),
          legend.text = element_text(size = FS$legend), legend.title = element_text(size = FS$legend),
          plot.title.position = "plot", plot.caption.position = "plot",   # align titles with the figure edge
          plot.margin = margin(6, 12, 6, 6),
          legend.position = "top")
}

FIG_DIR <- "analysis/results/figures"
FIG_DPI <- 200

fig_width <- function(w) if (is.character(w)) unname(FIG_W[w]) else w   # "text" / "page" or inches

# ggplot / patchwork objects -> PDF + PNG. w = width class ("text" / "page"); h = height in inches
save_fig <- function(p, name, w, h) {
  w <- fig_width(w)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = SURFACE, device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = FIG_DPI, bg = SURFACE)
}

# grid-based drawings (ComplexHeatmap, UpSet): pass a function that draws the figure
save_drawing <- function(draw_fun, name, w, h) {
  w <- fig_width(w)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  cairo_pdf(file.path(FIG_DIR, paste0(name, ".pdf")), width = w, height = h); draw_fun(); dev.off()
  png(file.path(FIG_DIR, paste0(name, ".png")), width = w, height = h, units = "in", res = FIG_DPI); draw_fun(); dev.off()
}
