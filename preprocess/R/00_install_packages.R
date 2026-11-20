# Installs the R packages for the ENCODE eCLIP atlas project (Windows binaries only).
# Usage: Rscript preprocess/R/00_install_packages.R

options(repos = c(CRAN = "https://cloud.r-project.org"))
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("Installing into:", lib, "\n")

cran <- c("BiocManager", "dplyr", "tidyr", "readr", "purrr", "stringr", "ggplot2",
          "ggrepel", "patchwork", "uwot", "arrow", "renv")
bioc <- c("GenomicRanges", "rtracklayer", "txdbmaker", "GenomicFeatures", "DESeq2",
          "apeglm", "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db")

# Never compile code (no Rtools); safe to re-run: only installs what is missing.
options(install.packages.compile.from.source = "never")
missing <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]

# Annotation *data* packages ship as source only but contain no compiled code.
data_pkgs <- c("GO.db", "org.Hs.eg.db", "GenomeInfoDbData")

m <- missing(cran)
if (length(m)) install.packages(m, lib = lib, type = "binary")
m <- missing(data_pkgs)
if (length(m)) BiocManager::install(m, lib = lib, ask = FALSE, update = FALSE, type = "source")
m <- missing(setdiff(bioc, data_pkgs))
if (length(m)) BiocManager::install(m, lib = lib, ask = FALSE, update = FALSE, type = "binary")

# Report results
all_pkgs <- c(cran, bioc)
ok <- vapply(all_pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
cat("\n==== Install summary ====\n")
for (p in all_pkgs) cat(sprintf("%-16s %s\n", p, if (ok[[p]]) "OK" else "FAILED"))
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
if (!all(ok)) quit(status = 1)
