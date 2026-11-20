# Consistency and biology checks for the region annotation (steps 03-04). Read-only.
# Prints PASS/FAIL per check; exits with status 1 if any check fails.
#
# Usage: Rscript preprocess/checks/check_annotation.R   (from the project root)

suppressPackageStartupMessages({ library(GenomicRanges); library(dplyr) })
source("preprocess/R/lib_regions.R")
n_fail <- 0L; n_all <- 0L
check <- function(name, ok, detail = "") {
  n_all <<- n_all + 1L; if (!isTRUE(ok)) n_fail <<- n_fail + 1L
  cat(sprintf("%s  %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, if (nzchar(detail)) paste0("  (", detail, ")") else ""))
}

feat <- readRDS("data/processed/annot_features.rds")
flat <- readRDS("data/processed/annot_flat.rds")
check("features: priority column matches REGIONS order", all(feat$priority == match(feat$region, REGIONS)))
check("features: canonical chromosomes only", all(as.character(seqnames(feat)) %in% CANONICAL))
check("flattened annotation is disjoint in each cell line", all(vapply(flat, isDisjoint, logical(1))),
      paste(names(flat), collapse = ", "))

# proximal + distal intron pieces exactly tile each transcript's introns
txdb <- AnnotationDbi::loadDb("data/ref/gencode_v19.sqlite")
tx_keep <- unique(feat$tx_id[feat$region %in% c("proximal_intron", "distal_intron")])
intr <- unlist(GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)[tx_keep])
intr <- intr[as.character(seqnames(intr)) %in% CANONICAL]
ints <- feat[feat$region %in% c("proximal_intron", "distal_intron")]
cov_i <- sum(as.numeric(width(reduce(intr)))); cov_f <- sum(as.numeric(width(reduce(ints))))
check("proximal + distal intron cover exactly the intronic bases", cov_i == cov_f,
      sprintf("%.1f Mb vs %.1f Mb", cov_i / 1e6, cov_f / 1e6))
prox <- feat[feat$region == "proximal_intron"]
check("proximal intron pieces are at most 500 nt wide", max(width(prox)) <= PROXIMAL_NT, sprintf("max %d", max(width(prox))))

# biology spot-checks: region ranks from log2_vs_all_rbps (peak share relative to the average RBP)
enr <- read.delim("preprocess/results/tables/region_enrichment.tsv")
rank_of <- function(rbp, region) {
  e <- enr |> filter(target == rbp) |> group_by(experiment) |>
    mutate(rk = rank(-log2_vs_all_rbps)) |> filter(region == !!region)
  if (nrow(e) == 0) return(NA_integer_)
  as.integer(min(e$rk))                                # best rank over that RBP's experiments
}
spot <- list(c("RBFOX2", "proximal_intron", 1), c("U2AF2", "proximal_intron", 2), c("PUM2", "3UTR", 2),
             c("DKC1", "small_ncRNA", 1), c("SRSF1", "CDS", 3))
for (s in spot) {
  rk <- rank_of(s[1], s[2])
  check(sprintf("biology: %s ranks %s in top %s", s[1], s[2], s[3]), !is.na(rk) && rk <= as.integer(s[3]),
        sprintf("best rank %s", rk))
}

cat(sprintf("\n%d/%d checks passed\n", n_all - n_fail, n_all))
quit(status = if (n_fail > 0) 1L else 0L)
