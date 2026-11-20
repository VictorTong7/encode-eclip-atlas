# Lift the 27 GRCh38-only eCLIP IDR peak sets over to hg19.
#
# Input : data/processed/manifest.parquet (rows with needs_liftover == TRUE)
#         data/raw/eclip_peaks_grch38/*.bed.gz
#         data/ref/hg38ToHg19.over.chain   (unzipped by preprocess/python/02_download.py)
# Output: data/raw/eclip_peaks_hg19/<accession>_lifted.bed.gz   (narrowPeak, hg19)
#         preprocess/results/tables/liftover_stats.tsv
#
# A peak is kept only if it maps 1:1 (a single interval) onto a canonical chromosome
# (chr1-22, X, Y, M), on the same strand, with width within 10% of the original.
#
# Usage: Rscript preprocess/R/02b_liftover_peaks.R   (run from the project root)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
})

root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "data/processed/manifest.parquet")))
canonical <- paste0("chr", c(1:22, "X", "Y", "M"))
np_cols <- c("chrom", "start", "end", "name", "score", "strand",
             "signalValue", "pValue", "qValue", "peak")

man <- arrow::read_parquet(file.path(root, "data/processed/manifest.parquet")) |>
  filter(needs_liftover)
chain <- import.chain(file.path(root, "data/ref/hg38ToHg19.over.chain"))
out_dir <- file.path(root, "data/raw/eclip_peaks_hg19")
dir.create(file.path(root, "preprocess/results/tables"), recursive = TRUE, showWarnings = FALSE)

lift_one <- function(i) {
  rec <- man[i, ]
  bed <- read.delim(file.path(root, rec$local_path), header = FALSE,
                    col.names = np_cols, colClasses = c("character", "integer", "integer",
                    "character", "integer", "character", "numeric", "numeric", "numeric", "integer"))
  gr <- GRanges(bed$chrom, IRanges(bed$start + 1L, bed$end), strand = bed$strand)
  mcols(gr) <- bed[, c("name", "score", "signalValue", "pValue", "qValue", "peak")]
  gr$idx <- seq_along(gr)

  lifted <- liftOver(gr, chain)
  n_pieces <- lengths(lifted)
  one <- n_pieces == 1L
  hit <- unlist(lifted[one])
  orig <- gr[one]
  keep <- as.character(seqnames(hit)) %in% canonical &
    as.character(strand(hit)) == as.character(strand(orig)) &
    abs(width(hit) - width(orig)) <= 0.1 * width(orig)
  hit <- sort(hit[keep])

  out <- data.frame(chrom = as.character(seqnames(hit)), start = start(hit) - 1L, end = end(hit),
                    name = hit$name, score = hit$score, strand = as.character(strand(hit)),
                    signalValue = hit$signalValue, pValue = hit$pValue, qValue = hit$qValue,
                    peak = hit$peak)
  dest <- file.path(out_dir, paste0(rec$accession, "_lifted.bed.gz"))
  con <- gzfile(dest, "w")
  write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(con)

  tibble(accession = rec$accession, experiment = rec$experiment, target = rec$target,
         cell_line = rec$cell_line, peaks_in = length(gr), unmapped = sum(n_pieces == 0L),
         split = sum(n_pieces > 1L), filtered = sum(!keep), peaks_out = length(hit),
         pct_dropped = round(100 * (1 - length(hit) / length(gr)), 2),
         output = sub(paste0(root, "/"), "", dest, fixed = TRUE))
}

stats <- bind_rows(lapply(seq_len(nrow(man)), lift_one))
write.table(stats, file.path(root, "preprocess/results/tables/liftover_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Lifted %d experiments: %d -> %d peaks; median dropped %.2f%%, max %.2f%% (%s)\n",
            nrow(stats), sum(stats$peaks_in), sum(stats$peaks_out),
            median(stats$pct_dropped), max(stats$pct_dropped),
            stats$target[which.max(stats$pct_dropped)]))
flag <- filter(stats, pct_dropped > 5)
if (nrow(flag) > 0) {
  cat("WARNING: experiments losing >5% of peaks:\n")
  print(as.data.frame(flag[, c("accession", "target", "cell_line", "pct_dropped")]))
}
