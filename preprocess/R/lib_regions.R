# Shared definitions for transcript-region annotation (sourced by R/03 and R/04).

REGIONS <- c("CDS", "3UTR", "5UTR", "small_ncRNA", "noncoding_exon",
             "proximal_intron", "distal_intron")            # priority order, highest first
REGION_LEVELS <- c(REGIONS, "intergenic")
SMALL_NCRNA <- c("snRNA", "snoRNA", "miRNA", "misc_RNA", "rRNA", "Mt_rRNA", "Mt_tRNA")
CANONICAL <- paste0("chr", c(1:22, "X", "Y", "M"))
PROXIMAL_NT <- 500L
EXPRESSED_TPM <- 1

# For each query in `hits` (a Hits object: query -> feature), keep the single best feature:
# highest region priority, then highest expression of the feature's gene.
# Returns a data.frame(query, feature) with one row per query that had any hit.
pick_best <- function(hits, feat_priority, feat_tpm) {
  q <- S4Vectors::queryHits(hits)
  f <- S4Vectors::subjectHits(hits)
  tpm <- feat_tpm[f]
  tpm[is.na(tpm)] <- -1
  o <- order(q, feat_priority[f], -tpm)
  q <- q[o]; f <- f[o]
  keep <- !duplicated(q)
  data.frame(query = q[keep], feature = f[keep])
}

# Transcript landmarks (1-nt, strand-aware) for the given transcripts: 5' and 3' splice sites
# (first / last intron base) and start / stop codons (first / last CDS base).
build_landmarks <- function(txdb, tx_keep) {
  intr <- unlist(GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)[tx_keep])
  cds_r <- unlist(range(GenomicFeatures::cdsBy(txdb, "tx", use.names = TRUE)))
  cds_r <- cds_r[names(cds_r) %in% tx_keep]
  one <- function(x, fix) {
    y <- unique(GenomicRanges::granges(GenomicRanges::resize(x, 1L, fix = fix)))
    y[as.character(GenomicRanges::seqnames(y)) %in% CANONICAL]
  }
  list(`5′ splice site` = one(intr, "start"), `3′ splice site` = one(intr, "end"),
       `Start codon` = one(cds_r, "start"), `Stop codon` = one(cds_r, "end"))
}
