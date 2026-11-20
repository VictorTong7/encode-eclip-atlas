# Build the transcript-region annotation (GENCODE V19, hg19) and per-cell-line expression.
#
# Inputs : data/ref/gencode.v19.annotation.gtf.gz, data/ref/hg19.fa.fai,
#          data/processed/manifest.parquet + data/raw/ctrl_counts/*.tsv
# Outputs: data/ref/gencode_v19.sqlite              TxDb
#          data/processed/annot_features.rds        unflattened region features (for peak assignment)
#          data/processed/annot_flat.rds            per cell line: flattened, priority-resolved regions
#                                                   of expressed genes (for region lengths)
#          data/processed/expression_by_cell.parquet  median control TPM per gene x cell line
#          preprocess/results/tables/region_lengths.tsv
#
# Transcript set: protein-coding genes -> protein_coding transcripts tagged "basic";
# other genes -> "basic" transcripts; small ncRNA genes -> all transcripts; pseudogenes excluded.
#
# Usage: Rscript preprocess/R/03_build_annotation.R   (from the project root)

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomicFeatures); library(dplyr)
})
source("preprocess/R/lib_regions.R")
set.seed(1)
gtf_path <- "data/ref/gencode.v19.annotation.gtf.gz"
dir.create("preprocess/results/tables", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Transcript metadata (parsed directly so repeated "tag" attributes are kept) --------
message("Reading transcript metadata ...")
lines <- readLines(gzfile(gtf_path))
lines <- lines[grepl("\ttranscript\t", lines, fixed = TRUE)]
attr_of <- function(key) sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", lines)
tx <- data.frame(tx_id = attr_of("transcript_id"), gene_id = attr_of("gene_id"),
                 gene_name = attr_of("gene_name"), gene_type = attr_of("gene_type"),
                 transcript_type = attr_of("transcript_type"),
                 basic = grepl('tag "basic"', lines, fixed = TRUE),
                 chrom = sub("\t.*", "", lines))
rm(lines)
tx <- tx |>
  filter(chrom %in% CANONICAL, !grepl("pseudogene", gene_type), !grepl("pseudogene", transcript_type)) |>
  filter(case_when(gene_type %in% SMALL_NCRNA ~ TRUE,
                   gene_type == "protein_coding" ~ basic & transcript_type == "protein_coding",
                   TRUE ~ basic))
message(sprintf("  transcript set: %d transcripts, %d genes", nrow(tx), n_distinct(tx$gene_id)))

# ---- 2. TxDb and region features ------------------------------------------------------------
db_path <- "data/ref/gencode_v19.sqlite"
if (!file.exists(db_path)) {
  message("Building TxDb (several minutes) ...")
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf",
                                                     dataSource = "GENCODE V19", organism = "Homo sapiens",
                                                     taxonomyId = 9606))
  AnnotationDbi::saveDb(txdb, db_path)
}
txdb <- AnnotationDbi::loadDb(db_path)

grab <- function(grl, region) {
  grl <- grl[names(grl) %in% tx$tx_id]
  gr <- unlist(grl, use.names = TRUE)
  gr <- gr[as.character(seqnames(gr)) %in% CANONICAL]
  mcols(gr) <- NULL
  gr$tx_id <- names(gr); names(gr) <- NULL
  gr$region <- region
  gr
}
message("Extracting regions ...")
coding_tx <- tx$tx_id[tx$transcript_type == "protein_coding"]
nc_tx <- tx$tx_id[tx$transcript_type != "protein_coding"]
exons <- exonsBy(txdb, "tx", use.names = TRUE)
nc_exons <- grab(exons[names(exons) %in% nc_tx], "noncoding_exon")
nc_exons$region[nc_exons$tx_id %in% tx$tx_id[tx$gene_type %in% SMALL_NCRNA]] <- "small_ncRNA"
introns <- grab(intronsByTranscript(txdb, use.names = TRUE), "intron")
w <- width(introns)
prox <- c(resize(introns, pmin(w, PROXIMAL_NT), fix = "start"),
          resize(introns[w > PROXIMAL_NT], pmin(w[w > PROXIMAL_NT] - PROXIMAL_NT, PROXIMAL_NT), fix = "end"))
prox$region <- "proximal_intron"
dist <- introns[w > 2 * PROXIMAL_NT]
dist <- narrow(dist, start = PROXIMAL_NT + 1L, end = -(PROXIMAL_NT + 1L))
dist$region <- "distal_intron"
feat <- c(grab(cdsBy(txdb, "tx", use.names = TRUE), "CDS"),
          grab(threeUTRsByTranscript(txdb, use.names = TRUE), "3UTR"),
          grab(fiveUTRsByTranscript(txdb, use.names = TRUE), "5UTR"),
          nc_exons, prox, dist)
m <- match(feat$tx_id, tx$tx_id)
feat$gene_id <- tx$gene_id[m]; feat$gene_name <- tx$gene_name[m]; feat$gene_type <- tx$gene_type[m]
feat$subtype <- ifelse(feat$region == "small_ncRNA", feat$gene_type, NA_character_)
feat$priority <- match(feat$region, REGIONS)
fai <- read.delim("data/ref/hg19.fa.fai", header = FALSE)
seqlevels(feat, pruning.mode = "coarse") <- CANONICAL
seqlengths(feat) <- fai$V2[match(CANONICAL, fai$V1)]
# drop identical features shared by several transcripts of the same gene (keep distinct genes/regions)
key <- paste(seqnames(feat), start(feat), end(feat), strand(feat), feat$region, feat$gene_id)
feat <- feat[!duplicated(key)]
saveRDS(feat, "data/processed/annot_features.rds")
message(sprintf("  %d region features", length(feat)))
print(table(feat$region)[REGIONS])

# ---- 3. Expression per cell line from control RNA-seq -------------------------------------
message("Computing control expression ...")
man <- arrow::read_parquet("data/processed/manifest.parquet") |> filter(category == "ctrl_counts")
tpm <- lapply(seq_len(nrow(man)), function(i) {
  d <- read.delim(man$local_path[i], colClasses = c(gene_id = "character"))[, c("gene_id", "TPM")]
  d <- d[startsWith(d$gene_id, "ENSG"), ]
  data.frame(gene_id = d$gene_id, TPM = d$TPM, cell_line = man$cell_line[i], file = man$accession[i])
}) |> bind_rows()
expr <- tpm |> group_by(gene_id, cell_line) |>
  summarise(median_tpm = median(TPM), n_samples = n(), .groups = "drop") |>
  mutate(expressed = median_tpm >= EXPRESSED_TPM)
# Small ncRNAs (<200 nt: snRNA, snoRNA, miRNA, ...) are barely measured by standard RNA-seq,
# so TPM cannot tell whether they are present. Count them as expressed in both cell lines.
small_genes <- unique(tx$gene_id[tx$gene_type %in% SMALL_NCRNA])
expr <- expr |> mutate(small_ncRNA = gene_id %in% small_genes,
                       expressed = expressed | small_ncRNA)
arrow::write_parquet(expr, "data/processed/expression_by_cell.parquet")
print(expr |> group_by(cell_line) |> summarise(samples = max(n_samples), expressed_genes = sum(expressed)))

# ---- 4. Flattened regions of expressed genes, per cell line -----------------------------------
flatten <- function(cell) {
  e <- expr[expr$cell_line == cell, ]
  ex_genes <- e$gene_id[e$expressed]
  f <- feat[feat$gene_id %in% ex_genes]
  tpm_f <- e$median_tpm[match(f$gene_id, e$gene_id)]
  dj <- disjoin(f, with.revmap = TRUE)          # strand-aware disjoint pieces
  hits <- S4Vectors::Hits(rep(seq_along(dj), lengths(dj$revmap)), unlist(dj$revmap),
                          length(dj), length(f), sort.by.query = TRUE)
  best <- pick_best(hits, f$priority, tpm_f)
  out <- dj[best$query]
  out$revmap <- NULL
  out$region <- f$region[best$feature]; out$gene_id <- f$gene_id[best$feature]
  out$tpm <- tpm_f[best$feature]
  # abundance-weighted background: floor small-ncRNA TPM at 1 (not measurable by RNA-seq)
  out$tpm[out$region == "small_ncRNA"] <- pmax(out$tpm[out$region == "small_ncRNA"], 1)
  out$cell_line <- cell
  out
}
message("Flattening regions per cell line ...")
flat <- list(K562 = flatten("K562"), HepG2 = flatten("HepG2"))
stopifnot(all(sapply(flat, function(g) isDisjoint(g))))
saveRDS(flat, "data/processed/annot_flat.rds")

lens <- bind_rows(lapply(flat, function(g)
  data.frame(cell_line = g$cell_line, region = g$region, width = width(g), tpm = g$tpm))) |>
  group_by(cell_line, region) |>
  summarise(length_nt = sum(as.numeric(width)), length_x_tpm = sum(as.numeric(width) * tpm), .groups = "drop") |>
  mutate(region = factor(region, REGIONS)) |> arrange(cell_line, region)
write.table(lens, "preprocess/results/tables/region_lengths.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(lens |> mutate(length_Mb = round(length_nt / 1e6, 2)) |> select(cell_line, region, length_Mb)))
message("Done.")
