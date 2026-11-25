# Aim 2, step 1: differential expression for each RBP knockdown vs its matched control.
#
# Selection: shRNA knockdowns only; one per RBP x cell line (strongest target knockdown by TPM).
# Model    : one DESeq2 model per shared control batch (control + every selected knockdown that uses
#            it; design ~ condition, control = reference); each knockdown is its own coefficient,
#            shrunk with apeglm.
# Genes    : RSEM expected_count (rounded), GENCODE V19 genes only (tRNA rows dropped), >= 10 counts
#            in >= 2 samples of the batch.
#
# Outputs: data/analysis/de_results.parquet       knockdown x gene results
#          data/analysis/common_response.parquet  per cell line x gene: fraction of knockdowns DE
#          analysis/results/tables/kd_selected.tsv, kd_summary.tsv
#          data/analysis/de_cache/<control>.rds   per-batch cache (re-runs resume; --recompute)
# Usage  : Rscript analysis/R/10_knockdown_de.R [--only=N] [--workers=8] [--recompute]

suppressPackageStartupMessages({ library(dplyr); library(BiocParallel) })
args <- commandArgs(TRUE)
arg_val <- function(name, default) {
  a <- args[startsWith(args, paste0("--", name, "="))]
  if (length(a)) sub(paste0("--", name, "="), "", a[1]) else default
}
ONLY <- as.integer(arg_val("only", NA)); WORKERS <- as.integer(arg_val("workers", 8))
RECOMPUTE <- "--recompute" %in% args
CACHE <- "data/analysis/de_cache"; dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

pairs <- arrow::read_parquet("data/processed/kd_pairs.parquet")
man <- arrow::read_parquet("data/processed/manifest.parquet")
feat <- readRDS("data/processed/annot_features.rds")
expr <- arrow::read_parquet("data/processed/expression_by_cell.parquet")

read_counts <- function(path) {
  d <- read.delim(path, colClasses = c(gene_id = "character"))[, c("gene_id", "expected_count", "TPM")]
  d[startsWith(d$gene_id, "ENSG"), ]
}
files_of <- function(exp, cat) man$local_path[man$category == cat & man$experiment == exp]

# ---- 1. target gene ids and TPM-based knockdown efficiency; selection ------------------------------
gene_names <- unique(data.frame(gene_name = feat$gene_name, gene_id = feat$gene_id))
target_gene <- function(sym, cell) {
  ids <- gene_names$gene_id[gene_names$gene_name == sym]
  if (length(ids) <= 1) return(if (length(ids)) ids else NA_character_)
  e <- expr[expr$cell_line == cell & expr$gene_id %in% ids, ]
  e$gene_id[which.max(e$median_tpm)]
}
pairs$target_gene_id <- mapply(target_gene, pairs$target, pairs$cell_line)
tpm_of <- function(exp, cat, gid) median(sapply(files_of(exp, cat), function(f) {
  d <- read_counts(f); sum(d$TPM[d$gene_id == gid]) }))
pairs$tpm_ratio <- mapply(function(k, c, g) (tpm_of(k, "kd_counts", g) + 0.5) / (tpm_of(c, "ctrl_counts", g) + 0.5),
                          pairs$kd_experiment, pairs$control_experiment, pairs$target_gene_id)
pairs <- pairs |> group_by(target, cell_line) |>
  mutate(selected = assay_title == "shRNA RNA-seq" &
           tpm_ratio == min(tpm_ratio[assay_title == "shRNA RNA-seq"])) |> ungroup() |>
  mutate(reason = case_when(selected ~ "selected", assay_title != "shRNA RNA-seq" ~ "CRISPRi (shRNA only)",
                            TRUE ~ "weaker duplicate shRNA knockdown"))
write.table(pairs |> select(kd_experiment, control_experiment, target, cell_line, assay_title,
                            target_gene_id, tpm_ratio, selected, reason) |> mutate(tpm_ratio = round(tpm_ratio, 4)),
            "analysis/results/tables/kd_selected.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
sel <- pairs |> filter(selected)
message(sprintf("Selected %d shRNA knockdowns (%d RBP x cell line) in %d control batches; dropped %d",
                nrow(sel), nrow(distinct(sel, target, cell_line)), n_distinct(sel$control_experiment), sum(!pairs$selected)))

# ---- 2-3. one DESeq2 model per control batch ------------------------------------------------------------
batches <- split(sel, sel$control_experiment)
if (!is.na(ONLY)) batches <- head(batches, ONLY)
run_batch <- function(b, man, cache, recompute) {
  suppressPackageStartupMessages({ library(DESeq2); library(dplyr) })
  ctl <- b$control_experiment[1]
  out_f <- file.path(cache, paste0(ctl, ".rds"))
  if (file.exists(out_f) && !recompute) return(readRDS(out_f))
  rc <- function(path) { d <- read.delim(path, colClasses = c(gene_id = "character"))[, c("gene_id", "expected_count")]
                         d[startsWith(d$gene_id, "ENSG"), ] }
  samp <- rbind(data.frame(file = man$local_path[man$category == "ctrl_counts" & man$experiment == ctl], condition = "control"),
                do.call(rbind, lapply(b$kd_experiment, function(k)
                  data.frame(file = man$local_path[man$category == "kd_counts" & man$experiment == k], condition = k))))
  mats <- lapply(samp$file, rc)
  genes <- mats[[1]]$gene_id
  stopifnot(all(sapply(mats, function(m) identical(m$gene_id, genes))))
  cnt <- round(sapply(mats, `[[`, "expected_count")); rownames(cnt) <- genes
  cnt <- cnt[rowSums(cnt >= 10) >= 2, , drop = FALSE]
  cd <- data.frame(condition = factor(make.names(samp$condition), levels = c("control", make.names(b$kd_experiment))))
  dds <- DESeqDataSetFromMatrix(cnt, cd, ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- bind_rows(lapply(seq_len(nrow(b)), function(i) {
    lev <- make.names(b$kd_experiment[i]); coef <- paste0("condition_", lev, "_vs_control")
    r <- results(dds, name = coef)
    s <- suppressMessages(lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE))
    data.frame(kd_experiment = b$kd_experiment[i], target = b$target[i], cell_line = b$cell_line[i],
               control_experiment = ctl, gene_id = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange,
               lfcSE = r$lfcSE, lfc_shrunk = s$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
  }))
  saveRDS(res, out_f)
  res
}
message("Running DESeq2 on ", length(batches), " control batches with ", WORKERS, " workers ...")
t0 <- Sys.time()
bp <- if (WORKERS > 1) SnowParam(workers = WORKERS, progressbar = FALSE) else SerialParam()
de <- bind_rows(bplapply(batches, run_batch, man = man, cache = CACHE, recompute = RECOMPUTE, BPPARAM = bp))
message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- 4-5. efficiency, QC, common response ------------------------------------------------------------------
summ <- de |> inner_join(sel |> select(kd_experiment, target_gene_id, tpm_ratio), by = "kd_experiment") |>
  group_by(kd_experiment, target, cell_line, control_experiment, tpm_ratio) |>
  summarise(n_genes = n(), n_de = sum(padj < 0.05, na.rm = TRUE),
            target_lfc = log2FC[gene_id == target_gene_id[1]][1],
            target_lfc_shrunk = lfc_shrunk[gene_id == target_gene_id[1]][1],
            target_lfc_percentile = mean(log2FC <= log2FC[gene_id == target_gene_id[1]][1], na.rm = TRUE),
            .groups = "drop") |>
  mutate(kd_efficiency = 1 - 2^target_lfc, weak_kd = kd_efficiency < 0.5)
# correlation of LFC with other knockdowns sharing the control (control-artefact check)
w <- de |> select(kd_experiment, control_experiment, gene_id, log2FC)
summ$mean_lfc_cor_same_control <- sapply(seq_len(nrow(summ)), function(i) {
  o <- w |> filter(control_experiment == summ$control_experiment[i], kd_experiment != summ$kd_experiment[i])
  if (!nrow(o)) return(NA_real_)
  me <- w |> filter(kd_experiment == summ$kd_experiment[i])
  mean(sapply(split(o, o$kd_experiment), function(x) {
    j <- match(me$gene_id, x$gene_id); cor(me$log2FC, x$log2FC[j], use = "complete.obs") }))
})
write.table(summ |> mutate(across(where(is.numeric), ~ round(.x, 4))), "analysis/results/tables/kd_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
common <- de |> group_by(cell_line, gene_id) |>
  summarise(n_kd = n(), frac_de = mean(padj < 0.05, na.rm = TRUE), .groups = "drop")
if (is.na(ONLY)) {
  arrow::write_parquet(de, "data/analysis/de_results.parquet")
  arrow::write_parquet(common, "data/analysis/common_response.parquet")
}

cat(sprintf("\n%d knockdowns: median DE genes %d (IQR %d-%d); median knockdown efficiency %.0f%%; <50%%: %d\n",
            nrow(summ), as.integer(median(summ$n_de)), as.integer(quantile(summ$n_de, .25)), as.integer(quantile(summ$n_de, .75)),
            100 * median(summ$kd_efficiency, na.rm = TRUE), sum(summ$weak_kd, na.rm = TRUE)))
cat(sprintf("DESeq2 vs TPM efficiency, Spearman: %.3f | target in bottom 1%% of LFCs: %d/%d\n",
            cor(summ$target_lfc, log2(summ$tpm_ratio), method = "spearman", use = "complete.obs"),
            sum(summ$target_lfc_percentile <= 0.01, na.rm = TRUE), nrow(summ)))
cat(sprintf("Mean LFC correlation with knockdowns sharing a control: median %.2f\n",
            median(summ$mean_lfc_cor_same_control, na.rm = TRUE)))
