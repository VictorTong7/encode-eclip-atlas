# Shared helpers for the binding-vs-expression analyses (analysis/R/13, 18, 19).
#
# load_kd_data()         per analysed knockdown, the gene universe with binding labels (13, 18, 19).
#   universe = genes tested by DESeq2 for that knockdown (data/analysis/de_results.parquet, non-missing shrunken LFC)
#              AND expressed in that cell line (nuclear TPM >= 1; data/analysis/expression_nuclear.parquet),
#              minus the knocked-down RBP's own gene (it changes by design).
#   bound    = >= 1 peak of the knockdown's eCLIP experiment (kd_groups.tsv: the RBP x cell line experiment
#              with most peaks, as in analysis/R/11).
# prep_kd_genes()        shared setup of the figure set (13, 18, 19): genes with padj, changed flag, expression deciles.
# label_hub_categories() number of other RBPs binding each gene, hub x binding category (19; 18 hub check).

suppressPackageStartupMessages({ library(dplyr) })

load_kd_data <- function(kd) {
  sel <- read.delim("analysis/results/tables/kd_selected.tsv") |> dplyr::select(kd_experiment, target_gene_id)
  kd <- kd |> left_join(sel, by = "kd_experiment")
  expr <- arrow::read_parquet("data/analysis/expression_nuclear.parquet") |> filter(expressed) |>
    dplyr::select(cell_line, gene_id, ntpm = median_tpm)
  de <- arrow::read_parquet("data/analysis/de_results.parquet",
                            col_select = c("kd_experiment", "gene_id", "lfc_shrunk", "padj")) |>
    filter(kd_experiment %in% kd$kd_experiment, !is.na(lfc_shrunk))
  bound <- arrow::read_parquet("data/processed/peaks_annotated.parquet", col_select = c("experiment", "gene_id")) |>
    filter(experiment %in% kd$eclip_experiment, !is.na(gene_id)) |> distinct(experiment, gene_id)
  lapply(split(kd, kd$kd_experiment), function(k) {
    d <- de |> filter(kd_experiment == k$kd_experiment) |>
      inner_join(expr |> filter(cell_line == k$cell_line), by = "gene_id") |>
      filter(gene_id != k$target_gene_id) |>
      mutate(de = !is.na(padj) & padj < 0.05, abs_lfc = abs(lfc_shrunk),
             bound = gene_id %in% bound$gene_id[bound$experiment == k$eclip_experiment])
    list(info = k, d = d)
  })
}

# Shared setup for the sequential figure set (13_bound_vs_unbound.R figures a, d; 18_change_rate.R figures b, c):
# the knockdown's universe minus genes with missing padj (status unknown), changed = padj < 0.05, and N_DEC deciles of
# log1p control nuclear TPM within the knockdown. Every analysis in the set uses exactly these genes and deciles.
prep_kd_genes <- function(d, n_dec = 10) {
  d |> filter(!is.na(padj)) |> mutate(changed = padj < 0.05, dec = ntile(log1p(ntpm), n_dec))
}

# hub x binding categories per knockdown (19_hub_attribution.R; hub check in 18_change_rate.R).
# n_other = RBPs binding the gene other than the knocked-down RBP (gene_rbp_counts.tsv minus this RBP's own binding).
label_hub_categories <- function(kdl, hub_hi = 20) {
  nr <- read.delim("analysis/results/tables/gene_rbp_counts.tsv") |> dplyr::select(cell_line, gene_id, n_rbps, gene_span)
  lapply(kdl, function(x) {
    x$d <- x$d |> inner_join(nr |> filter(cell_line == x$info$cell_line) |> dplyr::select(-cell_line), by = "gene_id") |>
      mutate(n_other = n_rbps - bound,
             category = case_when(n_other >= hub_hi & !bound ~ "hub only", n_other >= hub_hi & bound ~ "both",
                                  bound ~ "specific", TRUE ~ "neither"))
    x
  })
}
