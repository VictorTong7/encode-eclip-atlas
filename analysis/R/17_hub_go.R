# Binding vs expression: what are the hub genes? GO enrichment (COMPUTE).
#
# Hub genes = bound by >= 20 profiled RBPs in that cell line (hub_groups.tsv from 19_hub_attribution.R; universe =
# expressed genes tested in most knockdowns). Per cell line, GO enrichment (clusterProfiler enrichGO, org.Hs.eg.db,
# ENSEMBL ids without version; BH; q < 0.05; gene sets of 15-500 genes) against two backgrounds:
#   all       all genes of the universe (what distinguishes hub genes from expressed genes in general)
#   matched   hub genes + N_MATCH expression-matched non-hub genes per hub gene (20 quantile bins of log1p nuclear TPM):
#             what distinguishes hub genes from non-hub genes of the same expression. Main result, because hub genes are
#             ~10x more expressed than genes bound by < 5 RBPs and highly expressed genes are mostly housekeeping.
# Ontologies: BP (biological process) and CC (cellular component); redundant terms collapsed with clusterProfiler::simplify.
#
# Outputs: analysis/results/tables/hub_go.tsv, hub_gene_profile.tsv; data/analysis/hub_go_plotdata.rds
# Usage  : Rscript analysis/R/17_hub_go.R   (~3 min; needs 19_hub_attribution.R)

suppressPackageStartupMessages({ library(dplyr); library(clusterProfiler); library(org.Hs.eg.db) })
HUB <- "≥ 20"; N_MATCH <- 3; N_BIN <- 20; SEED <- 17

genes <- read.delim("analysis/results/tables/hub_groups.tsv") |>
  mutate(ens = sub("\\..*$", "", gene_id), hub = hub_group == HUB)
feat <- readRDS("data/processed/annot_features.rds")
gtype <- distinct(data.frame(gene_id = feat$gene_id, gene_type = feat$gene_type))
u3 <- feat[feat$region == "3UTR"]
utr <- sapply(split(u3, u3$gene_id), function(g) sum(GenomicRanges::width(GenomicRanges::reduce(g))))
genes <- genes |> left_join(gtype, by = "gene_id") |>
  left_join(data.frame(gene_id = names(utr), utr3_len = as.numeric(utr)), by = "gene_id") |>
  mutate(utr3_len = coalesce(utr3_len, 0))

set.seed(SEED)
matched_bg <- function(g) {
  br <- quantile(log1p(g$ntpm), seq(0, 1, length.out = N_BIN + 1))
  bin <- findInterval(log1p(g$ntpm), br[-c(1, N_BIN + 1)]) + 1
  need <- tabulate(bin[g$hub], N_BIN) * N_MATCH
  pick <- unlist(lapply(seq_len(N_BIN), function(b) { p <- which(!g$hub & bin == b)
    if (!length(p) || !need[b]) return(integer(0)); p[sample.int(length(p), min(need[b], length(p)))] }))
  g[pick, ]
}

go_one <- function(hub_ens, universe, ont) {
  e <- enrichGO(hub_ens, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = ont, universe = universe,
                pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 15, maxGSSize = 500, readable = TRUE)
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  e <- simplify(e, cutoff = 0.7, by = "p.adjust", select_fun = min)
  as.data.frame(e) |>
    mutate(k = as.numeric(sub("/.*", "", GeneRatio)), n = as.numeric(sub(".*/", "", GeneRatio)),
           K = as.numeric(sub("/.*", "", BgRatio)), N = as.numeric(sub(".*/", "", BgRatio)),
           fold_enrichment = (k / n) / (K / N))
}

res <- list(); prof <- list()
for (cl in c("K562", "HepG2")) {
  g <- genes |> filter(cell_line == cl)
  mb <- matched_bg(g)
  prof[[cl]] <- bind_rows(g |> filter(hub) |> mutate(set = "hub (≥ 20 RBPs)"),
                          mb |> mutate(set = "expression-matched non-hub"),
                          g |> filter(!hub) |> mutate(set = "all non-hub")) |>
    group_by(set) |>
    summarise(n = n(), median_tpm = median(ntpm), median_span = median(gene_span), median_utr3 = median(utr3_len),
              pct_protein_coding = 100 * mean(gene_type == "protein_coding", na.rm = TRUE), .groups = "drop") |>
    mutate(cell_line = cl, .before = 1)
  hub_ens <- g$ens[g$hub]
  for (bg in c("matched", "all")) {
    uni <- if (bg == "matched") unique(c(hub_ens, mb$ens)) else g$ens
    for (ont in c("BP", "CC")) {
      message(sprintf("%s · %s background · %s (%d hub genes, universe %d) ...", cl, bg, ont, length(hub_ens), length(uni)))
      r <- go_one(hub_ens, uni, ont)
      if (!is.null(r)) res[[paste(cl, bg, ont)]] <- r |> mutate(cell_line = cl, background = bg, ontology = ont, .before = 1)
    }
  }
}
go <- bind_rows(res)
prof <- bind_rows(prof)

# terms found in both cell lines (matched background)
both <- go |> filter(background == "matched") |> count(ontology, ID, Description) |> filter(n == 2)

write.table(go |> dplyr::select(cell_line, background, ontology, ID, Description, k, n, K, N, fold_enrichment, pvalue, p.adjust,
                                qvalue, geneID) |> mutate(across(c(fold_enrichment, pvalue, p.adjust, qvalue), ~ signif(.x, 4))),
            "analysis/results/tables/hub_go.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prof |> mutate(across(where(is.double), ~ signif(.x, 4))), "analysis/results/tables/hub_gene_profile.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
top_hubs <- genes |> filter(hub) |> group_by(gene_name) |> filter(n() == 2) |>
  summarise(frac = mean(frac_rbps), tpm = mean(ntpm), .groups = "drop") |> arrange(desc(frac))
saveRDS(list(go = go, prof = prof, both = both, top_hubs = top_hubs, params = list(N_MATCH = N_MATCH)),
        "data/analysis/hub_go_plotdata.rds")

# ---- summary ---------------------------------------------------------------------------------------------------------
cat("\nGene profiles:\n"); print(as.data.frame(prof), digits = 3)
cat("\nSignificant terms (after simplify):\n"); print(go |> count(cell_line, background, ontology) |> as.data.frame())
for (bg in c("matched", "all")) for (ont in c("BP", "CC")) {
  cat(sprintf("\n== %s background, %s: top terms ==\n", bg, ont))
  print(go |> filter(background == bg, ontology == ont) |> group_by(cell_line) |> slice_min(p.adjust, n = 12, with_ties = FALSE) |>
          ungroup() |> transmute(cell_line, Description = substr(Description, 1, 60), k, fold = round(fold_enrichment, 2),
                                 q = signif(p.adjust, 2)) |> as.data.frame())
}
cat(sprintf("\nTerms significant in both cell lines (matched background): %d\n", nrow(both)))
cat("\nTop 25 hub genes (both lines, by mean share of RBPs):\n"); print(head(top_hubs, 25) |> as.data.frame(), digits = 2)
message("Wrote hub_go.tsv, hub_gene_profile.tsv, hub_go_plotdata.rds")
