library(ggplot2)
library(ggrepel)

# 1. 先在完整 object 上算 centrality（不要 subset，避免 subsetCellChat 的 bug）
object.list <- lapply(object.list, function(x) {
  netAnalysis_computeCentrality(x, slot.name = "netP")
})

# 2. 提取每个 object 的 outgoing/incoming strength + interaction count
extractSignalingRoleDF <- function(object) {
  centr <- object@netP$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  num.link <- rowSums(object@net$count) + colSums(object@net$count) - diag(object@net$count)
  data.frame(
    celltype = rownames(outgoing),
    x = rowSums(outgoing),
    y = rowSums(incoming),
    Count = num.link
  )
}

# 3. 合并所有 condition 成一个长表，并打上 condition 标签
df.all <- do.call(rbind, lapply(names(object.list), function(nm) {
  df <- extractSignalingRoleDF(object.list[[nm]])
  df$condition <- nm
  df
}))

# 4. 只保留你要的 cell type
celltypes.use <- c("CellTypeA", "CellTypeB", "CellTypeC")
df.sel <- df.all[df.all$celltype %in% celltypes.use, ]

# 5. 点大小范围（用全体数据算，保证跨图/跨条件可比；如果只想按选中的算，用 df.sel 代替 df.all）
weight.MinMax <- c(min(df.all$Count), max(df.all$Count))

# 6. 画在同一张图上，颜色映射 condition
ggplot(df.sel, aes(x = x, y = y, colour = condition)) +
  geom_point(aes(size = Count), alpha = 0.75) +
  scale_size_continuous(range = c(3, 10), limits = weight.MinMax, name = "Strength") +
  geom_text_repel(aes(label = celltype), size = 3, max.overlaps = Inf,
                   show.legend = FALSE) +
  theme_bw() +
  labs(x = "Outgoing interaction strength",
       y = "Incoming interaction strength",
       colour = "Condition") +
  theme(legend.position = "right")

# ============================================================
# GSEA analysis for Fib2 vs Fib1
# Positive avg_log2FC / NES = enriched in Fib2
# Negative avg_log2FC / NES = enriched in Fib1
# ============================================================

# ----------------------------
# 0. Load packages
# ----------------------------
library(Seurat)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)

# ----------------------------
# 1. Prepare DEG table
# ----------------------------
# If gene names are stored as row names:
deg_fib2 <- deg_fib %>%
  rownames_to_column(var = "gene")

# Check required columns
stopifnot(
  all(c("gene", "avg_log2FC") %in% colnames(deg_fib2))
)

# Remove missing or infinite values
deg_fib2 <- deg_fib2 %>%
  filter(
    !is.na(gene),
    !is.na(avg_log2FC),
    is.finite(avg_log2FC)
  )

# ----------------------------
# 2. Build ranked gene list
# ----------------------------
# Do NOT filter by p-value or logFC for GSEA
rank_df <- deg_fib2 %>%
  select(gene, avg_log2FC) %>%
  distinct(gene, .keep_all = TRUE)

# Convert SYMBOL to ENTREZID
gene_map <- bitr(
  rank_df$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Join logFC with ENTREZ IDs
rank_entrez <- rank_df %>%
  inner_join(gene_map, by = c("gene" = "SYMBOL"))

# If multiple symbols map to the same ENTREZID,
# retain the value with the largest absolute logFC
rank_entrez <- rank_entrez %>%
  group_by(ENTREZID) %>%
  slice_max(
    order_by = abs(avg_log2FC),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Create named numeric vector
geneList <- rank_entrez$avg_log2FC
names(geneList) <- rank_entrez$ENTREZID

# GSEA requires decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# Inspect ranked genes
head(geneList)
tail(geneList)
length(geneList)

# ============================================================
# 3. Hallmark GSEA
# ============================================================

# Obtain Hallmark gene sets
hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  select(gs_name, ncbi_gene) %>%
  distinct()

# Run GSEA
gsea_hallmark <- GSEA(
  geneList = geneList,
  TERM2GENE = hallmark_df,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = FALSE,
  seed = TRUE
)

# View results
hallmark_result <- gsea_hallmark@result %>%
  arrange(p.adjust)

head(
  hallmark_result %>%
    select(
      ID,
      Description,
      setSize,
      enrichmentScore,
      NES,
      pvalue,
      p.adjust,
      qvalue,
      core_enrichment
    ),
  20
)

# Save Hallmark result
write.csv(
  hallmark_result,
  file = "Fib2_vs_Fib1_Hallmark_GSEA.csv",
  row.names = FALSE
)

# ----------------------------
# 4. Hallmark dot plot
# ----------------------------
p_hallmark <- dotplot(
  gsea_hallmark,
  showCategory = 20,
  split = ".sign"
) +
  facet_grid(. ~ .sign) +
  labs(
    title = "Hallmark GSEA: Fib2 vs Fib1",
    subtitle = "Positive NES = Fib2; Negative NES = Fib1"
  ) +
  theme_classic()

p_hallmark

ggsave(
  filename = "Fib2_vs_Fib1_Hallmark_GSEA_dotplot.pdf",
  plot = p_hallmark,
  width = 11,
  height = 7
)

# ----------------------------
# 5. Bar plot using NES
# ----------------------------
hallmark_plot_df <- hallmark_result %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    group = ifelse(NES > 0, "Fib2", "Fib1"),
    Description = gsub("^HALLMARK_", "", Description),
    Description = gsub("_", " ", Description)
  ) %>%
  arrange(NES) %>%
  slice_head(n = 10) %>%
  bind_rows(
    hallmark_result %>%
      filter(p.adjust < 0.05) %>%
      mutate(
        group = ifelse(NES > 0, "Fib2", "Fib1"),
        Description = gsub("^HALLMARK_", "", Description),
        Description = gsub("_", " ", Description)
      ) %>%
      arrange(desc(NES)) %>%
      slice_head(n = 10)
  ) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(
    Description = factor(
      Description,
      levels = Description[order(NES)]
    )
  )

p_nes <- ggplot(
  hallmark_plot_df,
  aes(x = NES, y = Description, fill = group)
) +
  geom_col() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Hallmark pathways enriched in Fib2 versus Fib1",
    subtitle = "Positive NES = Fib2; Negative NES = Fib1",
    x = "Normalized enrichment score",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 12)

p_nes

ggsave(
  filename = "Fib2_vs_Fib1_Hallmark_GSEA_NES.pdf",
  plot = p_nes,
  width = 9,
  height = 7
)

# ----------------------------
# 6. Plot selected Hallmark pathways
# ----------------------------
selected_hallmarks <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_APOPTOSIS"
)

available_hallmarks <- intersect(
  selected_hallmarks,
  gsea_hallmark@result$ID
)

available_hallmarks

# Example: plot EMT if available
if (
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" %in%
    gsea_hallmark@result$ID
) {
  
  p_emt <- gseaplot2(
    gsea_hallmark,
    geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    title = "Hallmark EMT: Fib2 vs Fib1",
    pvalue_table = TRUE
  )
  
  p_emt
  
  ggsave(
    filename = "Fib2_vs_Fib1_Hallmark_EMT_GSEA.pdf",
    plot = p_emt,
    width = 8,
    height = 6
  )
}

# Plot several selected pathways together
if (length(available_hallmarks) > 0) {
  
  p_selected <- gseaplot2(
    gsea_hallmark,
    geneSetID = available_hallmarks,
    title = "Selected Hallmark pathways: Fib2 vs Fib1",
    pvalue_table = TRUE
  )
  
  p_selected
}

# ============================================================
# 7. GO Biological Process GSEA
# ============================================================

gsea_go <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = FALSE,
  seed = TRUE
)

go_result <- gsea_go@result %>%
  arrange(p.adjust)

head(
  go_result %>%
    select(
      ID,
      Description,
      setSize,
      enrichmentScore,
      NES,
      pvalue,
      p.adjust,
      qvalue,
      core_enrichment
    ),
  20
)

write.csv(
  go_result,
  file = "Fib2_vs_Fib1_GO_BP_GSEA.csv",
  row.names = FALSE
)

# ----------------------------
# 8. GO GSEA dot plot
# ----------------------------
p_go <- dotplot(
  gsea_go,
  showCategory = 20,
  split = ".sign"
) +
  facet_grid(. ~ .sign) +
  labs(
    title = "GO Biological Process GSEA: Fib2 vs Fib1",
    subtitle = "Positive NES = Fib2; Negative NES = Fib1"
  ) +
  theme_classic()

p_go

ggsave(
  filename = "Fib2_vs_Fib1_GO_BP_GSEA_dotplot.pdf",
  plot = p_go,
  width = 12,
  height = 8
)

# ============================================================
# 9. Extract Fib2- and Fib1-enriched pathways
# ============================================================

hallmark_fib2 <- hallmark_result %>%
  filter(
    NES > 0,
    p.adjust < 0.05
  ) %>%
  arrange(desc(NES))

hallmark_fib1 <- hallmark_result %>%
  filter(
    NES < 0,
    p.adjust < 0.05
  ) %>%
  arrange(NES)

go_fib2 <- go_result %>%
  filter(
    NES > 0,
    p.adjust < 0.05
  ) %>%
  arrange(desc(NES))

go_fib1 <- go_result %>%
  filter(
    NES < 0,
    p.adjust < 0.05
  ) %>%
  arrange(NES)

# Inspect top pathways
head(
  hallmark_fib2 %>%
    select(Description, NES, p.adjust),
  15
)

head(
  hallmark_fib1 %>%
    select(Description, NES, p.adjust),
  15
)

head(
  go_fib2 %>%
    select(Description, NES, p.adjust),
  15
)

head(
  go_fib1 %>%
    select(Description, NES, p.adjust),
  15
)

# Save separate pathway tables
write.csv(
  hallmark_fib2,
  "Fib2_enriched_Hallmark_pathways.csv",
  row.names = FALSE
)

write.csv(
  hallmark_fib1,
  "Fib1_enriched_Hallmark_pathways.csv",
  row.names = FALSE
)

write.csv(
  go_fib2,
  "Fib2_enriched_GO_BP_pathways.csv",
  row.names = FALSE
)

write.csv(
  go_fib1,
  "Fib1_enriched_GO_BP_pathways.csv",
  row.names = FALSE
)
