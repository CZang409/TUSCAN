# ============================================================
# GeoMx DSP 空间转录组数据分析 + 细胞类型 Deconvolution
# ------------------------------------------------------------
# 场景假设:
#   1. 作者提供的矩阵已做过 normalization + batch effect correction
#      (这类矩阵常见于ComBat等批次矫正方法输出,可能包含负值/非线性
#       变换,不适合直接喂给count-based的deconvolution模型)
#   2. 你手头还有原始 raw counts 矩阵 (基因 x ROI/AOI)
#   3. 你有自己的 scRNA-seq 数据集,将作为deconvolution的参考
#      (用于构建 cell profile / signature matrix)
#
# 核心deconvolution工具: SpatialDecon (Danaher & Kim, Nat Commun 2022)
#   - 专为 NanoString GeoMx 数据设计的log-normal回归反卷积算法
#   - 内置基于negative probe的background估计
#   - 支持用自定义scRNA-seq衍生的signature matrix替代内置safeTME
#
# 请搜索 "## >>> 修改这里" 找到需要根据你实际数据调整的地方
# ============================================================


## ============ 0. 安装与加载R包 ============
pkgs_cran <- c("dplyr", "tidyr", "ggplot2", "pheatmap", "Seurat", "matrixStats")
pkgs_bioc <- c("SpatialDecon", "GeomxTools", "NanoStringNCTools", "Biobase")

installed <- rownames(installed.packages())
for (p in pkgs_cran) if (!p %in% installed) install.packages(p, dependencies = TRUE)
if (!"BiocManager" %in% installed) install.packages("BiocManager")
for (p in pkgs_bioc) if (!p %in% installed) BiocManager::install(p, update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(SpatialDecon)
  library(Seurat)
})

set.seed(20260729)

## ============ 1. 路径与关键参数设置 ============
## >>> 修改这里: 换成你自己的文件路径和列名
path_raw_counts    <- "data/raw_counts_matrix.csv"        # 原始counts: 基因(含NegProbe) x ROI/AOI
path_norm_matrix   <- "data/author_normalized_matrix.csv" # 作者提供的normalized+batch-corrected矩阵
path_segment_meta  <- "data/segment_annotation.csv"       # ROI/AOI注释(slide/scan/segment/group等)
path_scref_rds     <- "data/scRNA_reference.rds"          # 你自己的scRNA-seq参考数据(Seurat对象)
scref_celltype_col <- "cell_type"                         # scRNA-seq meta.data中细胞类型所在列名
group_col          <- "group"                             # segment注释中用于分组比较的列名


## ============ 2. 读取GeoMx数据 ============

# 2.1 原始counts (基因 x ROI),行=gene/probe(含negative probe),列=ROI/AOI segment ID
raw_counts <- as.matrix(read.csv(path_raw_counts, row.names = 1, check.names = FALSE))

# 2.2 作者提供的 normalized & batch-corrected 矩阵(基因 x ROI),仅用于下游常规分析
norm_author <- as.matrix(read.csv(path_norm_matrix, row.names = 1, check.names = FALSE))

# 2.3 ROI/AOI层面的注释,行名需与样本ID一致
segment_meta <- read.csv(path_segment_meta, row.names = 1, check.names = FALSE)

# 2.4 样本对齐:三份数据的样本ID必须完全对应
common_samples <- Reduce(intersect, list(colnames(raw_counts), colnames(norm_author), rownames(segment_meta)))
cat(sprintf("三份数据共有 %d 个共同的ROI/AOI样本\n", length(common_samples)))
stopifnot(length(common_samples) > 0)

raw_counts   <- raw_counts[, common_samples, drop = FALSE]
norm_author  <- norm_author[, common_samples, drop = FALSE]
segment_meta <- segment_meta[common_samples, , drop = FALSE]

# 2.5 识别negative probe行(通常行名含"NegProbe")
neg_idx <- grepl("negprobe", rownames(raw_counts), ignore.case = TRUE)
has_negprobe <- any(neg_idx)
cat(sprintf("检测到 %d 个negative probe行\n", sum(neg_idx)))


## ============ 3. 为什么deconvolution要单独重新normalize ============
## SpatialDecon基于log-normal回归模型,要求输入是线性尺度、非负的标准化矩阵
## (推荐Q3 normalization),而不是log值,也不是batch校正后可能出现负值/非线性
## 变换的矩阵(如ComBat常引入负值)。因此:
##   - norm_author  -> 用于聚类/差异表达/PCA等下游分析(已处理批次效应)
##   - 用raw counts重新做的Q3 normalization -> 专门喂给SpatialDecon做deconvolution
if (any(norm_author < 0, na.rm = TRUE)) {
  message("提示: 作者提供的normalized矩阵中存在负值(常见于批次矫正后的输出),",
          "不建议直接用于SpatialDecon。脚本将用raw counts重新做Q3 normalization。")
}

q3_normalize <- function(counts_mat) {
  q3 <- apply(counts_mat, 2, function(x) quantile(x[x > 0], 0.75, na.rm = TRUE))
  q3_geomean <- exp(mean(log(q3[q3 > 0])))
  norm_factors <- q3 / q3_geomean
  list(norm = sweep(counts_mat, 2, norm_factors, "/"), factors = norm_factors)
}

raw_genes  <- raw_counts[!neg_idx, , drop = FALSE]
q3_out     <- q3_normalize(raw_genes)
norm_q3    <- q3_out$norm


## ============ 4. Background 估计(SpatialDecon核心输入之一) ============
if (has_negprobe) {
  neg_raw  <- raw_counts[neg_idx, , drop = FALSE]
  neg_norm <- sweep(neg_raw, 2, q3_out$factors, "/")   # 用同一套归一化因子缩放阴性探针

  norm_for_bg <- rbind(norm_q3, neg_norm)
  probepool <- rep(1, nrow(norm_for_bg))
  ## >>> 修改这里: 若你的数据混合了多个probe pool(如CTA+Panel Plus),
  ## 需要按 fData 中的 Module 信息为每个基因分别指定 pool 编号,而不是全部设为1

  bg <- derive_GeoMx_background(norm = norm_for_bg, probepool = probepool, negnames = rownames(neg_norm))
  bg <- bg[rownames(norm_q3), , drop = FALSE]
} else {
  message("未检测到negative probe原始counts,使用简化背景近似(低分位数估计)。",
          "建议尽量补充阴性探针数据以提高deconvolution准确性。")
  bg <- matrix(rep(apply(norm_q3, 2, quantile, probs = 0.1), each = nrow(norm_q3)),
               nrow = nrow(norm_q3), dimnames = dimnames(norm_q3))
}


## ============ 5. 用你自己的scRNA-seq数据构建cell profile matrix ============
scref <- readRDS(path_scref_rds)   ## >>> 修改这里: 假设是Seurat对象
stopifnot(scref_celltype_col %in% colnames(scref@meta.data))

# 5.1 提取raw counts (gene x cell)。create_profile_matrix需要未标准化的counts
sc_counts <- tryCatch(
  as.matrix(GetAssayData(scref, layer = "counts")),   # Seurat v5写法
  error = function(e) as.matrix(GetAssayData(scref, slot = "counts"))  # Seurat v4写法
)

# 5.2 (可选但推荐) 统一基因symbol大小写,避免GeoMx与scRNA-seq基因名不匹配
# rownames(sc_counts) <- toupper(rownames(sc_counts))
# rownames(norm_q3)   <- toupper(rownames(norm_q3))
# rownames(bg)        <- toupper(rownames(bg))

# 5.3 细胞注释表(需要 CellID + cellType 两列)
cell_annots <- data.frame(
  CellID   = colnames(sc_counts),
  cellType = as.character(scref@meta.data[[scref_celltype_col]]),
  stringsAsFactors = FALSE
)
print(table(cell_annots$cellType))
## >>> 建议在这里检查是否有需要手动剔除的低质量/不明确细胞类型标签

# 5.4 生成cell profile matrix(基因 x 细胞类型 的signature matrix)
custom_profile_mtx <- create_profile_matrix(
  mtx              = sc_counts,
  cellAnnots       = cell_annots,
  cellTypeCol      = "cellType",
  cellNameCol      = "CellID",
  matrixName       = "custom_scRNA_profile",
  outDir           = NULL,     # 填路径可将矩阵另存为.RData
  normalize        = FALSE,    # 输入已是raw counts
  minCellNum       = 20,       # 每种细胞类型所需最少细胞数   ## >>> 按数据量调整
  minGenes         = 10,
  scalingFactor    = 5,        # 官方vignette推荐值,与GeoMx Q3-normalized量级匹配
  discardCellTypes = TRUE      # 自动剔除doublet/mitotic/unknown等标签
)
cat("最终纳入的细胞类型:\n"); print(colnames(custom_profile_mtx))


## ============ 6. 运行 SpatialDecon ============
shared_genes <- Reduce(intersect, list(rownames(norm_q3), rownames(bg), rownames(custom_profile_mtx)))
cat(sprintf("GeoMx数据与scRNA-seq参考矩阵共有 %d 个基因可用于deconvolution\n", length(shared_genes)))
if (length(shared_genes) < 100) {
  warning("共有基因数偏少,请检查基因symbol命名是否一致(如大小写、Ensembl vs Symbol)")
}

decon_res <- spatialdecon(
  norm        = norm_q3,
  bg          = bg,
  X           = custom_profile_mtx,
  raw         = raw_genes,     # 提供raw counts用于估计基因权重,降低低counts基因影响
  align_genes = TRUE
)

## decon_res主要输出:
##   $beta          细胞丰度得分 (细胞类型 x ROI),未强制归一化,可跨样本比较相对量级
##   $prop_of_all   每个ROI内各细胞类型占比(按行加总=1),适合展示组成比例
##   $p / $t        每个细胞类型-ROI估计的显著性


## ============ 7. 整理并保存结果 ============
cell_fraction <- as.data.frame(t(decon_res$prop_of_all))
cell_fraction$SampleID <- rownames(cell_fraction)

cell_fraction_annot <- cell_fraction %>%
  left_join(segment_meta %>% mutate(SampleID = rownames(segment_meta)), by = "SampleID")

write.csv(cell_fraction_annot, "GeoMx_deconvolution_cell_proportions.csv", row.names = FALSE)
cat("细胞组成比例已保存至 GeoMx_deconvolution_cell_proportions.csv\n")


## ============ 8. 可视化 ============

# 8.1 堆叠柱状图:每个ROI的细胞类型组成
plot_df <- cell_fraction %>%
  pivot_longer(-SampleID, names_to = "CellType", values_to = "Proportion")

p_stack <- ggplot(plot_df, aes(x = SampleID, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  labs(x = "ROI / AOI", y = "细胞类型占比", fill = "Cell type",
       title = "GeoMx ROI 细胞类型组成 (SpatialDecon)")
ggsave("deconvolution_stacked_barplot.pdf", p_stack, width = 12, height = 6)

# 8.2 热图:细胞丰度(beta)在样本间的分布,按分组列做列注释
if (group_col %in% colnames(segment_meta)) {
  annotation_col <- segment_meta[colnames(decon_res$beta), group_col, drop = FALSE]
  pheatmap(decon_res$beta, scale = "row", annotation_col = annotation_col,
           show_colnames = FALSE, main = "Cell abundance (beta) heatmap",
           filename = "deconvolution_beta_heatmap.pdf")
}

# 8.3 按分组比较某个细胞类型的比例(示例取第一个细胞类型,可替换为你关心的类型)
if (group_col %in% colnames(cell_fraction_annot)) {
  target_cell <- colnames(decon_res$prop_of_all)[1]   ## >>> 修改这里: 换成你关心的细胞类型名
  p_box <- ggplot(cell_fraction_annot, aes(x = .data[[group_col]], y = .data[[target_cell]],
                                            fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    theme_minimal(base_size = 12) +
    labs(y = paste0(target_cell, " proportion"), x = "Group",
         title = paste0(target_cell, " abundance by group"))
  ggsave(paste0("boxplot_", gsub("[^A-Za-z0-9]", "_", target_cell), "_by_group.pdf"),
         p_box, width = 6, height = 5)

  # 简单的组间统计检验(两组用wilcox.test,多组用kruskal.test)
  n_groups <- length(unique(cell_fraction_annot[[group_col]]))
  if (n_groups == 2) {
    test_res <- wilcox.test(cell_fraction_annot[[target_cell]] ~ cell_fraction_annot[[group_col]])
  } else if (n_groups > 2) {
    test_res <- kruskal.test(cell_fraction_annot[[target_cell]] ~ cell_fraction_annot[[group_col]])
  }
  print(test_res)
}


## ============ 9. (可选) 基于作者normalized矩阵的常规下游分析 ============
## 聚类/差异表达/PCA等建议基于 norm_author (已处理批次效应),
## deconvolution则依赖上面单独重新计算的 norm_q3。
pca_res <- prcomp(t(log2(pmax(norm_author, 0) + 1)), scale. = FALSE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$SampleID <- rownames(pca_df)
pca_df <- left_join(pca_df, cell_fraction_annot, by = "SampleID")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA (author-normalized matrix)")
ggsave("PCA_author_normalized.pdf", p_pca, width = 6, height = 5)

cat("全部分析完成。\n")
