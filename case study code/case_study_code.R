setwd("./tuscan project/code")

library(tidyverse)
library(Seurat)
library(slingshot)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(phangorn)

################################################################################
######### calculate and visualize cnv burden ##########
load('breast_CNV_r3.RData')
breast_label <- read.csv("breast_labels.csv")

cnv.mat <- t(breast_cnv_obj@expr.data)
breast_label <- breast_label[match(rownames(cnv.mat),breast_label$barcode),]

cnv.mat.df <- as.data.frame(cnv.mat)
cnv.mat.df$label <- breast_label$clone3
cnv.mat.other <- cnv.mat.df %>% filter(label != 'Normal')

ref_idx = unlist(breast_cnv_obj@reference_grouped_cell_indices)
vals = breast_cnv_obj@expr.data[,ref_idx]
mean_ref_vals = mean(vals)
mean_ref_sd <- mean(apply(vals, 2, function(x) sd(x, na.rm=TRUE)))

para_up = mean_ref_vals + mean_ref_sd*2
para_low = mean_ref_vals - mean_ref_sd*2

cnv_status_other <- ifelse(cnv.mat.other > para_up, 1, 
                     ifelse(cnv.mat.other < para_low, -1, 0))

count_cnv.gain <- apply(cnv_status_other, 1, function(x) sum(x>0))
count_cnv.loss <- apply(cnv_status_other, 1, function(x) sum(x<0))

count.all <- data.frame(CNV_gain=count_cnv.gain, CNV_loss=count_cnv.loss, Group=factor(cnv.mat.other$label, 
                                                                                       levels = c('Tumor clone1','Tumor clone2','Tumor clone3')))

data_long <- pivot_longer(count.all, cols = c(CNV_gain, CNV_loss), names_to = "State", values_to = "Value")

p <- ggplot(data_long, aes(x = Group, y = Value, fill = State)) +
  geom_boxplot() +
  theme_minimal()+ 
  labs(y = "CNV burdens", title = "") +
  scale_fill_manual(values = c("CNV_gain" = "#D56365", "CNV_loss"="#048ABF"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
p

################################################################################
######### summarise cnv status ##########

cnv_matrix <- breast_cnv_obj@expr.data
cnv_status <- ifelse(cnv_matrix > para_up, 1, 
                     ifelse(cnv_matrix < para_low, -1, 0))

cnv_status.df <- as.data.frame(cnv_status)
cnv_status.df$gene <- row.names(cnv_status.df)

gene.band <- read.csv("gene_band_annotation.csv")
cnv_status.merge <- merge(cnv_status.df, gene.band, by='gene')
spot.names <- colnames(cnv_status.merge)[2:3799]

## average cnv status
cnv_by_band_avg <- cnv_status.merge %>% group_by(loc) %>% summarise(across(spot.names, mean, na.rm = TRUE))
band <- cnv_by_band_avg$loc
cnv_by_band_avg2 <- cnv_by_band_avg[,-1]
rownames(cnv_by_band_avg2) <- band
cnv_by_band_avg2 <- as.data.frame(t(cnv_by_band_avg2))

breast_label <- breast_label[match(row.names(cnv_by_band_avg2),breast_label$barcode),]
cnv_by_band_avg2$subclones <- breast_label$clone6

cnv_by_band_clone2 <- cnv_by_band_avg2 %>% group_by(subclones) %>% summarise(across(1:278, mean, na.rm = TRUE))
col.id <- cnv_by_band_clone2$subclones
cnv_by_band_clone2 <- cnv_by_band_clone2[,-1]
cnv_by_band_clone2 <- as.data.frame(t(cnv_by_band_clone2))
colnames(cnv_by_band_clone2) <- col.id

## summary cnv by band
rows_clone1_gain <- which(cnv_by_band_clone2$A > 0.2)
clone1_gain <- row.names(cnv_by_band_clone2[rows_clone1_gain, ])
rows_clone1_loss <- which(cnv_by_band_clone2$A < -0.2)
clone1_loss <- row.names(cnv_by_band_clone2[rows_clone1_loss, ])

rows_clone2_gain <- which(cnv_by_band_clone2$B > 0.2)
clone2_gain <- row.names(cnv_by_band_clone2[rows_clone2_gain, ])
rows_clone2_loss <- which(cnv_by_band_clone2$B < -0.2)
clone2_loss <- row.names(cnv_by_band_clone2[rows_clone2_loss, ])

rows_clone31_gain <- which(cnv_by_band_clone2$C > 0.2)
clone31_gain <- row.names(cnv_by_band_clone2[rows_clone31_gain, ])
rows_clone31_loss <- which(cnv_by_band_clone2$C < -0.2)
clone31_loss <- row.names(cnv_by_band_clone2[rows_clone31_loss, ])

rows_clone32_gain <- which(cnv_by_band_clone2$D > 0.2)
clone32_gain <- row.names(cnv_by_band_clone2[rows_clone32_gain, ])
rows_clone32_loss <- which(cnv_by_band_clone2$D < -0.2)
clone32_loss <- row.names(cnv_by_band_clone2[rows_clone32_loss, ])

rows_clone33_gain <- which(cnv_by_band_clone2$E > 0.2)
clone33_gain <- row.names(cnv_by_band_clone2[rows_clone33_gain, ])
rows_clone33_loss <- which(cnv_by_band_clone2$E < -0.2)
clone33_loss <- row.names(cnv_by_band_clone2[rows_clone33_loss, ])

rows_clone34_gain <- which(cnv_by_band_clone2$F > 0.2)
clone34_gain <- row.names(cnv_by_band_clone2[rows_clone34_gain, ])
rows_clone34_loss <- which(cnv_by_band_clone2$F < -0.2)
clone34_loss <- row.names(cnv_by_band_clone2[rows_clone34_loss, ])

clone_gain <- unique(c(clone1_gain,clone2_gain,clone31_gain,clone32_gain,clone33_gain,clone34_gain))
clone_loss <- unique(c(clone1_loss,clone2_loss,clone31_loss,clone32_loss,clone33_loss,clone34_loss))
clone_gain <- sort(clone_gain)
clone_loss <- sort(clone_loss)
clone_var <- c(clone_gain,clone_loss)

cnv_clone.df <- cnv_by_band_clone2[clone_var,]
cnv_clone.df2 <- apply(cnv_clone.df, c(1,2), function(x) ifelse(x > 0.4, 2, 
                                                                ifelse(x > 0.2, 1, 
                                                                       ifelse(x < -0.4, -2,
                                                                              ifelse(x < -0.2, -1, 0)))))
cnv_clone.df2 <- as.data.frame(cnv_clone.df2)
write.csv(cnv_clone.df2, "summary_5states_final.csv", row.names = T)
## this file will be used to construct tumor clone tree

#################################################################################
### DEG analysis
breast_ST_obj <- Load10X_Spatial(
  data.dir="./breastA",
  filename = "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")

breast_ST_obj <- SCTransform(breast_ST_obj, assay = "Spatial", verbose = FALSE)
breast_ST_obj <- RunPCA(breast_ST_obj, assay = "SCT", verbose = FALSE)
#breast_ST_obj <- FindNeighbors(breast_ST_obj, reduction = "pca", dims = 1:30)
#breast_ST_obj <- FindClusters(breast_ST_obj, verbose = FALSE)
breast_ST_obj <- RunUMAP(breast_ST_obj, reduction = "pca", dims = 1:30)

# set my own cluster
id <- breast_ST_obj@assays$Spatial@counts@Dimnames[[2]]
breast_label <- breast_label[match(id, breast_label$barcode),]
breast_cluster <- breast_label$clone3
breast_cluster2 <- breast_label$clone6
breast_ST_obj$clone3 <- breast_cluster
breast_ST_obj$clone6 <- breast_cluster2
Idents(breast_ST_obj) <- breast_ST_obj$clone3

# visualize
SpatialDimPlot(breast_ST_obj, group.by='clone3',label = FALSE, label.size = 3, cols = c('Normal'='#F2DCC9',
                                                                                      'Tumor clone1'='#B10026',
                                                                                      'Tumor clone2'='#E31A1C',
                                                                                      'Tumor clone3'='#FC4E2A'))

SpatialDimPlot(breast_ST_obj, group.by='clone6',label = FALSE, label.size = 3, cols = c('Normal'='darkgray',
                                                                                        'A'='#B10026',
                                                                                        'B'='#E31A1C',
                                                                                        'C'='#FC4E2A',
                                                                                        'D'='#FD8D3C',
                                                                                        'E'='#FED976',
                                                                                        'F'='#FFFFCC'))

## DEG analysis: all clusters
markers <- FindAllMarkers(breast_ST_obj, only.pos = TRUE)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> topgenes

DoHeatmap(breast_ST_obj, features = topgenes$gene)+NoLegend()

## DEGs between pairwise clones
markers_c12 <- FindMarkers(breast_ST_obj, ident.1 = "Tumor clone1", ident.2 = "Tumor clone2",
                                  logfc.threshold = 0.05,
                                  min.pct = 0.1)

markers_c13 <- FindMarkers(breast_ST_obj, ident.1 = "Tumor clone1", ident.2 = "Tumor clone3",
                                  logfc.threshold = 0.05,
                                  min.pct = 0.1)

markers_c23 <- FindMarkers(breast_ST_obj, ident.1 = "Tumor clone2", ident.2 = "Tumor clone3",
                                  logfc.threshold = 0.05,
                                  min.pct = 0.1)

write.csv(markers_c12.markers,"markers_c1_c2.csv",row.names = T)
write.csv(markers_c13.markers,"markers_c1_c3.csv",row.names = T)
write.csv(markers_c23,"markers_c2_c3.csv",row.names = T)

#################################################################################
### GSEA
# Read the Hallmark gene set GMT file
gmt_file <- "h.all.v2023.2.Hs.symbols.gmt"
hallmark_genesets <- read.gmt(gmt_file)

# reading in DEG data
df = read.csv("markers_c1_c2.csv", header=TRUE)

# we want the log2 fold change 
logFC <- df$avg_log2FC

# name the vector
symbol <- df$X
data <- data.frame(SYMBOL=symbol, logFC=logFC)
geneList <- data$logFC
names(geneList) <- data$SYMBOL

# Make sure geneList is sorted from high to low according to expression change value
geneList <- sort(geneList, decreasing = TRUE)

# GSEA analysis
## clone1 vs clone2
gsea_result_12 <- GSEA(geneList, 
                       nPerm = 10000,
                       TERM2GENE = hallmark_genesets,
                       minGSSize = 10,  
                       maxGSSize = 500,  
                       pAdjustMethod = "BH",  
                       pvalueCutoff = 0.05,  
                       verbose = FALSE,
                       seed=1234)  

# see GSEA result
head(gsea_result_12)
dotplot(gsea_result_12, showCategory=10, split=".sign") + facet_wrap(.~.sign, scales = "free")+ggtitle("Tumor subclone 1 vs. subclone 2")

## clone1 vs clone3
df = read.csv("markers_c1_c3.csv", header=TRUE)
logFC <- df$avg_log2FC
symbol <- df$X
data <- data.frame(SYMBOL=symbol, logFC=logFC)
geneList <- data$logFC
names(geneList) <- data$SYMBOL
geneList <- sort(geneList, decreasing = TRUE)
gsea_result_13 <- GSEA(geneList, 
                       nPerm = 10000,
                       TERM2GENE = hallmark_genesets,
                       minGSSize = 10,  
                       maxGSSize = 500,  
                       pAdjustMethod = "BH",  
                       pvalueCutoff = 0.05,  
                       verbose = FALSE,
                       seed=1234)  
head(gsea_result_13)
dotplot(gsea_result_13, showCategory=10, split=".sign") + facet_wrap(.~.sign, scales = "free")+ggtitle("Tumor subclone 1 vs. subclone 3")

# clone2 vs clone3
df = read.csv("markers_c2_c3.csv", header=TRUE)
logFC <- df$avg_log2FC
symbol <- df$X
data <- data.frame(SYMBOL=symbol, logFC=logFC)
geneList <- data$logFC
names(geneList) <- data$SYMBOL
geneList <- sort(geneList, decreasing = TRUE)

gsea_result_23 <- GSEA(geneList, 
                       nPerm = 10000,
                       TERM2GENE = hallmark_genesets,
                       minGSSize = 10,  
                       maxGSSize = 500,  
                       pAdjustMethod = "BH",  
                       pvalueCutoff = 0.05,  
                       verbose = FALSE,
                       seed=1234)  
head(gsea_result_23)
dotplot(gsea_result_23, showCategory=10, split=".sign") + facet_wrap(.~.sign, scales = "free")+ggtitle("Tumor subclone 1 vs. subclone 3")

########################
### plot
df_c1_c2 <- as.data.frame(gsea_result_12)
df_c1_c3 <- as.data.frame(gsea_result_13)
df_c2_c3 <- as.data.frame(gsea_result_23)

# To distinguish in the dot plot, add a comparison group column to each data frame
df_c1_c2$Comparison <- 'c1 vs c2'
df_c1_c3$Comparison <- 'c1 vs c3'
df_c2_c3$Comparison <- 'c2 vs c3'

# combine dataframes
combined_df <- rbind(df_c1_c2, df_c1_c3, df_c2_c3)

ggplot(combined_df, aes(x = Comparison, y = ID, size = -log10(p.adjust), color = NES)) +
  geom_point() +
  theme_classic() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "GSEA Results", x = "Comparison", y = "Gene Set", size = "-log10(p.adjust)", color = "NES")

################################################################################
######### tumor clonal tree reconstruction ##########

cnv_summary <- read.csv("summary_5states_final.csv")
phy_dat <- phyDat(t(as.matrix(cnv_summary)), type="USER", levels=c(-2,-1,0,1,2))
tree <- pratchet(phy_dat, trace=0) |> acctran(phy_dat)
parsimony(tree, phy_dat)
tre.pars <- optim.parsimony(tree, phy_dat)

anc.acctran <- ancestral.pars(tree, phy_dat, "ACCTRAN")
anc.mpr <- ancestral.pars(tree, phy_dat, "MPR")

plotAnc(tree, anc.mpr, col="white", pos = NULL)
title("Phylogenetic tree based on 5 CNV status")

################################################################################
######### trajectory analysis ##########
### create object
breast_cnv_mat <- breast_cnv_obj@expr.data
breast_seurat_obj <- CreateSeuratObject(counts = breast_cnv_mat, project = "human breast cancer", min.cells = 0)
# skip scale step
#breast_seurat_obj <- SetAssayData(breast_seurat_obj, assay = 'RNA', layer = "scale.data", new.data = breast_cnv_mat)
breast_seurat_obj <- SetAssayData(breast_seurat_obj, assay = 'RNA', slot = "scale.data", new.data = breast_cnv_mat)

all_genes <- rownames(breast_seurat_obj)
breast_seurat_obj <- RunPCA(breast_seurat_obj, features = all_genes)

#GRAPH based clustering (adjust resolution parameter)
breast_seurat_obj <- FindNeighbors(breast_seurat_obj, dims = 1:10)
breast_seurat_obj <- FindClusters(breast_seurat_obj, resolution = 0.3)

# run UMAP
breast_seurat_obj <- RunUMAP(breast_seurat_obj, dims = 1:15)

### my own cluster
## 3 clones
id <- colnames(breast_cnv_mat)
breast_label <- breast_label[match(id, breast_label$barcode),]
breast_cluster <- breast_label$clone3

breast_seurat_obj$my_cluster <- breast_cluster
DimPlot(breast_seurat_obj, reduction = "umap", group.by='my_cluster', label = 'TRUE')

## 6 clones
breast_cluster2 <- breast_label$clone6
breast_seurat_obj$my_cluster2 <- breast_cluster2
DimPlot(breast_seurat_obj, reduction = "umap", group.by='my_cluster2', label = 'TRUE')

sce <- as.SingleCellExperiment(breast_seurat_obj)
sce

sce_slingshot1 <- slingshot(sce,      
                            reducedDim = 'UMAP', 
                            clusterLabels = sce$my_cluster,  
                            #start.clus = 'Normal',      
                            approx_points = 150)

SlingshotDataSet(sce_slingshot1) 

## Extract the matrix of pseudotime values or cells' weights along each lineage.
pse_time <- as.data.frame(slingPseudotime(sce_slingshot1))
breast_label <- breast_label[match(row.names(pse_time),breast_label$barcode),]
breast_label$time <- pse_time$Lineage1

library(ggplot2)
ggplot(breast_label, aes(x = width, y = height, color =time)) +
  geom_point(size=2) +
  scale_y_reverse()+
  scale_color_gradient(low = "#F1ECDC", high = "#C0765A") + 
  labs(title = "Cell Pseudotime", color = "Pseudotime") +
  theme_void() +
  coord_equal()+
  theme(plot.title = element_text(hjust = 0.5, size =10),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6))




