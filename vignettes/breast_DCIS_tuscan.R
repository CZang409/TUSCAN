library(TUSCAN)
library(Seurat)
library(magick)
library(dplyr)

# load data
count.data <- Seurat::Read10X_h5("Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5")
loc.data <- read.csv("tissue_positions_list.csv")
loc.data <- loc.data %>% filter (X0==1)
cell.id <- loc.data[,1]
loc.data <- data.frame(x=loc.data[,6],y=loc.data[,5])
row.names(loc.data) <- cell.id
image_path <- "Visium_FFPE_Human_Breast_Cancer_image.tif"
img <- image_read(image_path)
loc.data <- loc.data[match(colnames(count.data),row.names(loc.data)),]

# create TUSCAN object
obj <- createTUSCANObject(counts = count.data, 
                          location = loc.data, 
                          img = img, 
                          project="TUSCAN", 
                          r=24)

# dimension reduction
obj <- reduceDim(object = obj, nhvg=2000,d=20,epochs=50,batch_size=64, seed = 1234)

# pre-clustering
obj <- PreCluster(object = obj, method = 'louvain', res=1.5, seed = 1234)
plotSpatial(obj = obj, class = "cluster", size = 2)

# find normal reference
FindNormalCluster(obj,w=0.3)

# infer copy number
cnv_obj <- runCNV(object = obj, 
                  ref_group_names = c("5"),
                  min_ave_counts_per_gene=0.01,
                  min_spots_per_gene=3,
                  chr_exclude=c('chrX', 'chrY', 'chrM'))

plotCNVheatmap(object = cnv_obj, title = "./TUSCAN/breast_heatmap", output_format = "png")

# classify tumor
tumor <- TumorCluster(object = cnv_obj, seed = 1234)
plotSpatial(obj = tumor, class = "tumor", size = 2)

