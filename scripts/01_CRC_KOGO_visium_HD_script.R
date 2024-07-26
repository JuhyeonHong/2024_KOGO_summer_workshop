# Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)
library(arrow)

# Set path
setwd("~/2024_KOGO_visium/")


## Visium V2
# Load the data
crc.vs <- Load10X_Spatial(data.dir = "data/visium_v2/")
DefaultAssay(crc.vs) <- "Spatial"

# Explore the data
crc.vs
vln.plot <- VlnPlot(crc.vs, features = "nCount_Spatial", pt.size = 0, raster = FALSE) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(crc.vs, features = "nCount_Spatial", pt.size.factor = 5) + theme(legend.position = "right") 
vln.plot | count.plot

# Preprocessing
crc.vs <- SCTransform(crc.vs, assay = "Spatial", verbose = FALSE)
crc.vs <- RunPCA(crc.vs, assay = "SCT", verbose = FALSE)
crc.vs <- FindNeighbors(crc.vs, reduction = "pca", dims = 1:30)
crc.vs <- FindClusters(crc.vs, verbose = FALSE, resolution = 0.5)
crc.vs <- RunUMAP(crc.vs, reduction = "pca", dims = 1:30)
crc.vs

# Visualize
library(CellChat)
nlevels(crc.vs) # 10
colors = scPalette(nlevels(crc.vs))
names(colors) = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
p1 <- DimPlot(crc.vs, reduction="umap",cols=colors, label=F) + theme(aspect.ratio = 1)
p2 <- SpatialDimPlot(crc.vs, label=F,repel=T,cols=colors, pt.size.factor = 5) + theme(aspect.ratio=1)
p1|p2
SpatialFeaturePlot(crc.vs, features = c("CEACAM6"), pt.size.factor = 5)


## Visium HD
# Subset the visium HD data
# DO NOT RUN #
# crc.vshd <- Load10X_Spatial(data.dir = "data/visium_hd/binned_outputs/square_008um/") # takes about 1min
# DefaultAssay(crc.vshd) <- "Spatial"
# coords <- GetTissueCoordinates(crc.vshd)
# head(coords) %>% View()
# max_x <- max(coords$x)
# max_y <- max(coords$y)
# min_x <- min(coords$x)
# min_y <- min(coords$y)
# x_min_threshold <- min_x + (max_x - min_x) * 0.1
# y_min_threshold <- min_y + (max_y - min_y) * 0.7
# x_max_threshold <- min_x + (max_x - min_x) * 0.3
# y_max_threshold <- min_y + (max_y - min_y) * 0.9
# selected_spots <- coords[coords$x > x_min_threshold & coords$x < x_max_threshold & coords$y > y_min_threshold & coords$y < y_max_threshold, ]
# print(selected_spots)
# selected_barcodes <- rownames(selected_spots)
# crc.vshd_subset <- subset(crc.vshd, cells = selected_barcodes)

## Preprecessing (takes large memory to analyze)
# DO NOT RUN # 
# crc.vshd_subset <- subset(crc.vshd, nCount_Spatial > 0)
# crc.vshd_subset <- SCTransform(crc.vshd_subset, assay = "Spatial", verbose = FALSE)
# crc.vshd_subset <- RunPCA(crc.vshd_subset, assay = "SCT", verbose = FALSE)
# crc.vshd_subset <- FindNeighbors(crc.vshd_subset, reduction = "pca", dims = 1:30)
# crc.vshd_subset <- FindClusters(crc.vshd_subset, verbose = FALSE, resolution = 0.3)
# crc.vshd_subset <- RunUMAP(crc.vshd_subset, reduction = "pca", dims = 1:30)
# saveRDS(crc.vshd_subset, 'results/visium_hd/visium_hd_subset_preprocessed.rds')


# Bring the subset Visium HD data and explore
crc.vshd <- readRDS('results/visium_hd/visium_hd_subset_preprocessed.rds')
crc.vshd
SpatialFeaturePlot(crc.vshd, features = "nCount_Spatial", pt.size.factor = 15) + theme(legend.position = 'right')

# Visualize
SpatialFeaturePlot(crc.vshd, features = c("CEACAM6", "PIGR", "COL1A1"), pt.size.factor = 15)

nlevels(crc.vshd) # 14
colors = scPalette(nlevels(crc.vshd))
names(colors) = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13')
p1 <- DimPlot(crc.vshd, reduction="umap",cols=colors, label=F) + theme(aspect.ratio = 1)
p2 <- SpatialDimPlot(crc.vshd, label=F,repel=T,cols=colors, pt.size.factor = 15) + theme(aspect.ratio=1)
p1|p2

###------------------------------------------------------------------------------------------
## Deconvolution with RCTD using scRNA-seq reference data
## Load single cell RNA-seq data 
# DO NOT RUN #
# crcsc_count <- Read10X(data.dir = "data/flexsc/filtered_feature_bc_matrix/")
# MetaData <- read.csv("data/metadata/SingleCell_MetaData.csv")
# MetaData_new <- subset(MetaData, Level1 != "QC_Filtered")
# table(MetaData_new$Level1)
# Spacexr restriction
# Spacexr restriction, Clusters with > 25 cells
# KpIdents<-names(which(table(MetaData_new$Level1)>25))
# MetaData_new<-MetaData_new[MetaData_new$Level1%in%KpIdents,]
# crcsc_count<-crcsc_count[,MetaData_new$Barcode]
# 
# Fix cell type labels as spacexr doesn't allow special characters (i.e. spaces)
# CTRef<-MetaData_new$Level1
# CTRef<-gsub("/","_",CTRef)
# CTRef<-as.factor(CTRef)
# names(CTRef)<-MetaData_new$Barcode
# Build reference object
# reference <- Reference(crcsc_count[,names(CTRef)], CTRef , colSums(crcsc_count)) # spacexr function
#saveRDS(reference, 'results/visium_hd/20240724_rctd_sc_reference.rds')

reference <- readRDS('results/visium_hd/visium_hd_subset_rctd_sc_reference.rds')

## Load Visium HD count and coordinate data
counts <- crc.vshd@assays$Spatial$counts
coords <- GetTissueCoordinates(crc.vshd)
coords<-coords[,c("x", "y")]
head(rownames(coords))
nUMI <- colSums(counts)
puck <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(puck@counts)

# DO NOT RUN #
# RCTD <- create.RCTD(puck, reference, max_cores = 24)
# RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet')
# saveRDS(RCTD,'results/visium_hd/visium_hd_subset_rctd_result.rds')
RCTD <- readRDS('results/visium_hd/visium_hd_subset_rctd_result.rds')
RCTD@results$results_df %>% View()

# Clear unnecessary variables
rm(list=setdiff(ls(), c("crc.vshd", "RCTD")))

# Add deconvolution metadata from RCTD results
crc.vshd <- AddMetaData(crc.vshd, metadata = RCTD@results$results_df)
crc.vshd$first_type <- as.character(crc.vshd$first_type)
sum(is.na(crc.vshd$first_type))
table(crc.vshd$first_type)
crc.vshd$first_type[is.na(crc.vshd$first_type)] <- "Undetermined"
crc.vshd$first_type <- as.factor(crc.vshd$first_type)
crc.vshd@active.ident <- crc.vshd$first_type
head(crc.vshd@meta.data)

# Visualize the annotated visium HD data
nlevels(crc.vshd)
colors = scPalette(nlevels(crc.vshd))
names(colors) = c("B cells", "Endothelial", "Fibroblast", "Intestinal Epithelial", "Myeloid", "Neuronal", "Smooth Muscle", "T cells", "Tumor", "Undetermined")

SpatialDimPlot(crc.vshd, label=F, repel=T, cols=colors, pt.size.factor = 15)


# Identify diverse cell types
SpatialDimPlot(crc.vshd,
               cells.highlight = CellsByIdentities(object = crc.vshd,
                                                   idents = c("Tumor", "Fibroblast", "Myeloid")),
               facet.highlight = T, ncol=3, pt.size.factor = 15)

# Marker genes tumor, goblet cells and enterocytes, fibroblasts, macrophage
SpatialFeaturePlot(crc.vshd, features = c("CEACAM6", "PIGR", "COL1A1", "C1QC"), ncol = 2, pt.size.factor = 15)


