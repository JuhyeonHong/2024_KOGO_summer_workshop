# Download R packages to be used in analysis #
options(timeout = 600000000) ### set this to avoid timeout error
install.packages("hdf5r") 
install.packages("dplyr")
install.packages('Seurat')
install.packages("ggplot2")
install.packages("remotes")
install.packages("anndata")
install.packages("cowplot")
install.packages("devtools")
devtools::install_github("thomasp85/patchwork")
devtools::install_github('satijalab/seurat-data')
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
remotes::install_github("sqjin/CellChat")
remotes::install_github("jbergenstrahle/STUtility")

# Load packages #
library(hdf5r)
library(dplyr)
library(Seurat)
library(ggplot2)
library(remotes)
library(anndata)
library(cowplot)
library(devtools)
library(patchwork)
library(SeuratData)
library(spacexr)
library(SeuratData)
library(CellChat)
library(STutility)

# Chapter 4. Integrative Analysis of Spatial Datasets: Seurat --------


## Load and processing Visium --------

### Load Brain Data --------
setwd("2023_KOGO_summer_workshop_data")
options(timeout = 600000000)
InstallData("stxBrain")
brain = LoadData("stxBrain", type = "anterior1")

### Explore Seurat object --------
brain

### Explore metadata --------
brain@meta.data %>% head(3)

### Explore coordinates --------
brain@images$anterior1@coordinates %>% head(3)
GetTissueCoordinates(brain) %>% head(3)
SpatialPlot(brain) + ggplot2::theme_minimal()

### nCount violin and spatial plot  --------
plot1 = VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 = SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

### Pre-processing: SCT normalize --------
brain = SCTransform(brain, assay = "Spatial", verbose = FALSE)

### Visualization of gene expression --------
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

### Dimension reduction and clustering --------
brain = RunPCA(brain, assay = "SCT", verbose = FALSE)
brain = FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain = FindClusters(brain, verbose = FALSE)
brain = RunUMAP(brain, reduction = "pca", dims = 1:30)

### Visualization of each cluster --------
p1 = DimPlot(brain, reduction = "umap", label = TRUE)
p2 = SpatialDimPlot(brain, label = TRUE, label.size = 3)
wrap_plots(p1, p2)

### Spatial Dimplot of some cluster --------
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 5, 8)), 
               facet.highlight = TRUE, ncol = 3)

### Identification of spatially variable features --------
de_markers = FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)
saveRDS(brain, "./KOGO_visium_brain_processed.rds")

### Subset out anatomical regions --------
cortex = subset(brain, idents = c(1, 2, 3, 4, 6, 7))

### now remove additional cells, use SpatialDimPlots to visualize what to remove
cortex = subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex = subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex = subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 = SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 = SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2

### After subsetting, we renormalize cortex --------
cortex = SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
saveRDS(cortex, "./KOGO_visium_cortex_processed.rds")


## Load and proprocessing single cell --------

### Load reference single cell data --------

allen_reference= readRDS("./KOGO_sc_allen_cortex.rds")
allen_reference = SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% 
  RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
allen_reference

### the annotation is stored in the 'subclass' column of object metadata --------
DimPlot(allen_reference, group.by = "subclass", label = TRUE)


## Find TransferAnchors --------
### Find TransferAnchors --------
anchors = FindTransferAnchors(reference = allen_reference, query = cortex, 
                              normalization.method = "SCT")
anchors@anchors %>% head(15)
anchors

### Acquire prediction assay --------
predictions.assay = TransferData(anchorset = anchors, 
                                 refdata = allen_reference$subclass, 
                                 prediction.assay = TRUE, 
                                 weight.reduction = cortex[["pca"]], dims = 1:30)
predictions.assay[,1:3] %>% t()

### Insert prediction assay into the cortex object --------
cortex[["predictions"]] = predictions.assay
DefaultAssay(cortex) = "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE)

### Spatial Feature Plot of predicted cell type proportion --------
pdf('./integreated_spatialfeautrePlot.pdf', width = 10, height = 5)
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT"), 
                   pt.size.factor = 1, ncol = 5, crop = FALSE, alpha = c(0.1, 1))
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT"), 
                   pt.size.factor = 1, ncol = 5, crop = TRUE, alpha = c(0.1, 1))
dev.off()


## Transfer annotation --------
### Transfer annotation --------
predictions = TransferData(anchorset = anchors, refdata = allen_reference$subclass,
                           weight.reduction = cortex[["pca"]], dims = 1:30)
cortex = AddMetaData(cortex, metadata = predictions)
cortex$predicted.id = factor(cortex$predicted.id)
cortex = SetIdent(cortex, value = "predicted.id")
SpatialDimPlot(cortex, label = T, label.size = 3)

### Remove all objects --------
rm(allen_reference)
rm(brain)
rm(cortex)
rm(anchors)
rm(predictions);
rm(predictions.assay);
gc()



# Chapter 5. Deconvolution Analysis - SpaceXR --------


## Process single cell for RCTD input --------

### Load RCTD input reference dataset (single cell) --------
allen_reference = readRDS("./KOGO_sc_allen_cortex.rds")
allen_reference = subset(allen_reference, subclass != c("CR"))

### Prepare single cell dataset for RCTD input --------
counts_sc = allen_reference$RNA@counts
meta_sc = allen_reference@meta.data
allen_reference = SetIdent(allen_reference, value = "subclass")

allen_reference = RenameIdents(allen_reference, "L2/3 IT" = "L2_3 IT") #Rename the cluster due to prohibited characters “/”
allen_reference$subclass = allen_reference@active.ident

annotation_sc  = allen_reference@meta.data$subclass
names(annotation_sc) = rownames(meta_sc)
annotation_sc = as.factor(annotation_sc)
nUMI_sc = meta_sc$nCount_RNA
names(nUMI_sc) = rownames(meta_sc)
reference = Reference(counts_sc, annotation_sc, nUMI_sc)


## Process spatial for RCTD input --------

### Process spatial dataset for RCTD input --------
brain = readRDS("./KOGO_visium_brain_processed.rds")
counts_visium = brain@assays$Spatial@counts
nUMI_visium = colSums(counts_visium)
coords_visium = brain@images$anterior1@coordinates[,c("col","row")]
query = SpatialRNA(coords_visium, counts_visium, nUMI_visium)


## Run RCTD --------

### Run RCTD in doublet mode --------
#RCTD = create.RCTD(query, reference, max_cores = 8)
#RCTD = run.RCTD(RCTD, doublet_mode = 'doublet')
#brain = AddMetaData(brain, metadata = RCTD@results$results_df)
#write.csv(RCTD@results$results_df,"./KOGO_visium_brain_RCTD.csv")


## Process output for deconvolution --------

### Process RCTD decomposed file --------
RCTD = read.table("./KOGO_visium_brain_RCTD.csv", sep=",", header=TRUE)
meta_RCTD = RCTD[,c("X", "first_type")]
colnames(meta_RCTD) = c("barcodes", "spacexr_first_type")
rownames(meta_RCTD) = meta_RCTD$barcodes; meta_RCTD$barcodes = NULL
brain = AddMetaData(brain, meta_RCTD)
brain = SetIdent(brain, value="spacexr_first_type")

### Brain spatial map of predicted cell type by RCTD --------
SpatialDimPlot(brain, label = T, label.size = 4, repel = T)

### Cortex spatial map of predicted cell type by RCTD --------
cortex = subset(brain, seurat_clusters %in% c(1, 2, 3, 4, 6, 7))
cortex = subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex = subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex = subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
SpatialDimPlot(cortex, label = T, label.size = 4, repel = T)

saveRDS(cortex, "./KOGO_visium_cortex_annotated.rds") # This file will be used in the Squidpy and CellChat Analysis

### Remove all objects before starting next chapter --------
gc()



# Chapter 6. Neighborhood Analysis of Co-occurrence - Squidpy --------


## Convert / load data --------

### Convert data from Seurat object to Anndata --------
cortex = readRDS("./KOGO_visium_cortex_annotated.rds")
cortex = AddMetaData(cortex, cortex@images$anterior1@coordinates)
cortex = UpdateSeuratObject(cortex)
SaveH5Seurat(cortex, filename = "./KOGO_visium_cortex_anndata.h5Seurat")
Convert("./KOGO_visium_cortex_anndata.h5Seurat", dest = "h5ad")



# Chapter 7. Cell-cell interaction: Cellchat  --------


## Load dataset  --------

### Load data --------
visium.brain = readRDS('./KOGO_visium_cortex_annotated.rds')
visium.brain$spacexr_first_type  = factor(visium.brain$spacexr_first_type, 
                                          levels = c("Astro", "L2_3 IT", "L4", "L5 IT", "L5 PT", 
                                                   "L6 CT", "L6 IT", "L6b", "Lamp5", "Meis2", 
                                                   "Oligo", "Peri", "Sncg", "Sst", "VLMC"))
Idents(visium.brain) = visium.brain$spacexr_first_type
colors = scPalette(nlevels(visium.brain))
names(colors) = c("Astro", "L2_3 IT", "L4", "L5 IT", "L5 PT", 
                  "L6 CT", "L6 IT", "L6b", "Lamp5", "Meis2", 
                  "Oligo", "Peri", "Sncg", "Sst", "VLMC")
SpatialDimPlot(visium.brain, label = T, label.size = 3, cols = colors)

### Prepare input data for CellChat analysis --------
data.input = GetAssayData(visium.brain, slot = "data", assay = "SCT")
meta = data.frame(labels = Idents(visium.brain), 
                  row.names = names(Idents(visium.brain)))

unique(meta$labels) # check the cell labels

### Load spatial imaging information --------
spatial.locs = GetTissueCoordinates(visium.brain, scale = NULL, 
                                    cols = c("imagerow", "imagecol")) 
scale.factors = jsonlite::fromJSON(txt = "./KOGO_scalefactors_json.json")
scale.factors = list(spot.diameter=65, 
                     spot = scale.factors$spot_diameter_fullres, 
                     fiducial = scale.factors$fiducial_diameter_fullres, 
                     hires = scale.factors$tissue_hires_scalef, 
                     lowres = scale.factors$tissue_lowres_scalef)

### Create a CellChat object --------
cellchat = createCellChat(object = data.input, 
                          meta = meta, 
                          group.by = "labels", 
                          datatype = "spatial", 
                          coordinates = spatial.locs, 
                          scale.factors = scale.factors)
cellchat

### Set the ligand-receptor interaction database --------
CellChatDB = CellChatDB.mouse
cellchat@DB = CellChatDB


## Preprocessing --------

### Preprocessing the expression data for cell-cell communication analysis --------
cellchat = subsetData(cellchat)
future::plan("multisession", workers = 4) # do parallel
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)


## Infer cell-cell communication network --------

### Compute the communication probability and infer cellular communication network --------
cellchat = computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                             distance.use = TRUE, interaction.length = 200, 
                             scale.distance = 0.01)
cellchat = filterCommunication(cellchat, min.cells = 10)

saveRDS(cellchat, file = './KOGO_visium_cortex_prob_cellchat.rds')
#cellchat <- readRDS("./KOGO_visium_cortex_prob_cellchat.rds")

### Infer the cell-cell communication at a signaling pathway level --------
cellchat = computeCommunProbPathway(cellchat)

### Calculate the aggregated cell-cell communication network --------
cellchat = aggregateNet(cellchat)

### Visualization of the aggregated cell-cell communication network --------
groupSize = as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = rowSums(cellchat@net$count), 
                 weight.scale = T, label.edge = F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = rowSums(cellchat@net$weight), 
                 weight.scale = T, label.edge = F, 
                 title.name = "Interaction weights/strength")


## Visualization --------
### Identify ligand-receptor pairs between L6b and Oligo --------
CellChat::netVisual_bubble(cellchat, 
                           sources.use =  c(8), 
                           targets.use = c(11), 
                           remove.isolate = FALSE, angle.x = 90, 
                           thresh = 0.05,
                          ) + coord_flip()

### Visualize the inferred communication network of signaling pathways - Circle plot --------
pathways.show = c("PSAP")

par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    layout = "circle")

### Visualize the inferred communication network of signaling pathways - Spatial plot --------
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, 
                    vertex.label.cex = 3.5)

### Compute the network centrality scores --------
cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")

par(mfrow = c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)

### Visualize the network centrality scores - Spatial plot --------
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    layout = "spatial", edge.width.max = 2, 
                    alpha.image = 0.2, vertex.weight = "outgoing", 
                    vertex.size.max = 3, vertex.label.cex = 3.5)



# Chapter 8. Visualization of blended spatial plot of several features: STUtility --------


## Construct STUtility object --------

### Construct STUtility object --------
setwd("./KOGO_visium_STUtility_input_files/")

infoTable = data.frame(samples = c("./V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5"),
                       spotfiles = c("./spatial/tissue_positions_list.csv"),
                       imgs = c("./spatial/tissue_hires_image.png"),
                       json = c("./spatial/scalefactors_json.json"), 
                       stringsAsFactors = FALSE)

Mouse_Brain_STUtility = InputFromTable(infotable = infoTable, platform = "Visium")

### Function used to read HE images in jpeg or png format --------
Mouse_Brain_STUtility = LoadImages(Mouse_Brain_STUtility, time.resolve = FALSE, verbose = TRUE)


## Blended Spatial Feature Plot --------

### Feature plot of Ligand and Receptor and Blended Spatial Feature Plot --------
plot_grid(
    ST.FeaturePlot(Mouse_Brain_STUtility, features = c("Psap"), cols = c("white","green")),
    ST.FeaturePlot(Mouse_Brain_STUtility, features = c("Gpr37l1"), cols = c("white","red")),
    ST.FeaturePlot(Mouse_Brain_STUtility, features = c("Psap","Gpr37l1"), 
                   blend = TRUE, channels.use = c("green","red")), ncol = 3)
saveRDS(Mouse_Brain_STUtility, "../KOGO_visium_cortex_STUtility.rds")



