# Download R packages to be used in analysis #
options(timeout = 600000000) ### set this to avoid timeout error
remotes::install_github("sqjin/CellChat")

# Load packages #
library(CellChat)
library(Seurat)
library(dplyr)
library(spacexr)
library(ggpubr)
library(SoupX)

## Single cell pre-processing and annotation process are provided at singlecell_preprocess_script.R
## Visium pre-processing is provided at visium_preprocess_script.R

# Chapter 2. Deconvolution Analysis - SpaceXR --------
## Process single cell for RCTD input --------

### Load RCTD input reference dataset (single cell) --------

## Process single cell for RCTD input --------
### Process single cell dataset for RCTD input --------
setwd("")
Breast_sc <- readRDS("./R_object/Bioinfo2023_Breast_singlecell.rds")

### Prepare single cell dataset for RCTD input --------
counts_sc = Breast_sc$RNA@counts

annotation_sc  = Breast_sc$annotation
names(annotation_sc) = rownames(Breast_sc@meta.data)
annotation_sc = as.factor(annotation_sc)

nUMI_sc = Breast_sc@meta.data$nCount_RNA
names(nUMI_sc) = rownames(Breast_sc@meta.data)

reference = Reference(counts_sc, annotation_sc, nUMI_sc)
gc()

## Process spatial for RCTD input --------

### Process spatial dataset for RCTD input --------
breast_visium = readRDS("./R_object/Bioinfo2023_Breast_visium.rds")

coords_visium = breast_visium@images$slice1@coordinates[,c("col","row")]

counts_visium = breast_visium@assays$Spatial@counts

nUMI_visium = colSums(counts_visium)

query = SpatialRNA(coords_visium, counts_visium, nUMI_visium)
gc()


## Run RCTD --------

### Run RCTD in doublet mode --------
RCTD = create.RCTD(query, reference, max_cores = 8)
RCTD = run.RCTD(RCTD, doublet_mode = 'doublet')
# RCTD = readRDS("./R_object/Bioinfo2023_Breast_RCTD.rds")
RCTD_results = RCTD@results$results_df
breast_visium = AddMetaData(breast_visium, metadata = RCTD_results)

## Visualize decomposed cell types --------

### Process RCTD decomposed file --------
breast_visium = SetIdent(breast_visium, value="first_type")
table(breast_visium$first_type)

### breast_visium spatial map of predicted cell type by RCTD --------
### Remove spots with lower than 100 nUMI
breast_visium <- subset(breast_visium, first_type %in% names(table(breast_visium$first_type)))
breast_visium$first_type <- factor(breast_visium$first_type)
SpatialDimPlot(breast_visium)
SpatialDimPlot(breast_visium, cells.highlight = CellsByIdentities(object = breast_visium,
                                                                  idents = c('DCIS #1','DCIS #2','Invasive','Mixed')), facet.highlight = TRUE)
SpatialDimPlot(breast_visium, cells.highlight = CellsByIdentities(object = breast_visium,
                                                                  idents = c('CD4+ T cell', 'CD8+ T cell', 'B cell', 'Myeloid cell',  'Plasma cell', 'Stromal cell')), facet.highlight = TRUE, ncol = 3)

#saveRDS(breast_visium, "./R_object/Bioinfo2023_Breast_visium_final.rds")

### Remove all objects before starting next chapter --------
gc()



## Microenvironment construction of tumor boundary is provided at visium_cottrazm_script.R

# Chapter 3. Cell-cell interaction: Cellchat --------


## Load data --------

### Load cell type annotated visium data and visiualization --------
visium.breast = readRDS("./R_object/Bioinfo2023_Breast_visium_final.rds")
visium.breast$first_type = factor(visium.breast$first_type, 
                                  levels = c("B cell", "CD4+ T cell", "CD8+ T cell", "DCIS #1", "DCIS #2", "Mixed", "Myeloid cell", "Plasma cell", "Stromal cell", "Invasive"))
Idents(visium.breast) = visium.breast$first_type
colors = scPalette(nlevels(visium.breast))
names(colors) = c("B cell", "CD4+ T cell", "CD8+ T cell", "DCIS #1", "DCIS #2", "Mixed", "Myeloid cell", "Plasma cell", "Stromal cell", "Invasive")
SpatialDimPlot(visium.breast, label = F, cols = colors)

### Load tumor annotated visium data --------
visium.tumor = readRDS('./R_object/Bioinfo2023_Breast_visium_TumorST.rds')
visium.tumor = subset(visium.tumor, nCount_Spatial > 100)
table(visium.tumor@meta.data$tumor_annotation)

### Subset only tumor boundaries --------
visium.breast@meta.data$tumor_annotation = visium.tumor@meta.data$tumor_annotation
Idents(visium.breast)= visium.breast@meta.data$tumor_annotation
visium.boundary = subset(visium.breast, idents = "Bdy")

Idents(visium.boundary)= visium.boundary@meta.data$tumor_annotation
names(colors) = "Bdy" 
SpatialDimPlot(visium.boundary, label = F, cols = colors)

### Load layer annotated LoupeBrowser data --------
loupe.data = read.csv("./Raw_file/LoupeBrowser/layer_annotation.csv", header = T, row.names = 1)
table(loupe.data$layer_annotation)

### Subset only tumor boundaries and layers --------
visium.breast@meta.data = merge(visium.breast@meta.data, loupe.data, by = "row.names")
row.names(visium.breast@meta.data) = visium.breast@meta.data$Row.names
Idents(visium.breast) = visium.breast@meta.data$layer_annotation
tumor.bdy = subset(visium.breast, nCount_Spatial > 100, idents = c("Bdy", "Layer"))

### Visualization of our data --------
Idents(tumor.bdy) = tumor.bdy$first_type
names(colors) = c("B cell", "CD4+ T cell", "CD8+ T cell", "DCIS #1", "DCIS #2", "Mixed", "Myeloid cell", "Plasma cell", "Stromal cell", "Invasive")
SpatialDimPlot(tumor.bdy, label = F, cols = colors)

### Prepare input data for CellChat analysis --------
data.input = GetAssayData(tumor.bdy, slot = "data", assay = "SCT")
meta = data.frame(labels = Idents(tumor.bdy), row.names = names(Idents(tumor.bdy)))
unique(meta$labels) # check the cell labels

### Load spatial imaging information --------
spatial.locs = GetTissueCoordinates(tumor.bdy, scale = NULL, cols = c("imagerow", "imagecol")) 
scale.factors = jsonlite::fromJSON(txt = "./Raw_file/visium/spatial/scalefactors_json.json")
scale.factors = list(spot.diameter = 65, 
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


## Set the ligand-receptor interaction database --------

### Load CellChat database --------
CellChatDB = CellChatDB.human
cellchat.gene = as.data.frame(CellChatDB.human$geneInfo$Symbol)
colnames(cellchat.gene) = "gene"

### Load Xenium gene panel --------
xenium.gene = read.csv("./Raw_file/xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/")

overlap.gene = merge(xenium.gene, cellchat.gene, by = "gene")
CellChatDB$interaction = CellChatDB$interaction[CellChatDB$interaction$ligand %in% overlap.gene$gene &
                                                  CellChatDB$interaction$receptor %in% overlap.gene$gene,]
cellchat@DB = CellChatDB


## Preprocessing --------

### Preprocessing the expression data for cell-cell communication analysis --------
cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)


## Infer cell-cell communication network --------

### Compute the communication probability and infer cellular communication network --------
cellchat = computeCommunProb(cellchat, type = "triMean", distance.use = TRUE, interaction.length = 200, scale.distance = 0.1)
cellchat@meta$labels %>% table() %>% sort()
cellchat = filterCommunication(cellchat, min.cells = 10)

### Infer the cell-cell communication at a signaling pathway level --------
cellchat = computeCommunProbPathway(cellchat)

### Calculate the aggregated cell-cell communication network --------
cellchat = aggregateNet(cellchat)

#saveRDS(cellchat, file = './R_object/Bioinfo2023_Breast_CellChat.rds')
#cellchat = readRDS('./R_object/Bioinfo2023_Breast_CellChat.rds')


## Visualization --------

### Visualization of the aggregated cell-cell communication network --------
groupSize = as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),weight.scale = T, 
                 sources.use =  c("Stromal cell", "Myeloid cell"), 
                 targets.use = c("DCIS #1", "DCIS #2", "Invasive"),
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, 
                 sources.use =  c("Stromal cell", "Myeloid cell"), 
                 targets.use = c("DCIS #1", "DCIS #2", "Invasive"),
                 title.name = "Interaction weights/strength")

### Identify ligand-receptor pairs between cell types --------
CellChat::netVisual_bubble(cellchat, 
                           sources.use =  c("Stromal cell", "Myeloid cell"), 
                           targets.use = c("DCIS #1", "DCIS #2", "Invasive"), 
                           remove.isolate = F, angle.x = 90, thresh = 0.05,
                           sort.by.target = T, font.size = 15) + coord_flip()

### Compute the network centrality scores --------
cellchat = netAnalysis_computeCentrality(cellchat, net.name = "CXCL12-CXCR4")

par(mfrow = c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = "CXCL12-CXCR4", width = 8, height = 2.5, font.size = 10)


# Chapter 4. Cell-cell interaction: Cellchat --------

## Load the Xenium data --------
s1r1 = LoadXenium('./Raw_file/xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs', fov = 'fov')

## Preprocessing --------
# remove cells with 0 counts
s1r1 = subset(s1r1, subset = nCount_Xenium > 0)

# Add metadata
s1r1@meta.data$cells = 'cells'

## Visualize the expression level of CXCL12 and CXCR4 --------
# Plot the positions of CXCL12 and CXCR4
ImageDimPlot(s1r1, fov = "fov", molecules = c("CXCL12", "CXCR4"), group.by = 'cells', nmols = 20000)

# Visualize the expression level of CXCL12 and CXCR4
ImageFeaturePlot(s1r1, features = c("CXCL12", "CXCR4"), max.cutoff = c(15, 3), size = 0.5, cols = c("white", "red"))

## Visualize the expression level of PTN and SDC4 --------
# Plot the positions of PTN and SDC4
ImageDimPlot(s1r1, fov = "fov", molecules = c("PTN", "SDC4"), group.by = 'cells', nmols = 20000)

# Visualize the expression level of PTN and SDC4
ImageFeaturePlot(s1r1, features = c("PTN", "SDC4"), max.cutoff = c(8, 8), size = 0.5, cols = c("white", "red"))

# Increase your RAM usage (8GB)
options(future.globals.maxSize = 8000 * 1024^2)

# Define cropped area
cropped.coords = Crop(s1r1[["fov"]], x = c(3850, 4900), y = c(6150, 7000), coords = "plot")
s1r1[["zoom"]] = cropped.coords

# Visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(s1r1[["zoom"]]) = "segmentation"
ImageDimPlot(s1r1, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("PTN", "SDC4"), nmols = 10000, group.by = 'cells')

#saveRDS(s1r1, './R_object/Bioinfo2023_Breast_xenium.rds')

