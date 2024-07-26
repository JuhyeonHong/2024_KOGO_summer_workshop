# This is for users who are using Server 1
# Package CellChat version 2.1.2

# Install packages
install.packages('tidyverse')
install.packages('wesanderson')
install.packages('furrr')
install.packages('paletteer')
install.packages('distances')
install.packages('readbitmap')

# Load packages
library(Seurat)
library(scattermore)
library(tidyverse)
library(data.table)
library(wesanderson)
library(patchwork)
library(RColorBrewer)
library(furrr)
library(paletteer)
library(arrow)
library(pheatmap)
library(RColorBrewer)
library(distances)
library(spatial)


setwd("/BiO/home/user_ID/2024_KOGO_visium/")

# Load the RDS file and the spatial coordinates data
visium.hd <- readRDS("./results/CCI/visium_hd_subset_for_CCI.rds")

source("./scripts/02_GetNearbySpots_functions.R")
source("./scripts/02_calculate_distance_functions.R")
source("./scripts/02_SelectPeripheryDiscrete_functions.R")
source("./scripts/02_GetSlice_functions.R")


coords <- as.data.frame(GetTissueCoordinates(visium.hd))
coords$type <- Idents(visium.hd)
head(coords) 


# Identify Tumor spots
tumor_spots <- coords %>% filter(type == 'Tumor')
coords$label[coords$type == 'Tumor'] <- 'Tumor'
head(coords) 


# Calculate distances and identify Periphery spots using the updated units
coords <- coords %>%
  mutate(imagecol = x, imagerow = y, DeconvolutionLabel1 = type, barcode = cell)  
head(coords)


# PATH to scalefactors_json.json
PATH <- "./data/visium_hd"  # Change this to your actual path

distance_in_units <- 50


# Use the GetNearbySpots function to identify spots within 50um distance
tumor_spots <- coords %>% filter(type == "Tumor") %>% pull(barcode)
nearby_spots <- GetNearbySpots(tumor_spots, distance_in_units, coords, PATH)


# Update labels to ‘Periphery’ within 50um distance from Tumor
coords$label[coords$barcode %in% nearby_spots & coords$type != "Tumor"] <- "Periphery"


# Update labels for Tumor spots based on the presence of Periphery within the distance
for (i in 1:nrow(coords)) {
  if (coords$type[i] == "Tumor") {
    nearby_spots <- GetNearbySpots(coords$barcode[i], distance_in_units, coords, PATH)
    nearby_labels <- coords$label[coords$barcode %in% nearby_spots]
    
    if (any(nearby_labels == "Periphery")) {
      coords$label[i] <- "bdy"
    } else {
      coords$label[i] <- "Tumor"
    }
  }
}


# Update the Seurat object’s metadata with the new labels
visium.hd <- AddMetaData(visium.hd, metadata = coords$label, col.name = 'Region')

# Set the identities to the new Region column
Idents(visium.hd) <- visium.hd$Region

# Define the colors
colors <- c("Tumor" = "red", "Periphery" = "blue", "bdy" = "green")

# Generate the SpatialDimPlot with the tissue image
SpatialDimPlot(visium.hd, 
               label = FALSE, 
               repel = TRUE, 
               pt.size.factor = 20,
               cols = colors) +
  theme(legend.text = element_text(size = 15))

saveRDS(visium.hd, "./results/CCI/visium_hd_tumor_periphery.rds") # optional





### CellChat ###

library(CellChat)
library(spacexr)
library(ggpubr)
library(Cairo)
library(readbitmap)

visium.hd <- read_rds("./results/CCI/visium_hd_tumor_periphery.rds")

# Subset only tumor boundaries and layers
Idents(visium.hd) <- visium.hd@meta.data$Region
visium.hd <- subset(visium.hd, 
                    subset = nCount_Spatial > 100 & (Region == "Periphery" | Region == "bdy"))


# check the cell labels
visium.hd$Region %>% table



# Load cell type annotated visium data and set levels for the downstream analysis.
data.input <- GetAssayData(visium.hd, layer = "data", assay = "SCT")

Idents(visium.hd) <- visium.hd@meta.data$first_type

meta <- data.frame(labels = Idents(visium.hd), 
                   row.names = names(Idents(visium.hd)))

# check the cell labels
meta$labels %>% table


# Load spatial imaging information to get the spot information
spatial.locs <- GetTissueCoordinates(visium.hd, scale = NULL)[,c(1:2)] 

scale.factors <- jsonlite::fromJSON(txt = "./data/visium_hd/binned_outputs/square_008um/spatial/scalefactors_json.json")

spot.size <- 8 # the theoretical spot size (um) in 10X Visium
ratio <-  spot.size/scale.factors$spot_diameter_fullres
tol <- spot.size/2

spatial.factors <- list(spot.diameter=8,
                        ratio = ratio,
                        tol = tol)


# Create a CellChat object for the downstream analysis.
cellchat <- createCellChat(object = data.input,
                           meta = meta, 
                           group.by = "labels",
                           datatype = "spatial", 
                           coordinates = spatial.locs,
                           spatial.factors = spatial.factors)

cellchat





# Load CellChat DB
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB


# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)

# do parallel
future::plan("multisession", workers = 4)

# Identify over-expressed ligands or receptors in one cell group
cellchat <- identifyOverExpressedGenes(cellchat)

# Identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
cellchat <- identifyOverExpressedInteractions(cellchat)


# Infers the biologically significant cell-cell communication with permutation test
cellchat = computeCommunProb(cellchat,
                             type = "triMean", 
                             distance.use = TRUE, 
                             interaction.length = 200, 
                             scale.distance = 0.1)

cellchat <- filterCommunication(cellchat, min.cells = 10)


saveRDS(cellchat, file = './results/CCI/visium_hd_cellchat.rds') # (optional)





# read RDS file
cellchat <- readRDS("./results/CCI/visium_hd_cellchat.rds")

# Computes the communication probablity on signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)


# Visualization of the aggregated cell-cell communication network
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, 
                 vertex.weight = rowSums(cellchat@net$count), 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight,
                 vertex.weight = rowSums(cellchat@net$weight), 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")



# Identify ligand-receptor pairs between cell types
CellChat::netVisual_bubble(cellchat, 
                           sources.use = c("Fibroblast", "Endothelial", "Intestinal Epithelial", "Myeloid", "Neuronal","Smooth Muscle", "T cells", "B cells"),
                           targets.use = "Tumor",
                           remove.isolate = FALSE, angle.x = 90, thresh = 0.05) + coord_flip()




# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, net.name = "COL1A1-CD44")


# Visualize the centrality score
netAnalysis_signalingRole_network(cellchat, signaling = "COL1A1-CD44", 
                                  width = 8, height = 2.5, font.size = 10)

