# Function to select the barcodes that are within a given distance from a given cell type (cluster)

SelectPeripheryDiscrete <- function(bcs, CellType, distance = 50, PATH) {
  SelectedBCs <- bcs %>% filter(DeconvolutionLabel1 == CellType)
  Result <- GetSlice(SelectedBCs$barcode, distance, bcs, PATH, CellT = CellType)
  
  if (length(distance) > 1) {
    for (jj in 1:length(Result)) {
      Result[[jj]] <- Result[[jj]][!(Result[[jj]] %in% SelectedBCs$barcode)]
    }
    return(Result)
  } else {
    Result <- Result[!(Result %in% SelectedBCs$barcode)]
    return(Result)
  }
}