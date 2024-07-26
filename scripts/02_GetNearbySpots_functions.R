# Function to get nearby spots within a given distance

GetNearbySpots <- function(Spot, SizeMicrons, BarcodeDF, PATH, size = "008um") {
  path_scales <- paste0(PATH, "/binned_outputs/square_", size, "/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = path_scales)
  Scale <- (SizeMicrons * scales$spot_diameter_fullres) / as.numeric(unlist(strsplit(size, "um"))[1])
  
  Index <- match(Spot, BarcodeDF$barcode)
  Result <- vector("list", length = length(Index))
  
  for (jj in 1:length(Index)) {
    Distance <- sqrt(((BarcodeDF$imagecol - BarcodeDF$imagecol[Index[jj]])^2) + ((BarcodeDF$imagerow - BarcodeDF$imagerow[Index[jj]])^2))
    BarcodeDF$Distance <- Distance
    Result[[jj]] <- BarcodeDF$barcode[BarcodeDF$Distance < Scale]
  }
  Result <- unique(unlist(Result))
  return(Result)
}