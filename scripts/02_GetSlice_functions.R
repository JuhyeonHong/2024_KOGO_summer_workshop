# Equivalent to GetSquare but for a circle instead

GetSlice <- function(Spot, SizeMicrons, BarcodeDF, PATH, CellT = NA, size = "008um") {
  path_scales <- paste0(PATH, "/binned_outputs/square_", size, "/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = path_scales)
  Scale <- (SizeMicrons * scales$spot_diameter_fullres) / as.numeric(unlist(strsplit(size, "um"))[1])
  
  Index <- match(Spot, BarcodeDF$barcode)
  Result <- vector("list", length = length(Index))
  
  for (jj in 1:length(Index)) {
    Distance <- sqrt(((BarcodeDF$imagecol - BarcodeDF$imagecol[Index[jj]])^2) + ((BarcodeDF$imagerow - BarcodeDF$imagerow[Index[jj]])^2))
    BarcodeDF$Distance <- Distance
    
    if (!is.na(CellT)) {
      ValTh <- sum(BarcodeDF$DeconvolutionLabel1[BarcodeDF$Distance < min(Scale)] == CellT, na.rm = TRUE)
      if (ValTh < 25) {
        next
      }
    }
    
    if (length(Scale) > 1) {
      Result[[jj]] <- lapply(Scale, function(X) { return(BarcodeDF$barcode[BarcodeDF$Distance < X]) })
    } else {
      Result[[jj]] <- BarcodeDF$barcode[BarcodeDF$Distance < Scale]
    }
  }
  
  if (length(Scale) > 1) {
    Rxx <- vector("list", length = length(Scale))
    names(Rxx) <- as.character(SizeMicrons)
    
    for (ii in 1:length(Scale)) {
      Rxx[[ii]] <- lapply(Result, function(X) { return(X[[ii]]) })
      Rxx[[ii]] <- unique(unlist(Rxx[[ii]]))
    }
    return(Rxx)
  } else {
    Result <- unique(unlist(Result))
    return(Result)
  }
}