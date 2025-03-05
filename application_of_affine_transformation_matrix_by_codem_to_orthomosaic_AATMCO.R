# Set working directory
setwd("E:/Data/Geysernaya/Lava/")



# Read and process transformation matrix from registration file
trmat <- readLines("registration_2025-02-16_17-09-34/registration.txt")
trmat <- trmat[grep("ICP REGISTRATION", trmat)+3:6]  # Extract relevant lines
clean <- gsub("[^0-9.eE+-]", " ", trmat)  # Remove non-numeric characters
clean <- gsub(" +", " ", clean)         # Collapse multiple spaces
clean <- trimws(clean)                  # Trim whitespace
split <- strsplit(clean, " ")           # Split into numeric components
numbers <- lapply(split, as.numeric)    # Convert to numeric values
matrix_4x4 <- matrix(unlist(numbers), nrow = 4, byrow = TRUE)  # Create 4x4 transformation matrix

# Load required library and input data
# install.packages("terra")
library(terra)
(dem <- rast("dem21.tif"))  # Load DEM raster
(opp <- rast("opp21.tif"))  # Load orthophoto raster

# Resolution matching options
# Option 1: Degrade orthophoto resolution to match DEM
#dem <- terra::aggregate(dem, 5)  # Reduce DEM resolution if you need (commented out)
opp <- resample(opp, dem)         # Resample orthophoto to match
# Option 2: Enhance DEM resolution to match orthophoto (commented out)
# dem <- resample(dem, opp)

# Configure tiling parameters
n_parts <- 4    # Create 4x4 grid (16 tiles). Use 1 part to do processing without raster split (if you have enough memory)
overlap <- 2    # 2-meter overlap between tiles

# Calculate raster extent
x <- ext(opp)[1:2]  # xmin/xmax
y <- ext(opp)[3:4]  # ymin/ymax

# Calculate tile dimensions
dx <- (x[2] - x[1]) / n_parts  # Tile width
dy <- (y[2] - y[1]) / n_parts  # Tile height

# Generate tile extents with overlap
tiles <- list()
for (i in 1:n_parts) {
  for (j in 1:n_parts) {
    # Calculate tile coordinates with overlap
    xmin <- x[1] + (i-1)*dx - overlap
    xmax <- x[1] + i*dx + overlap
    ymin <- y[1] + (j-1)*dy - overlap
    ymax <- y[1] + j*dy + overlap
    
    # Clamp to original raster boundaries
    xmin <- max(xmin, x[1])
    xmax <- min(xmax, x[2])
    ymin <- max(ymin, y[1])
    ymax <- min(ymax, y[2])
    
    tiles[[paste0("tile_", i, "_", j)]] <- ext(xmin, xmax, ymin, ymax)
  }
}

# Write tiles to disk
for (e in 1:length(tiles)) {
  writeRaster(crop(opp, tiles[[e]]), paste0("opp_tile_",e,".tif"), overwrite= T)
  writeRaster(crop(dem, tiles[[e]]), paste0("dem_tile_",e,".tif"), overwrite= T)
}

# Initialize registration processing
opp_registered <- list()

# Process each tile
for (i in 1:length(tiles)) {
  # Load tile data
  temp_dem <- rast(paste0("dem_tile_",i,".tif"))
  temp_opp <- rast(paste0("opp_tile_",i,".tif"))
  
  # Convert DEM to data frame
  demdf <- as.data.frame(temp_dem, xy = TRUE, na.rm = F)
  demdf <- data.frame(demdf, 1)  # Add homogeneous coordinate
  
  # Apply transformation matrix
  demdf_trans <- apply(na.omit(demdf), 1, function(x) matrix_4x4 %*% x)
  
  # Extract original colors
  oldcolours <- extract(temp_opp, na.omit(demdf)[,1:2])[,-1]
  
  # Process color channels
  demdf_trans2 <- cbind(demdf_trans[,-c(3:4)], oldcolours)
  
  # Rasterize each channel
  R <- rasterize(as.matrix(demdf_trans2[,1:2]), temp_dem, demdf_trans2[,3])
  G <- rasterize(as.matrix(demdf_trans2[,1:2]), temp_dem, demdf_trans2[,4])
  B <- rasterize(as.matrix(demdf_trans2[,1:2]), temp_dem, demdf_trans2[,5])
  alpha <- rasterize(as.matrix(demdf_trans2[,1:2]), temp_dem, demdf_trans2[,6])
  
  # Combine channels and store result
  opp_corrected <- c(R, G, B, alpha)
  opp_registered[[i]] <- opp_corrected
}

# Mosaic tiles and save result
opp_registered_mosaic <- do.call(terra::mosaic, opp_registered)
writeRaster(opp_registered_mosaic, "orthomosaic_registered.tif", NAflag = 0, overwrite = TRUE)