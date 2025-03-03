library(raster)
library(caret)
library(tidyverse)
library(whitebox)
wbt_init()
# Visually interpreted bedrock presence/absence
#n <- read.csv("~/R/bedrock2020/n_factor.csv")
n <- read.csv("n_factor.csv")

#coordinates(plots) <- ~easting+northing
#proj4string(plots) <- CRS("+init=epsg:26919")

# Make BDR a factor
n$bdr_new = with(n, ifelse(bdr == 0, 'No', 
                                ifelse(bdr == 1, 'Yes', 0)))
n$bdr_new <- as.factor(n$bdr_new)


# Import topo metrics
mxslope <- raster("HBValley/slope_rad_5m.tif")
hbtpi15 <- raster("HBValley/tpi15saga.tif")
hbtpi200 <- raster("HBValley/tpi200saga.tif")

#Stack and rename raster grids
SPC <- stack(mxslope, hbtpi15, hbtpi200)
names(SPC) <- c('maxslope', 'tpi15', 'tpi200')

#extract topo metric values to n_factor points
#make n spatial
n_spat <- dplyr::select(n, lon = x_coord, lat = y_coord) %>%
  SpatialPoints(proj4string = CRS("+proj=longlat +datum=WGS84")) %>%
  spTransform(CRS("+init=epsg:32610"))

extracted <- raster::extract(x = SPC, #raster 
                             y = n_spat, #points
                             method = "simple")

#extracted col names maxslope, tpi15, tpi200
n <- cbind(n, extracted)

# Prepare training scheme
control <- trainControl(method="cv", number=5, classProbs=TRUE)

# Train the GAM model
set.seed(7)
modelGAMs <- train(bdr_new ~ maxslope + tpi15 + tpi200,  metric = "Accuracy",
                   data = n, method = "gam", trControl = control)

# Predict the model (this takes a long time)
map.GAM.c <- predict(SPC, modelGAMs, type = "prob", dataType = "INT1U",
                     filename = "HBValley/hbmxslt15t200gam_new9.tif",
                     format = "GTiff", overwrite = T, progress = "text")
map.GAC.c.BU <- map.GAM.c
#classify BR areas 2 = BR 1 = Not BR
BRthreshold <- 0.45
map.GAM.c[map.GAM.c >= BRthreshold] <- 1
map.GAM.c[map.GAM.c < BRthreshold] <- 2

# Export the final predictions into a raster
writeRaster(map.GAM.c, filename = "HBValley/bossBRprediction_sagatpi.tif",
            format = "GTiff", progress="text", overwrite = TRUE) 
writeRaster(map.GAC.c.BU, filename = "HBValley/bossBRprobs_sagatpi.tif",
            format = "GTiff", progress="text", overwrite = TRUE) 

#clean small bits of BR
wbt_majority_filter(input = "HBValley/bossBRprediction_sagatpi.tif",
                      output = "HBValley/bossBRprediction_sagatpi_MF5.tif",
                      filterx = 5,
                      filtery = 5)
plot(raster("HBValley/bossBRprediction_sagatpi_MF5.tif"))

plot(raster("Original_Metrics/gamroumxvrmt15t25t30t100t140/bdr100k.tif"))

plot(raster("HBValley/bossBRprediction_sagatpi.tif"))

