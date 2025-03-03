library(raster)
library(caret)
library(tidyverse)
library(whitebox)
library(sf)
library(terra)

wbt_init()
# Visually interpreted bedrock presence/absence
yn <- read_csv("ConePond/Cone_training_BOSS.csv") 

sppoints <- st_read("ConePond/ConePrandompoints1.shp")

sppoints <- sppoints |> left_join(yn, by = "pointnum") |> drop_na(bdr)

#coordinates(plots) <- ~easting+northing
#proj4string(plots) <- CRS("+init=epsg:26919")

# Make BDR a factor
#n$bdr_new = with(n, ifelse(bdr == 0, 'No', 
#                    ifelse(bdr == 1, 'Yes', 0)))
#n$bdr_new <- as.factor(n$bdr_new)
sppoints <- sppoints |> mutate(bdr = as.factor(bdr))



# Import topo metrics
mxslope <- raster("ConePond/slope_rad_5m.tif")
hbtpi15 <- raster("ConePond/tpi15saga.tif")
hbtpi200 <- raster("ConePond/tpi200saga.tif")

#Stack and rename raster grids
SPC <- stack(mxslope, hbtpi15, hbtpi200)
names(SPC) <- c('maxslope', 'tpi15', 'tpi200')

#extract topo metric values to n_factor points
#make n spatial
#n_spat <- dplyr::select(n, lon = x_coord, lat = y_coord) %>%
#  SpatialPoints(proj4string = CRS("+proj=longlat +datum=WGS84")) %>%
#  spTransform(CRS("+init=epsg:32610"))

extracted <- terra::extract(x = SPC, #raster 
                             y = sppoints, #points
                             method = "simple")

#extracted col names maxslope, tpi15, tpi200
n <- cbind(sppoints, extracted)

#make bdr a factor with yes and no as factors
#this makes 0 no and 1 yes
levels(n$bdr) <- c("no", "yes")

# Prepare training scheme
control <- trainControl(method="cv", number=5, classProbs=TRUE)

# Train the GAM model
set.seed(7)
modelGAMs <- train(bdr ~ maxslope + tpi15 + tpi200,  metric = "Accuracy",
                   data = n, method = "gam", trControl = control)

# Predict the model (this takes a long time)
map.GAM.c <- predict(SPC, modelGAMs, type = "prob", dataType = "INT1U",
                     filename = "HBValley/20Feb2025_boss.tif",
                     format = "GTiff", overwrite = T, progress = "text")
map.GAC.c.BU <- map.GAM.c
#classify BR areas 2 = BR 1 = Not BR
BRthreshold <- 0.45
map.GAM.c[map.GAM.c >= BRthreshold] <- 1
map.GAM.c[map.GAM.c < BRthreshold] <- 2

# Export the final predictions into a raster
writeRaster(map.GAM.c, filename = "ConePond/bossBRprediction.tif",
            format = "GTiff", progress="text", overwrite = TRUE) 
writeRaster(map.GAC.c.BU, filename = "ConePond/bossBRprobabilities.tif",
            format = "GTiff", progress="text", overwrite = TRUE) 

#clean small bits of BR
wbt_majority_filter(input = "ConePond/bossBRprediction.tif",
                      output = "ConePond/bossBRprediction_MF5.tif",
                      filterx = 5,
                      filtery = 5)
plot(raster("ConePond/bossBRprediction_MF5.tif"))

plot(raster("ConePond/bossBRprediction.tif"))

plot(raster("Original_Metrics/gamroumxvrmt15t25t30t100t140/bdr100k.tif"))

plot(raster("HBValley/bossBRprediction_sagatpi.tif"))

