library(devtools)
install_version("maptools", version = "1.1-8", repos = "http://cran.us.r-project.org")
install_github("envirometrix/landmap")

url <- "https://download.r-forge.r-project.org/src/contrib/rgdal_1.6-7.tar.gz"
install.packages(url, type="source", repos=NULL)

#Load topo metric packages
set.seed(100)
library(tidyverse)
library(terra)
library(raster)
library(sp)
library(sf)
#library(rgeos)
#library(rgdal)
library(geoR)
#Load Ensemble Machine Learning packages
#library(landmap) #not on CRAN, install using devtool install_github
#library(plotKML) #not on CRAN
library(viridis) 
library(glmnet)
library(xgboost)
library(kernlab)
library(deepnet)
library(mlr)
library(rpart)
library(nnet)
library(tmap)

data_dir <- "mosaic_HB_CP_BP/"

#Add raster files (always from the main project folder otherwise things start to fail)
hbmaxslope <- raster(paste0(data_dir,"slope_per_qgis.tif"))
hbuaab     <- raster(paste0(data_dir,"uaab_norm2.tif")) #bad
twid       <- raster(paste0(data_dir,"hydem5m_TWId.tif")) #twi, check
mrvbf      <- raster(paste0(data_dir,"mrvbf.tif")) #check
tpi20      <- raster(paste0(data_dir,"tpi20saga.tif"))
tpi100     <- raster(paste0(data_dir,"tpi100saga.tif"))
tpi200     <- raster(paste0(data_dir,"tpi200saga.tif"))
#tpi500     <- raster(paste0(data_dir,"tpi500saga.tif"))
#tpi2000    <- raster(paste0(data_dir,"tpi2000saga.tif"))
EDb        <- raster(paste0(data_dir, "EDb.tif"))
#hbedbgam45 <- raster(paste0(data_dir,"hbedbgam45rec.tif"))
#bedrock    <- raster(paste0(data_dir, "bossBRprediction.tif"))

#NDVI       <- raster(paste0(data_dir,"S2aNDVI.tif"))
#twid <- raster("hydem5m_TWId.tif")


plot(hbmaxslope)
plot(hbuaab    )
plot(twid      )
plot(mrvbf     )
plot(tpi20     )
plot(tpi100    )
plot(tpi200    )
plot(EDb       )


#Load pedon shapefile

###################   JP    ###########################################################

#extract raster values to pedon points ##
#######################################################################################


# Add pedon shapefile
plots <- st_read("HBValley/HBEF_NEW_biscuit_Pedon_LAST/HBEF_NEW_biscuit_Pedon_LAST.shp")
plots <- as.data.frame(plots)

coordinates(plots) <- ~easting+northing
proj4string(plots) <- CRS("+init=epsg:26919")
#crs(plots) <- "+init=epsg:26919"
plots$hpu <- as.factor(plots$hpu)
levels(plots$hpu)

#OGbedrock <- raster("Original_Metrics/bedrock/bdr100k.tif")

#check points and rasters align
tmap_mode("view")
tm_shape((hbuaab))+
  tm_raster(style = "cont")+
  tm_shape(plots)+
  tm_dots(col = "red")

######################################################################################
# https://opengeohub.github.io/spatial-sampling-ml/generating-spatial-sampling.html
# Stack rasters and use naming convention from fields in pedon shapefile (this could be cleaned up in the future!)
SPC <- stack(tpi20, tpi100, tpi200, mrvbf, hbuaab, EDb, twid)
names(SPC) <- c('tpi20', 'tpi100', 'tpi200m', 'mrvbf', 'hbuaab', 'EDb', 'twid')

# PC transformation of a subset and then all covariates (if needed)
# spdf_spc = landmap::spc(spdf_all_layers)
spdf_all_layers <- as(SPC, "SpatialPixelsDataFrame")

############################# ENSEMBLE MACHING LEARNING  #############################
#sl.c <- c("classif.ranger", "classif.xgboost", "classif.nnTrain")
#SL.library <- c("classif.ranger", "classif.xgboost", "classif.multinom")
#SL.library <- c("classif.ranger", "classif.xgboost", "classif.multinom")

SL.library <- c("classif.ranger", "classif.svm", "classif.multinom")

mC <- train.spLearner(plots["hpu"], 
                      covariates = spdf_all_layers[,c('tpi20', 'tpi100', 'tpi200m',
                                                    'mrvbf', 'hbuaab', 'EDb',
                                                    'twid')],
                      SL.library = SL.library, 
                      super.learner = "classif.glmnet", 
                      parallel = FALSE, 
                      oblique.coords = TRUE)

saveRDS(mC, paste0(data_dir, "model", Sys.Date(), ".rds"))

# Predict. This will take a long time!
hb.hpu <- predict(mC)

r <- raster(hb.hpu$pred, "response")

writeRaster(r, 
            paste0(data_dir, "modelout", Sys.Date(), ".tif"), 
            overwrite=TRUE)

mC <- readRDS("HBValley/model2024-01-31.rds")

test <- raster(paste0(data_dir, "Model_Verified.tif"))
#########################################################################################################
# CONFUSION MATRIX
#########################################################################################################

newdata = mC@vgmModel$observations@data
sel.e = complete.cases(newdata[,mC@spModel$features])
newdata = newdata[sel.e, mC@spModel$features]
pred = predict(mC@spModel, newdata=newdata)
pred$data$truth = mC@vgmModel$observations@data[sel.e, "hpu"]
print(calculateConfusionMatrix(pred))

#print(calculateConfusionMatrix(pred$data$response))
table(as.data.frame(pred)$truth)
table(as.data.frame(pred)$response)

performance(pred, measures = mmce)
performance(pred, measures = list(mmce, acc))

#########################################################################################################
# PROBABILITIES
#########################################################################################################

pEE <- raster(hb.hpu$pred, "prob.E")
pEH <- raster(hb.hpu$pred, "prob.H")
pEO <- raster(hb.hpu$pred, "prob.O")
pEI <- raster(hb.hpu$pred, "prob.I")
pEBhs <- raster(hb.hpu$pred, "prob.Bhs")
pET <- raster(hb.hpu$pred, "prob.T")
pEBh <- raster(hb.hpu$pred, "prob.Bh")
pEBi <- raster(hb.hpu$pred, "prob.Bi")

writeRaster(pEE,  paste0(data_dir,'prob_E_', Sys.Date(),'.tif'))
writeRaster(pEH,  paste0(data_dir,'prob_H_', Sys.Date(),'.tif'))
writeRaster(pEO,  paste0(data_dir,'prob_O_', Sys.Date(),'.tif'))
writeRaster(pEI,  paste0(data_dir,'prob_I_', Sys.Date(),'.tif'))
writeRaster(pEBhs,paste0(data_dir,'prob_Bhs_', Sys.Date(),'.tif'))
writeRaster(pET,  paste0(data_dir,'prob_T_', Sys.Date(),'.tif'))
writeRaster(pEBh, paste0(data_dir,'prob_Bh_', Sys.Date(),'.tif'))
writeRaster(pEBi, paste0(data_dir,'prob_Bi_', Sys.Date(),'.tif'))
