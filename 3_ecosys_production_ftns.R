# -----------------------------------------------------------------------
# Calculation ecosystem service capacity
# -----------------------------------------------------------------------

# Mary Donovan

library(raster)


# Read in rasters ---------------------------------------------------------
home <- getwd()
setwd(file.path(home,'Spatial_Predictions','MauiNui'))

list.files(getwd())
names <- list.files(getwd(),pattern = "\\.tif$")
for(i in names){assign(unlist(strsplit(i,"[.]")),raster(i))} 
setwd(home)

habcomplex <- raster(file.path(home,'Predictor_Geotiffs','MauiNui','ReefFishCommunities_EnvironmentalPredictors_SeafloorTopography_MauiNui_Slope_RateofChange_Mean.tif'))

# Scale 0 to 1 ------------------------------------------------------------
coral <- (MauiNui_Mean_Predicted_Coral-cellStats(MauiNui_Mean_Predicted_Coral,"min"))/(cellStats(MauiNui_Mean_Predicted_Coral,"max")-cellStats(MauiNui_Mean_Predicted_Coral,"min")); summary(coral)

macro <- 1-(MauiNui_Mean_Predicted_Macro-cellStats(MauiNui_Mean_Predicted_Macro,"min"))/(cellStats(MauiNui_Mean_Predicted_Macro,"max")-cellStats(MauiNui_Mean_Predicted_Macro,"min")); summary(macro)

fishrich <- (MauiNui_Mean_Predicted_Fish_Richness-cellStats(MauiNui_Mean_Predicted_Fish_Richness,"min"))/(cellStats(MauiNui_Mean_Predicted_Fish_Richness,"max")-cellStats(MauiNui_Mean_Predicted_Fish_Richness,"min")); summary(fishrich)

coralrich <- (MauiNui_Mean_Predicted_Coral_rich-cellStats(MauiNui_Mean_Predicted_Coral_rich,"min"))/(cellStats(MauiNui_Mean_Predicted_Coral_rich,"max")-cellStats(MauiNui_Mean_Predicted_Coral_rich,"min")); summary(coralrich)

fishdens <- (MauiNui_Mean_Predicted_Fish_dens-cellStats(MauiNui_Mean_Predicted_Fish_dens,"min"))/(cellStats(MauiNui_Mean_Predicted_Fish_dens,"max")-cellStats(MauiNui_Mean_Predicted_Fish_dens,"min")); summary(fishdens)


# Create raster brick -----------------------------------------------------
health.brick <- brick(coral,macro,fishrich,coralrich,fishdens)
recreation.brick <- brick(coral,fishdens,fishrich,habcomplex)

# Raster math for EPF -----------------------------------------------------
health.calc <- overlay(health.brick, fun=function(a,b,c,d,e) 0.3*a+0.15*b+0.15*c+0.20*d+0.20*e)
plot(health.calc)

recreation <- overlay(recreation.brick, fun=function(a,b,c,d) 1*a+1*b+1*c+1*d)
plot(recreation)

# Export rasters ----------------------------------------------------------
writeRaster(health.calc, filename="cri_raster", format="GTiff", overwrite=TRUE)
writeRaster(recreation, filename="recreation_raster", format="GTiff", overwrite=TRUE)


