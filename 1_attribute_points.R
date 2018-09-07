
# -----------------------------------------------------------------------
#  Attribute points with predictor data 
# -----------------------------------------------------------------------

# Mary K. Donovan
# February 2018

rm(list=ls())
library(raster)
library(stringr)
library(plyr)


# Read in raster tiffs and attribute points -------------------------------

islands <- c('MauiNui')

for(k in 1:length(islands)){
  # 1. Read in and prepare rasters ------------------------------------------------------
  home <- getwd()
  setwd(file.path(getwd()),'predictor_rasters') 
  # this folder should be populated with the data accessed here:
  # https://www.nodc.noaa.gov/archive/arc0095/0155189/1.1/data/0-data/Marine_Biogeographic_Assessment_of_the_Main_Hawaiian_Islands/GIS_Data/GeoTiffs-Shapefiles/GIS-Data_Geotiffs-Shapefiles_Chapter-4_Fishes.zip
  
  list.files(getwd())
  rnames <- list.files(getwd(),pattern='\\.tif$')
  for(i in rnames){assign(unlist(strsplit(i,'[.]')),raster(i))}
  rnames <- str_replace_all(rnames,'.tif',"")
  
  setwd(home)
  
  # create a list of variable names
  vars <- rep(NA,length(rnames))
  for(i in 1:23){
    vars[i] <- unlist(strsplit(rnames[i],'ReefFishCommunities_EnvironmentalPredictors_'))[2]
    vars[i] <- unlist(strsplit(vars[i],'BenthicHabitatComposition_'))[2]
    vars[i] <- unlist(strsplit(vars[i],paste(islands[k],"_",sep="")))[2]
    # vars[i] <- str_replace_all(vars[i],'.tif',"")
  }
  for(i in 24:27){
    vars[i] <- unlist(strsplit(rnames[i],'ReefFishCommunities_EnvironmentalPredictors_'))[2]
    vars[i] <- unlist(strsplit(vars[i],'Geographic_'))[2]
    vars[i] <- unlist(strsplit(vars[i],paste(islands[k],"_",sep="")))[2]
    # vars[i] <- str_replace_all(vars[i],'.tif',"")
  }
  for(i in 28:35){
    vars[i] <- unlist(strsplit(rnames[i],'ReefFishCommunities_EnvironmentalPredictors_'))[2]
    vars[i] <- unlist(strsplit(vars[i],'PhysicalOceanography_'))[2]
    vars[i] <- unlist(strsplit(vars[i],paste(islands[k],"_",sep="")))[2]
    # vars[i] <- str_replace_all(vars[i],'.tif',"")
  }
  for(i in 36:length(rnames)){
    vars[i] <- unlist(strsplit(rnames[i],'ReefFishCommunities_EnvironmentalPredictors_'))[2]
    vars[i] <- unlist(strsplit(vars[i],'SeafloorTopography_'))[2]
    vars[i] <- unlist(strsplit(vars[i],paste(islands[k],"_",sep="")))[2]
    # vars[i] <- str_replace_all(vars[i],'.tif',"")
  }
  vars
  
  # 2. Read in point data ---------------------------------------------------
  setwd(file.path(getwd()),'response_data')
  resp <- read.csv('response_4marxan_030518.csv')
  pts <- read.csv('coords_transform_proj.csv')
  setwd(home)
  
  resp <- join(resp,pts[c('coord.join4','x_prj','y_prj')],by='coord.join4')
  resp <- subset(resp, Island.grp==islands[k])
  
  # create spatial object from points
  xy <- resp[c('x_prj','y_prj')]
  xy <- SpatialPointsDataFrame(coords=xy,data=resp)
  plot(xy)
  
  # 3. extract data ---------------------------------------------------------
  expl <- data.frame(coord.join4=resp$coord.join4)
  for(i in 1:length(rnames)){
    expl <- cbind(expl, extract(get(rnames[i]),xy))
    colnames(expl)[i+1] <- vars[i]
  }
  str(expl)
  expl <- join(resp,expl,by='coord.join4')
  expl <- expl[!is.na(expl['Rugosity_Mean']),]
  
  # 4. export data ----------------------------------------------------------
  setwd(file.path(getwd(),"Input_CSVs"))
  write.csv(expl, paste(islands[k],'_Data','.csv',sep=""),row.names=F)
  setwd(home)
  
  # 5. export new tif with short name ------------------------------------------
  dir.create(file.path(getwd(),"Predictor_geotiffs_new", islands[k]), showWarning=TRUE, recursive=TRUE)
  setwd(pastefile.path(getwd(),"Predictor_geotiffs_new", islands[k]))
  for(i in 1:length(rnames)){
    writeRaster(get(rnames[i]),paste(vars[i]),format='GTiff')
  }
  writeRaster(get(paste("ReefFishCommunities_EnvironmentalPredictors_Geographic_",islands[k],"_Latitude",sep="")),'y_prj',format='GTiff')
  writeRaster(get(paste("ReefFishCommunities_EnvironmentalPredictors_Geographic_",islands[k],"_Longitude",sep="")),'x_prj',format='GTiff')
  setwd(home)
}


# clean up column names ---------------------------------------------------

setwd(file.path(getwd(),'Input_CSVs'))
list.files()

data.names <- unique(names(read.csv('MauiNui/MauiNui_Data.csv')))
data.names <- data.names[c(28:length(data.names))]

fpred.names <- c(
  as.character(read.csv(paste(home,'/Input_CSVs/MauiNui/MauiNui_Final_Predictors.csv',sep=""),header=F)$V1))
fpred.names <- unique(fpred.names)

comb.names <- data.frame(data.name = data.names, fpred.names=NA)
comb.names$fpred.names[comb.names$data.name=="x_prj"] <- 'x'
comb.names$fpred.names[comb.names$data.name=="y_prj"] <- 'y'
comb.names$fpred.names[comb.names$data.name=='Contiguity_Index_Mean'] <- 'con_mn'
comb.names$fpred.names[comb.names$data.name=='Contiguity_Index_StdDev'] <- 'con_sd'
# comb.names$fpred.names[comb.names$data.name=='Edge_Density_CCA'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Edge_Density_Coral'] <- 'ed_cor'
comb.names$fpred.names[comb.names$data.name=='Edge_Density_Macroalgae'] <- 'ed_mac'
comb.names$fpred.names[comb.names$data.name=='Edge_Density_Soft_Bottom'] <- 'ed_sof'
comb.names$fpred.names[comb.names$data.name=='Edge_Density_Turf_Algae'] <- 'ed_tur'
comb.names$fpred.names[comb.names$data.name=='Fractal_Dimension_Mean'] <- 'fra_mn'
# comb.names$fpred.names[comb.names$data.name=='Fractal_Dimension_StdDev'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Patch_Shape_Index_CCA'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Patch_Shape_Index_Coral'] <- 'shp_cor'
# comb.names$fpred.names[comb.names$data.name=='Patch_Shape_Index_Macroalgae'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Patch_Shape_Index_Soft_Bottom'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Patch_Shape_Index_Turf_Algae'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Percent_CCA'] <- 'pls_cca'
comb.names$fpred.names[comb.names$data.name=='Percent_Coral'] <- 'pls_cor'
comb.names$fpred.names[comb.names$data.name=='Percent_Macroalgae'] <- 'pls_mac'
comb.names$fpred.names[comb.names$data.name=='Percent_Soft_Bottom'] <- 'pls_sof'
comb.names$fpred.names[comb.names$data.name=='Percent_Turf_Algae'] <- 'pls_tur'
comb.names$fpred.names[comb.names$data.name=='Proximity_Index_Mean'] <- 'pro_mn'
# comb.names$fpred.names[comb.names$data.name=='Proximity_Index_StdDev'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Shannons_Diversity_Index'] <- 'shdi'
# comb.names$fpred.names[comb.names$data.name=='Shannons_Evenness_Index'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Distance_to_Shore'] <- 'dist2shore'
# # comb.names$fpred.names[comb.names$data.name=='Latitude'] <- 'n_a'
# # comb.names$fpred.names[comb.names$data.name=='Longitude'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Proximity_to_Human_Population'] <- 'dist2pop'
# comb.names$fpred.names[comb.names$data.name=='Wave_Height_90th_Percentile'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Wave_Height_95th_Percentile'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Wave_Height_Mean'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Wave_Height_StdDev'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Wave_Power_90th_Percentile'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Wave_Power_95th_Percentile'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Wave_Power_Mean'] <- 'wav_Pmn'
comb.names$fpred.names[comb.names$data.name=='Wave_Power_StdDev'] <- 'wav_Psd'
comb.names$fpred.names[comb.names$data.name=='Aspect_Circular_Mean_Cosine'] <- 'cosasp'
comb.names$fpred.names[comb.names$data.name=='Aspect_Circular_Mean_Sine'] <- 'sinasp'
comb.names$fpred.names[comb.names$data.name=='Aspect_Circular_StdDev'] <- 'asp_sd'
comb.names$fpred.names[comb.names$data.name=='Bathymetric_Position_Index_120m'] <- 'bpi_120'
comb.names$fpred.names[comb.names$data.name=='Bathymetric_Position_Index_240m'] <- 'bpi_240'
comb.names$fpred.names[comb.names$data.name=='Bathymetric_Position_Index_60m'] <- 'bpi'
comb.names$fpred.names[comb.names$data.name=='Depth_Mean'] <- 'bathy'
# comb.names$fpred.names[comb.names$data.name=='Depth_StdDev'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Planform_Curvature_Mean'] <- 'plcurv'
comb.names$fpred.names[comb.names$data.name=='Planform_Curvature_StdDev'] <- 'plcurv_sd'
comb.names$fpred.names[comb.names$data.name=='Profile_Curvature_Mean'] <- 'prcurv'
# comb.names$fpred.names[comb.names$data.name=='Profile_Curvature_StdDev'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Rugosity_Max_240m'] <- 'rug_240'
# comb.names$fpred.names[comb.names$data.name=='Rugosity_Mean'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Rugosity_StdDev'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Slope_Max_240m'] <- 'n_a'
comb.names$fpred.names[comb.names$data.name=='Slope_Mean'] <- 'slp'
comb.names$fpred.names[comb.names$data.name=='Slope_RateofChange_Max_240m'] <- 'slp_240'
comb.names$fpred.names[comb.names$data.name=='Slope_RateofChange_Mean'] <- 'slpslp'
# comb.names$fpred.names[comb.names$data.name=='Slope_RateofChange_StdDev'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Slope_StdDev'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Terrain_Ruggedness_Max_240m'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Terrain_Ruggedness_Mean'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Terrain_Ruggedness_StdDev'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Total_Curvature_Max_240m'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Total_Curvature_Mean'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='Total_Curvature_StdDev'] <- 'n_a'
# comb.names$fpred.names[comb.names$data.name=='RateofChange_StdDev'] <- 'n_a'
# write.csv(comb.names, "combined_names.csv", row.names=F)


# read in each file and create a new ‘final predictor’ file  --------

island <- 'MauiNui'
fpred.names <- read.csv(file.path(getwd(),'Input_CSVs',island,paste(island,'_Final_Predictors.csv',sep="")),header=F)
colnames(fpred.names) <- 'fpred.names'
fpred.names <- join(fpred.names,comb.names,by="fpred.names")
write.csv(fpred.names, paste(island,'/',island,'_Final_Predictors.csv',sep=""))


