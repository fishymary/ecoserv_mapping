This repository includes data and analysis scripts to accompany:

# Mapping supply and demand of ecosystem services for effective spatial prioritization

### Authors: Mary K Donovan, Kirsten LL Oleson, Joey Lecky, Jade MS Delevaux, Alan M Friedlander, Matthew Poti, Kostantinos A Stamoulis, Susan Yee
### Journal: Nature Sustainability (submitted)
### Link: tbd

-----

### Description:
This work investigates components of ecosystem service flow from the capacity of ecosysytem to produce services to the ultimate benefits to society. We mapped ecological capacity of coral reefs in Hawaii for multiple services by translating ecological attributes (e.g., fish biomass) into ecosystem service capacity using ecosystem production functions. We then compared the spatial footprint of that capacity with human pressures and human demand to elucidate the intersection points for management, and to provide a framework for spatial prioritization of management actions.

### Contents:
#### Scripts:
* **1_attribute_points.R:** R script that takes spatial points and attributes data from predictor rasters
* **2_BRT_fit_predict.R:** R script that conducts BRT analysis and produces raster predictions
* **3_ecosys_production_ftns.R:** R script that translates ecological rasters into ecosystem service capacity rasters
* **4_quantile_overlaps.R:** R script that compares overlap between capacity, pressures, and demand of ecosystem services

#### Data:
Predictor shapefiles are available at <http://bit.ly/ch4_geotiffs>. Fish response variables are available at <http://bit.ly/ch4_geodatabase>. Benthic response variables are available at <http://bit.ly/ch3_geodatabase>.
