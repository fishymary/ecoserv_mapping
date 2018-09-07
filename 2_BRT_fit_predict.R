
# -----------------------------------------------------------------------
# Mapping reef metrics for main Hawaiian Islands
# -----------------------------------------------------------------------

# Methods follow Stamoulis et al. (2017) from Biogeographic Assessment of Main Hawaiian Islands - Chapter 4
# Code was developed by Matt Poti, Jade Delevaux, Kosta Stamoulis, and Mary Donovan in 2014-2016, with additional edits by Mary Donovan 2018

# This version was last modified in March 2018

island <- "MauiNui"   # island name/code
response_vars_all <- c('num','Resource','rich','cultsp.num','Coral','Macro','coral.rich')
output_name_all <- c('fish_dens','resource','fish_rich','cult_spp','coral','macro','coral_rich')
transformation_all <- c('square.root','fourth.root','square.root','fourth.root','square.root','fourth.root','square.root')


for(z in 1:length(response_vars_all)){
  response_var <- response_vars_all[z]		# response variable attribute field name
  output_name <- output_name_all[z]	# generic name of response variable for use in output file names
  transformation <- transformation_all[z] # options include "log" (log+1), "square.root", and "fourth.root"
  
  # Initialization ----------------------------------------------------------
  # rm(list=ls())
  tot.start <- Sys.time()
  
  library(dismo) # version 1.1.4
  library(parallel) # version 3.4.1
  library(geosphere) # version 1.5.5
  library(ape) # version 4.1
  
  # home <- getwd()
  home <- '/Users/mary/Dropbox/My backup/backup081612/Castle Fnd - healthy reefs/spm_formaxan'
  setwd(home)
  
  # Define Environments/Variables -------------------------------------------
  training_holdout <- 30   # percentage of data to be withheld for model validation
  num_bootstraps <- 25      # number of bootstrap samples to use to generate spatial predictions
  nCPUs <- detectCores()-1		# detect number of processing cores for parallel processing; can instead manually choose number
  set.seed(28)	
  
  # Define functions --------------------------------------------------------
  #Function to apply data transformation  --- Options include logarithm, square root, and fourth root
  Transform.Data <- function(d, transformation) {
    if (transformation == "log") {
      d_transformed <- log(d + 1)
    } else if (transformation == "square.root") {
      d_transformed <- sqrt(d)
    } else if (transformation == "fourth.root") {
      d_transformed <- d^(1/4)    
    }
  }
  
  #Function to fit BRT model with optimal number of boosting trees for a given set of model tuning parameters
  Fit.BRT.Model <- function(i) {
    library(dismo)
    gbm.step(data=calib_data,
             gbm.x = which(names(calib_data) %in% predictor_set),
             gbm.y = grep(response_var_transformed,names(calib_data)),
             tree.complexity=tuning_parameters[i,3],
             learning.rate=tuning_parameters[i,1],
             bag.fraction=tuning_parameters[i,2],
             n.folds=10,
             family="gaussian",
             plot.main=FALSE)
  }
  
  #Function to fit BRT model with fixed number of boosting trees for a given set of model tuning parameters using a bootstrap sample of the calibration data
  Fit.BRT.Model.Fixed.Bootstrap <- function(i) {
    library(dismo)
    gbm.fixed(data=calib_data[bootstrap_samples[[i]],],
              gbm.x = final_vars_ix,
              gbm.y = grep(response_var_transformed,names(calib_data)),
              tree.complexity=gbm_model_simp$gbm.call$tree.complexity,
              learning.rate=gbm_model_simp$gbm.call$learning.rate,
              bag.fraction=gbm_model_simp$gbm.call$bag.fraction,
              family="gaussian",
              n.trees=gbm_model_simp$gbm.call$best.trees)
  }
  
  #Function to generate a spatial (raster) prediction from a boosted regression tree model
  Spatial.Prediction <- function(i) {
    library(gbm)
    library(raster)
    library(rgdal)
    gbm_model <- gbm_models_final[[i]]
    predict(raster_stack, gbm_model, n.trees=gbm_model$gbm.call$best.trees, progress="text")
  }
  
  #Function to back-transform a spatial (raster) prediction
  Back.Transform.Prediction <- function(i) {
    library(raster)
    raster_prediction <- raster_predictions[[i]]
    if (transformation == "log") {
      exp(raster_prediction)-1
    } else if (transformation == "square.root") {
      raster_prediction^2
    } else if (transformation == "fourth.root") {
      raster_prediction^4
    }
  }
  
  #Function to back-transform data
  Back.Transform.Data <- function(d, transformation) {
    if (transformation == "log") {
      exp(d)-1
    } else if (transformation == "square.root") {
      d^2
    } else if (transformation == "fourth.root") {
      d^4
    }
  }
  
  # Import data and setup for analysis --------------------------------------
  # Import data frame with coordinates, response variables, and extracted environmental predictor values (all variables)
  data <- read.csv(file.path(home,"Input_CSVs", island, paste(island,"_Data.csv",sep="")))
  
  # Create list of potential environmental predictors (after highly correlated predictors have been dropped) 
  predictor_set <- as.vector(read.csv(file.path("Input_CSVs", island, paste(island,"_Final_Predictors.csv",sep="")), header=T)[,3])
  
  # Remove data surveys for which the response variable is NA
  data <- data[which(!is.na(data[,response_var])),]
  
  # Apply transformation to response variable and store as new variable in data frame
  response_var_transformed <- paste(transformation,".",response_var,sep="")
  data[[response_var_transformed]] <- Transform.Data(data[[response_var]], transformation) 
  
  # Subset data to include only the coordinates, the response variable, and the extracted values of the potential environment predictors for the specified island group 
  data <- data[,c("Long", "Lat", response_var, response_var_transformed, predictor_set)]
  
  # Create model calibration and model validation subsets by randomly selecting specified percentages of rows in "data"
  calib_ix <- sample(seq(1,nrow(data)), ((100-training_holdout)/100)*nrow(data))
  calib_data <- data[calib_ix,]
  valid_data <- data[-calib_ix,]
  
  # Write calibration and validation data to .csv files
  # write.csv(calib_data, file.path("Input_CSVs", island, paste(island,"_",output_name,"_Calibration.csv",sep="")))
  # write.csv(valid_data, file.path("Input_CSVs", island, paste(island,"_",output_name,"_Validation.csv",sep="")))
  
  # Model Tuning ------------------------------------------------------------
  
  #Cross-validation optimization of boosted regression tree model tuning parameters
  
  # Create lists of model tuning parameter options
  lr <- c(0.01,0.001,0.005) 	#list of options for learning rate
  bag <- c(0.5,0.75) 		#list of options for bag fraction
  tc <- c(2,3,4,5,10) 		#list of options for tree complexity
  
  # Create a data frame of all possible combinations of tuning parameters
  tuning_parameters <- expand.grid(lr,bag,tc)	
  names(tuning_parameters) <- c("learning.rate","bag.fraction","tree.complexity")
  
  # Create a data frame to store boosted regression tree model statistics
  model_tuning_outputs <- data.frame(mean.total.dev=rep(NA,nrow(tuning_parameters)),mean.resid.dev=rep(NA,nrow(tuning_parameters)),cv.mean.dev=rep(NA,nrow(tuning_parameters)),cv.se.dev=rep(NA,nrow(tuning_parameters)),perc.dev.expl=rep(NA,nrow(tuning_parameters)))
  
  #Loop through (in parallel) all possible model tuning parameter combinations, each time fitting a boosted regression tree model with the optimal number of boosting trees
  
  # Create directory for model calibration outputs; will give warning if folder already exists
  dir.create(file.path("Model_Calibration", island), showWarning=F, recursive=TRUE)
  
  # Create vector to store boosted regression tree models for each combination of model tuning parameters
  gbm_models_step <- vector("list", nrow(tuning_parameters))
  
  # Apply the function "Fit.BRT.Model" to each model tuning parameter combination
  t.start <- Sys.time()
  cl <- makeCluster(nCPUs)
  clusterExport(cl, list("response_var_transformed", "predictor_set", "calib_data", "tuning_parameters"))
  gbm_models_step <- parLapply(cl, seq(1,nrow(tuning_parameters)), Fit.BRT.Model)
  stopCluster(cl)
  (t.start - Sys.time())
  
  # For each boosted regression tree model from model parameter tuning, extract statistics
  for(i in 1:nrow(tuning_parameters)) {
    if(is.null(gbm_models_step[[i]])) next
    model_tuning_outputs[i,1] <- gbm_models_step[[i]]$self.statistics$mean.null				# mean total deviance
    model_tuning_outputs[i,2] <- gbm_models_step[[i]]$self.statistics$mean.resid			# mean residual deviance
    model_tuning_outputs[i,3] <- gbm_models_step[[i]]$cv.statistics$deviance.mean			# cross-validation mean residual deviance
    model_tuning_outputs[i,4] <- gbm_models_step[[i]]$cv.statistics$deviance.se				# cross-validation standard error residual deviance
    # Calculate percent deviance explained
    model_tuning_outputs[i,5] <- ((model_tuning_outputs[i,1] - model_tuning_outputs[i,3])/model_tuning_outputs[i,1])*100		
  }
  
  # Attach model statistics to data frame of tuning parameter options
  model_tuning_outputs <- cbind(tuning_parameters, model_tuning_outputs)
  
  # Write model tuning outputs table to csv file
  write.csv(model_tuning_outputs, file.path("Model_Calibration", island, paste(island,"_",output_name,"_BRT_Model_Tuning_Outputs.csv",sep="")))
  
  # Identify the optimal combination of model tuning parameters by identifying the model with the maximum percent deviance explained
  best <- which.max(model_tuning_outputs$perc.dev.expl)
  
  
  # Model Simplification ----------------------------------------------------
  # Use gbm.simplify to assess potential to drop lowest contributing predictors using k-fold cross validation
  ## gbm.simplify identifies sequence of variables to remove
  
  # Run gbm.simplify function on cross-validated model with optimal model tuning parameters 
  gbm_model_simp <- gbm.simplify(gbm_models_step[[best]], n.folds=10, n.drops="auto", alpha=1)
  
  # Create a plot (pdf) of change in predictive deviance as variables are dropped
  pdf(file.path("Model_Calibration", island, paste(island,"_",output_name,"_BRT_Model_Simplification.pdf",sep="")))
  plot(0:length(gbm_model_simp$deviance.summary$mean), c(0,gbm_model_simp$deviance.summary$mean), type="l", 
       xlab="Number of Predictors Dropped", ylab="Change in Predictive Deviance", 
       xaxt="n", ylim=c(min(gbm_model_simp$deviance.summary$mean - 1.96*gbm_model_simp$deviance.summary$se),max(gbm_model_simp$deviance.summary$mean + 1.96*gbm_model_simp$deviance.summary$se)))
  axis(1, at=seq(0,length(gbm_model_simp$deviance.summary$mean),2))
  lines(0:length(gbm_model_simp$deviance.summary$mean), c(0,gbm_model_simp$deviance.summary$mean - 1.96*gbm_model_simp$deviance.summary$se), lty=2)
  lines(0:length(gbm_model_simp$deviance.summary$mean), c(0,gbm_model_simp$deviance.summary$mean + 1.96*gbm_model_simp$deviance.summary$se), lty=2)
  abline(h=0, col="green", lty=2)
  abline(h=gbm_models_step[[best]]$cv.statistics$deviance.se, col="green", lty=2)
  if (length(which(apply(gbm_model_simp$deviance.matrix, 1, mean) < 0)) > 0) {
    abline(v=max(which(apply(gbm_model_simp$deviance.matrix, 1, mean) < 0)), col="red", lty=2)
  } else {
    abline(v=0, col="red", lty=2)
  }
  dev.off()
  
  # Identify number of predictors to drop (number that can be dropped without reduction in predictive deviance)
  if (length(which(apply(gbm_model_simp$deviance.matrix, 1, mean) < 0)) > 0) {
    num_predictors_dropped <- max(which(apply(gbm_model_simp$deviance.matrix, 1, mean) < 0))
  } else {
    num_predictors_dropped <- 0
  }
  
  # Create list of indices of predictors remaining after model simplification
  if (num_predictors_dropped > 0) {
    final_vars_ix <- gbm_model_simp$pred.list[[num_predictors_dropped]]
  } else {
    final_vars_ix <- which(names(calib_data) %in% predictor_set)
  }
  
  # Create list of the names of the predictors remaining after model simplification
  final_vars <- names(calib_data[,final_vars_ix])
  
  
  # Fit Final Model ---------------------------------------------------------
  gbm_model_final <- gbm.step(data=calib_data,
                              gbm.x = final_vars_ix,
                              gbm.y = grep(response_var_transformed,names(calib_data)),
                              tree.complexity=gbm_model_simp$gbm.call$tree.complexity,
                              learning.rate=gbm_model_simp$gbm.call$learning.rate,
                              bag.fraction=gbm_model_simp$gbm.call$bag.fraction,
                              family="gaussian",
                              n.trees=gbm_model_simp$gbm.call$best.trees,
                              max.trees=gbm_model_simp$gbm.call$best.trees)
  
  
  
  # Bootstrap final model ---------------------------------------------------
  #Using bootstrapping, fit a set of final models using the optimal model tuning parameters and simplified predictor set
  
  # Create a set of bootstrap samples
  bootstrap_samples <- vector("list", num_bootstraps)
  bootstrap_samples <- lapply(1:num_bootstraps, function(i) {sample(1:nrow(calib_data), size=nrow(calib_data), replace=TRUE)})
  
  # Create vector to store final models
  gbm_models_final <- vector("list", num_bootstraps)	
  
  # Using parallelization, generate a model for each bootstrap sample using the optimal model tuning parameters and simplified predictor set
  t.start <- Sys.time()
  cl <- makeCluster(nCPUs)
  clusterExport(cl, list("calib_data","bootstrap_samples","gbm_model_simp", "response_var_transformed", "final_vars_ix"))
  gbm_models_final <- parLapply(cl, seq(1,num_bootstraps), Fit.BRT.Model.Fixed.Bootstrap)
  stopCluster(cl)
  (t.start - Sys.time())
  
  
  # Spatial Prediction ------------------------------------------------------
  
  # Create directory for spatial prediction outputs; will give warning if folder already exists
  dir.create(file.path("Spatial_Predictions", island), showWarning=F, recursive=TRUE)
  
  # Create raster stack of final set of environmental predictors
  
  # Create a list of all the predictor geotiff files for the island
  raster_list <- list.files(file.path("Predictor_geotiffs_new", island), pattern="tif$")
  
  # Subset raster list to include only the predictors in the simplified predictor set
  raster_list <- raster_list[match(paste(final_vars,".tif",sep=""), raster_list)]
  
  # Create raster stack of RasterLayer objects from geotiffs of the simplified predictor set
  raster_stack <- stack(file.path("Predictor_geotiffs_new", island, raster_list))
  
  # Create vector to store raster predictions
  raster_predictions <- vector("list", num_bootstraps)
  
  # Using parallelization, generate a spatial prediction for each of the final models
  t.start <- Sys.time()
  cl <- makeCluster(nCPUs)
  clusterExport(cl, list("gbm_models_final", "raster_stack"))
  raster_predictions <- parLapply(cl, seq(1,num_bootstraps), Spatial.Prediction)
  stopCluster(cl)
  (t.start - Sys.time())
  
  # Back transform raster predictions
  
  # Create vector to store backtransformed raster predictions
  raster_predictions_backtransformed <- vector("list", num_bootstraps)
  
  # Using parallelization, back transform each of the spatial predictions
  cl <- makeCluster(nCPUs)
  clusterExport(cl, list("raster_predictions","transformation"))
  raster_predictions_backtransformed <- parLapply(cl, seq(1,num_bootstraps), Back.Transform.Prediction)
  stopCluster(cl)
  
  # Calculate mean and standard deviation of back-transformed raster predictions; export to geotiff
  mean_raster_prediction <- stackApply(stack(raster_predictions_backtransformed), indices=rep(1,num_bootstraps), fun=mean, na.rm=FALSE)
  proj4string(mean_raster_prediction) <- "+proj=omerc +lat_0=20.5713 +lonc=-157.1895 +alpha=-60 +k=0.9996 +x_0=0 +y_0=0 +gamma=0 +datum=WGS84 +units=m +no_uoff +no_defs +ellps=WGS84 +towgs84=0,0,0"
  writeRaster(mean_raster_prediction, filename=file.path("Spatial_Predictions", island, paste(island,"_Mean_Predicted_",output_name,sep="")), format="GTiff", overwrite=TRUE)
  sd_raster_prediction <- stackApply(stack(raster_predictions_backtransformed), indices=rep(1,num_bootstraps), fun=sd, na.rm=FALSE)
  proj4string(sd_raster_prediction) <- "+proj=omerc +lat_0=20.5713 +lonc=-157.1895 +alpha=-60 +k=0.9996 +x_0=0 +y_0=0 +gamma=0 +datum=WGS84 +units=m +no_uoff +no_defs +ellps=WGS84 +towgs84=0,0,0"
  writeRaster(sd_raster_prediction, filename=file.path("Spatial_Predictions", island, paste(island,"_StdDev_Predicted_",output_name,sep="")), format="GTiff", overwrite=TRUE)
  
  
  # Save R workspace --------------------------------------------------------
  dir.create(file.path("R_Workspaces", island), showWarning=F, recursive=TRUE)
  save(list=ls(all.names=TRUE), file=file.path("R_Workspaces", island, paste(island,"_",output_name,"_BRT_Model_Workspace.RData", sep="")))
  
  # Check final model residuals for spatial autocorrelation -----------------
  
  # Create data frame to store the coordinates of the calibration data surveys and the model residual values         vector to store p values for each calculation of Moran's I
  model_residuals <- data.frame(calib_data$Long, calib_data$Lat, gbm_model_final$residuals)
  names(model_residuals) <- c("Long","Lat","Resid")
  
  # Calculate a distance matrix containing the distances between each pair of survey data locations	
  pt_dists <- as.matrix(dist(cbind(model_residuals$Long, model_residuals$Lat)))
  
  # Convert distance matrix to an inverse distance matrix	
  pt_dists_inv <- 1 / pt_dists
  
  # For coincident survey data locations, convert inverse distance from Inf to Zero (otherwise, calculation of Moran's I will fail)
  pt_dists_inv[is.infinite(pt_dists_inv)] <- 0
  
  # Calculate Moran's I autocorrelation coefficent of the model residuals using the inverse distance matrix of weights; report p-value
  ## if p-value suggests that observed value of I is significantly different from the expected value, then the residuals are autocorrelated
  Morans_I_test <- Moran.I(model_residuals$Resid, pt_dists_inv)$p.value
  
  
  # Evaluate Predictor Variable Importance ----------------------------------
  
  # Create directory for model summary outputs; will give warning if folder already exists
  dir.create(file.path("Model_Summaries", island), showWarning=F, recursive=TRUE)
  
  # Calculate mean relative importance from set of final models      REQUIRES GBM PACKAGE
  relative_influence_df <- data.frame()
  for (i in seq(1,num_bootstraps)) {
    relative_influence <- summary(gbm_models_final[[i]], plotit=FALSE)
    relative_influence_df <- rbind(relative_influence_df, relative_influence[order(row.names(relative_influence)),]$rel.inf)
  }
  rm(relative_influence)
  names(relative_influence_df) <- final_vars
  mean_relative_influence <- colMeans(relative_influence_df)
  relative_influence_df <- relative_influence_df[,order(mean_relative_influence, decreasing=TRUE)]
  
  # Create boxplot of distribtions of relative influence by environmental predictor
  pdf(file.path("Model_Summaries", island, paste(island,"_",output_name,"_Predictor_Importance.pdf",sep="")))
  boxplot(relative_influence_df, xlab="", ylab="Relative Contribution (%)", ylim=c(0,100), las=2)
  dev.off()
  
  # Create barplot of mean relative influence by environmental predictor
  pdf(file.path("Model_Summaries", island, paste(island,"_",output_name,"_Mean_Predictor_Importance.pdf",sep="")))
  barplot(mean_relative_influence[order(mean_relative_influence, decreasing=TRUE)], xlab="", ylab="Mean Relative Contribution (%)", ylim=c(0,100), las=2)
  dev.off()
  
  
  # Create Partial Dependence Plots -----------------------------------------
  
  # Create directory for partial dependence plots; will give warning if folder already exists
  dir.create(file.path("Partial_Dependence_Plots", island), showWarning=F, recursive=TRUE)
  partial_dependence <- vector("list", length(final_vars_ix))
  for (i in seq(1,length(final_vars_ix))) {
    for (j in seq(1,num_bootstraps)) {
      if (j == 1) {
        partial_dependence[[i]] <- plot(gbm_models_final[[j]], i.var = i, return.grid=TRUE)
      } else {
        partial_dependence[[i]] <- cbind(partial_dependence[[i]], plot(gbm_models_final[[j]], i.var = i, return.grid=TRUE)[2])
      }
    }
    names(partial_dependence[[i]]) <- c(final_vars[i],paste("y",seq(1,num_bootstraps),sep=""))
    partial_dependence[[i]]$mean <- apply(partial_dependence[[i]][,c(paste("y",seq(1,num_bootstraps),sep=""))], 1, mean)
    partial_dependence[[i]]$sd <- apply(partial_dependence[[i]][,c(paste("y",seq(1,num_bootstraps),sep=""))], 1, sd)
    partial_dependence[[i]]$quant5 <- apply(partial_dependence[[i]][,c(paste("y",seq(1,num_bootstraps),sep=""))], 1, quantile, probs=c(0.05))
    partial_dependence[[i]]$quant95 <- apply(partial_dependence[[i]][,c(paste("y",seq(1,num_bootstraps),sep=""))], 1, quantile, probs=c(0.95))
    pdf(file.path("Partial_Dependence_Plots", island, paste(island,"_",output_name,"_Partial_Dependence_Plot_",final_vars[i],".pdf",sep="")))
    plot(partial_dependence[[i]][,1],partial_dependence[[i]]$mean,type="n",xlab="",ylab="")
    #polygon(c(partial_dependence[[i]][,1],rev(partial_dependence[[i]][,1])), c(partial_dependence[[i]]$quant5,rev(partial_dependence[[i]]$quant95)), col="grey80", border = NA)
    polygon(c(partial_dependence[[i]][,1],rev(partial_dependence[[i]][,1])), c((partial_dependence[[i]]$mean-partial_dependence[[i]]$sd),rev((partial_dependence[[i]]$mean+partial_dependence[[i]]$sd))), col="grey80", border = NA)
    lines(partial_dependence[[i]][,1],partial_dependence[[i]]$mean)
    dev.off()
  }
  
  
  # Evaluate Pairwise Interactions among the Predictors ---------------------
  
  #HOLD FOR LATER
  # Create directory for perspective (interaction) plots; will give warning if folder already exists
  #dir.create(file.path("Perspective_Plots", island), showWarning=TRUE, recursive=TRUE)
  
  
  # Model Validation --------------------------------------------------------
  
  # Calculate the percent deviance explained by the final model (full set of calibration data rather than from cross validation)
  ## (not to be used to evaluate model, just to demonstrate overfitting) ##
  percent_deviance_explained_train <- ((gbm_model_final$self.statistics$mean.null-gbm_model_final$self.statistics$mean.resid)/gbm_model_final$self.statistics$mean.null)*100
  
  # Calculate the mean and standard error percent deviance explained by the final cross validated model
  mean_percent_deviance_explained_cv <- ((gbm_model_final$self.statistics$mean.null-gbm_model_final$cv.statistics$deviance.mean)/gbm_model_final$self.statistics$mean.null)*100
  
  percent_deviance_explained_cv_upper <- ((gbm_model_final$self.statistics$mean.null-(gbm_model_final$cv.statistics$deviance.mean-gbm_model_final$cv.statistics$deviance.se))/gbm_model_final$self.statistics$mean.null)*100
  
  percent_deviance_explained_cv_lower <- ((gbm_model_final$self.statistics$mean.null-(gbm_model_final$cv.statistics$deviance.mean+gbm_model_final$cv.statistics$deviance.se))/gbm_model_final$self.statistics$mean.null)*100
  
  se_percent_deviance_explained_cv <- (percent_deviance_explained_cv_upper-percent_deviance_explained_cv_lower)/2
  
  # Calculate percent deviance explained by the validation data
  # Extract observed values from the validation data
  obsValues_test <- valid_data[,response_var_transformed]
  
  # Extract predicted (fitted) values at the validation data locations
  predValues_test <- predict(gbm_model_final, valid_data, n.trees=gbm_model_final$gbm.call$best.trees, type="response")
  
  # Calculate percent deviance explained
  total_deviance_test <- sum((obsValues_test-mean(obsValues_test))*(obsValues_test-mean(obsValues_test)))
  residual_deviance_test <- calc.deviance(obsValues_test, predValues_test, family="gaussian", calc.mean=FALSE)
  percent_deviance_explained_test <- ((total_deviance_test-residual_deviance_test)/total_deviance_test)*100
  
  #Calculate RMSE, MAE, and MAPE for the final model using the validation data		  ###In TRANSFORMED UNITS?###
  # Extracted observed values of untransformed response variable from the validation data
  # Back transform vectors of observed values and predicted values
  obsValues_test_untransformed <- valid_data[,response_var]
  predValues_test_backtransformed <- Back.Transform.Data(predValues_test, transformation)
  
  bias <- mean(obsValues_test_untransformed-predValues_test_backtransformed, na.rm=TRUE)
  MAE <- mean(abs(obsValues_test_untransformed-predValues_test_backtransformed), na.rm=TRUE)
  RMSE <- sqrt(mean((obsValues_test_untransformed-predValues_test_backtransformed)^2, na.rm=TRUE))
  
  
  # Summarize Model Outputs -------------------------------------------------
  
  model_summary <- c(island, output_name, response_var, transformation, 
                     gbm_model_final$gbm.call$tree.complexity, gbm_model_final$gbm.call$learning.rate, 
                     gbm_model_final$gbm.call$bag.fraction, gbm_model_final$gbm.call$best.trees,
                     percent_deviance_explained_train, mean_percent_deviance_explained_cv, se_percent_deviance_explained_cv,
                     percent_deviance_explained_test, bias, MAE, RMSE)
  model_summary_names <- c("Island", "Model Name", "Response Variable", "Data Transformation", "Tree Complexity", "Learning Rate", "Bag Fraction", "Number of Trees","Training Percent Deviance Explained", "Cross-Validation Mean Percent Deviance Explained", "Cross-Validation SE Percent Deviance Explained","Test Percent Deviance Explained", "Bias", "MAE", "RMSE")
  model_summary_df <- data.frame(model.summary.names=model_summary_names, model_summary=model_summary)
  write.table(model_summary_df, file.path("Model_Summaries", island, paste(island,"_",output_name,"_BRT_Model_Summary.csv",sep="")), row.names=FALSE, col.names=FALSE, sep=",")
  
  
  # Save R workspace --------------------------------------------------------
  save(list=ls(all.names=TRUE), file=file.path("R_Workspaces", island, paste(island,"_",output_name,"_BRT_Model_Workspace.RData", sep="")))
  
  (tot.time <- tot.start-Sys.time())
  rm('raster_predictions_backtransformed')
  rm('raster_predictions')
  rm('gbm_models_step')
  rm('gbm_models_final')
  rm('mean_raster_prediction')
  
}