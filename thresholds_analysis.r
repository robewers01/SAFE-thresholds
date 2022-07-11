
#Import functions
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

	require(safedata)
	set_safe_dir("C:/Users/rewers/OneDrive - Imperial College London/work/SAFE project - private/safedata",update = FALSE)

#Import data
	thresh.data <- readRDS("data/threshold_taxa_data.rds")	#Site x species matrix for each dataset, containing only taxa that are duplicated across >1 dataset
	lidar.data <- readRDS("data/lidar_percent.rds")		#Lidar data for all sites in full dataset

#Fit models and calculate summaries
	fitted_thresh <- fit.models(taxa_data = thresh.data, lidar = lidar.data, min.obs = 5,
		predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000'))
		saveRDS(fitted_mod, 'results/fitted_thresh.rds')
		


#Function to estimate turning points for fitted models

fitted_mod = fitted_thresh

function(fitted_mod){

	turnpoints <- matrix(NA, nrow = length(fitted_mod$models), ncol = 2)
	for(i in 1:length(fitted_mod$models)){			#For each fitted model....
		target = fitted_mod$models[[i]]
		#Check if fitted model was significant or not
		if(as.numeric(fitted_mod$coef$pval[i]) > 0.05){
			turnpoints[i, ] <- rep(NA, 2)		#No informative turning points
			}else{
				#Find turning points
				fits <- fitted.vals(target)		#Estimate fitted values and derivatives
				turns <- root.finder(y = d2.vals, x = predx, method = "pastecs")
				
				
				}
		
		
		
		}
	
	
	
	
	}


####
###

	
	
#Function to find roots and turning points
root.finder <- function(fitted_vals){
	#fitted_vals = output from fitted.vals
#	#y = vector of y-axis values to scan
#	#x = vector of x-axis values at which y-axis values are evaluated
	
	if (!require(pastecs)) install.packages("pastecs") && require(pastecs)   ## Check if required packages are installed
	 
	x <- fitted_vals$agb
	y <- round(fitted_vals$d2,8)
	
	past <- NA
	try(past <- turnpoints(y)$tp, silent = TRUE)
	try(peak <- turnpoints(y)$firstispeak, silent = TRUE)
	if(!is.na(past[1])){
		first <- x[past]
		}else{
			first <- numeric()
			peak <- NA
			}
		
	return(list(tps = min(first), peak = peak))		
	}
