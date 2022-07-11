
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

turns <- function(fitted_mod){
	#fitted_mod = output from fit.models

	turnpoints <- matrix(NA, nrow = length(fitted_mod$models), ncol = 2)
	for(i in 1:length(fitted_mod$models)){			#For each fitted model....
		target = fitted_mod$models[[i]]
		#Check if fitted model was significant or not
		if(as.numeric(fitted_mod$coef$pval[i]) < 0.05){	#If it was a significant model
			#Find turning points
			fits <- fitted.vals(target, mod_coef = fitted_mod$coefs[i , ])		#Estimate fitted values and derivatives
			turns <- root.finder(fitted_vals = fits)
			turnsdat <- unlist(extract.turns(turns))
			turnpoints[i ,] <- turnsdat
			}
		}
	}


####
###

	

