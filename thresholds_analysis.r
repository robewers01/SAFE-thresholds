
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
	fitted_mod <- fit.models(taxa_data = thresh.data, lidar = lidar.data, min.obs = 5,
		predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000'))
		
		



