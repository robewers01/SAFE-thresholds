
#Import functions
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

	require(safedata)
	set_safe_dir("C:/Users/rewers/OneDrive - Imperial College London/work/SAFE project - private/safedata",update = FALSE)

#Import data
	thresh.data <- readRDS("data/threshold_taxa_data.rds")	#Site x species matrix for each dataset, containing only taxa that are duplicated across >1 dataset
	lidar.data <- readRDS("data/lidar_percent.rds")			#Lidar data for all sites in full dataset
	
#Fit models and calculate summaries
#	fitted_thresh <- fit.models(taxa_data = thresh.data, lidar = lidar.data, min.obs = 5,
#		predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000'))
#		saveRDS(fitted_mod, 'results/fitted_thresh.rds')
#	turn_points <- turns(fitted_thresh)
#		saveRDS(turn_points, 'results/turn_points.rds')

#Read in pre-calculated versions
	fitted_thresh <- readRDS('results/fitted_thresh.rds')
	turn_points <- readRDS('results/turn_points.rds')

	taxa_summary <- summarise.taxa(full_data = thresh.data, min.obs = 5)
	taxa_cats <- assign.taxon(dataset = data.frame(taxon = taxa_summary$modelled.taxa))

#Summary data
	#Number of surveys
		length(thresh.data)
	#Number of taxa
		taxa_summary$num.all.taxa		#all taxa
		taxa_summary$num.modelled.taxa	#taxa with >= minimum number of occurrences
	#Higher order taxa
		length(unique(taxa_cats$matched.taxa$order))		#Number of orders
		length(unique(taxa_cats$matched.taxa$genus))		#Number of genera
		summary(factor(taxa_cats$dataset$TaxonType))		#Taxonomic categories
		length(which(taxa_cats$matched.taxa$order == 'Coleoptera'))		#Number of beetles
		length(which(taxa_cats$matched.taxa$order == 'Lepidoptera'))	#Number of leps
		length(which(taxa_cats$matched.taxa$family == 'Formicidae'))	#Number of ants
		

#Turnpoints
	#Number of taxa instantly impacted
		sum(turn_points$turn.point == 0, na.rm = TRUE)
		sum(turn_points$turn.point == 0, na.rm = TRUE) / nrow(turn_points)
	#Number of taxa with significant turnpoints
		sum(!is.na(turn_points$turn.point))
		sum(!is.na(turn_points$turn.point)) / nrow(turn_points)

	



