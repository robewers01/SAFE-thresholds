
#Import functions
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

#	require(safedata)
#	set_safe_dir("C:/Users/rewers/OneDrive - Imperial College London/work/SAFE project - private/safedata",update = FALSE)

#Import data
	thresh.data <- readRDS("data/threshold_taxa_data.rds")	#Site x species matrix for each dataset, containing only taxa that are duplicated across >1 dataset
	lidar.data <- readRDS("data/lidar_percent.rds")			#Lidar data for all sites in full dataset
	func.groups <- readRDS('data/functional_groups.rds')			#Functional groups
	
#Fit models and calculate summaries
#	fitted_thresh <- fit.models(full_data = thresh.data, lidar = lidar.data, min.observs = 5,
#		predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000'))
#		saveRDS(fitted_thresh, 'results/fitted_thresh.rds')
#	fitted_func <- fit.models(full_data = thresh.data, func_data = func.groups, lidar = lidar.data,
#		predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000'), min.observs = 5)
#	saveRDS(fitted_func, 'results/fitted_func.rds')
#	turn_points <- turns(fitted_mod = fitted_thresh)
#		saveRDS(turn_points, 'results/turn_points.rds')
#	func_points <- turns(fitted_mod = fitted_func)
#		saveRDS(func_points, 'results/func_points.rds')
#	break_points <- func.thresholds(func_groups = func.groups, turns_taxa = turn_points)
#		saveRDS(break_points, 'results/break_points.rds')

#Read in pre-calculated versions
	fitted_thresh <- readRDS('results/fitted_thresh.rds')
		fitted_thresh <- fitted_thresh[!is.na(fitted_thresh$num.occs), ]	#Remove taxa that weren't found for analyses
	fitted_func <- readRDS('results/fitted_func.rds')
		fitted_func <- fitted_func[!is.na(fitted_func$num.occs), ]	#Remove taxa that weren't found for analyses
	turn_points <- readRDS('results/turn_points.rds')
	func_points <- readRDS('results/func_points.rds')
	break_points <- readRDS('results/break_points.rds')
	func_points <- readRDS('results/func_points.rds')

	taxa_cats <- assign.taxon(dataset = fitted_thresh[fitted_thresh$num.occs >= 5, ])

#Summary data
	#Number of surveys
		length(thresh.data)
	#Number of taxa
		nrow(fitted_thresh)		#all taxa
		sum(as.numeric(fitted_thresh$num.occs) >= 5)	#taxa with >= minimum number of occurrences
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
	#Number of functional groups instantly impacted
		sum(func_points$turn.point == 0, na.rm = TRUE)
		sum(func_points$turn.point == 0, na.rm = TRUE) / nrow(func_points)
	#Number negative taxon responses
		sum(turn_points$slope[turn_points$pval < 0.05] < 0)		#are taxa responding negatively?
		sum(turn_points$slope[turn_points$pval < 0.05] < 0)	 / sum(as.numeric(turn_points$pval) < 0.05)
	#Number negative functional group responses
		sum(func_points$slope[func_points$pval < 0.05] < 0)		#are taxa responding negatively?
		sum(func_points$slope[func_points$pval < 0.05] < 0)	 / sum(as.numeric(func_points$pval) < 0.05)


#Ecological thresholds
	#Thresholds in taxa turning points
		taxa_turn <- turn_points$turn.point[!is.na(turn_points$turn.point)]
		break.points(density(taxa_turn)$x, density(taxa_turn)$y)
	#Thresholds in functional group turning points
		func_turn <- func_points$turn.point[!is.na(func_points$turn.point)]
		break.points(density(func_turn)$x, density(func_turn)$y)
	#Thresholds in taxa peak rate
		taxa_rate <- turn_points$maxrate[!is.na(turn_points$maxrate)]
		break.points(density(taxa_rate)$x, density(taxa_rate)$y)
	#Thresholds in functional group peak rate
		func_rate <- func_points$maxrate[!is.na(func_points$maxrate)]
		break.points(density(func_rate)$x, density(func_rate)$y)




#Taxonomic analysis
	#Number of taxa with significant turnpoints
		sum(!is.na(turn_points$turn.point))
		sum(!is.na(turn_points$turn.point)) / nrow(turn_points)



#Functional composition - example taxa
	#Habitat strata generalists
	find.egs(Function = 'StrataGeneralism_Categorised', function_qual = 'high', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Trophic
	find.egs(Function = 'TrophicGeneralism_Categorised', function_qual = 'high', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Endangered
	find.egs(Function = 'IUCNthreat', function_qual = 'Threatened', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Plants
	find.egs(Function = 'PlantPhotoCap_Categorised', function_qual = 'high', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Social animals
	find.egs(Function = 'Sociality', function_qual = 'social', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	find.egs(Function = 'Sociality', function_qual = 'pair', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Carnivores
	find.egs(Function = 'DietVert', function_qual = 'vertivore', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Subterranean mammals
	find.egs(Function = 'StraSub', function_qual = 'subterranean', taxtype = 'mammal', 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	






