
#Import functions
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

#	require(ape)
	require(betareg)
#	require(safedata)
#	set_safe_dir("C:/Users/rewers/OneDrive - Imperial College London/work/SAFE project - private/safedata",update = FALSE)

#Import data
	thresh.data <- readRDS("data/threshold_taxa_data.rds")	#Site x species matrix for each dataset, containing only taxa that are duplicated across >1 dataset
	lidar.data <- readRDS("data/lidar_percent.rds")			#Lidar data for all sites in full dataset
	func.groups <- readRDS('data/functional_groups.rds')	#Functional groups
#	taxa <- readRDS('data/taxon_table.rds')					#Full list of all taxa
#	map <- read.table('data/species_families_order_map.txt', sep = '-')	#Identified one family and example species per order that exists on TimeTree.org
#	tr <- read.tree("data/species.nwk")						#Imported phylogeny from TimeTree (www.timetree.org)
	bayes_results <- readRDS("data/bayes.rds")[[1]]			#Results from Bayesian slopes analysis (Replicability analysis)
	
#Fit models and calculate summaries
#	fitted_thresh <- fit.models(full_data = thresh.data, lidar = lidar.data, min.observs = 5,
#		predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000'))
#		saveRDS(fitted_thresh, 'results/fitted_thresh.rds')
	fitted_func <- fit.models(full_data = thresh.data, func_data = func.groups, lidar = lidar.data,
		predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000'), min.observs = 5)
	saveRDS(fitted_func, 'results/fitted_func.rds')
#	turn_points <- turns(fitted_mod = fitted_thresh)
#		saveRDS(turn_points, 'results/turn_points.rds')
#	func_points <- turns(fitted_mod = fitted_func)
#		saveRDS(func_points, 'results/func_points.rds')
#	break_points <- func.thresholds(func_groups = func.groups, turns_taxa = turn_points)
#		saveRDS(break_points, 'results/break_points.rds')

#Read in pre-calculated versions
	fitted_thresh <- readRDS('results/fitted_thresh.rds')
	fitted_func <- readRDS('results/fitted_func.rds')
		
##DELETE WHEN RE-RUN WITHOUT GHOST TAXA
		fitted_thresh <- fitted_thresh[!is.na(fitted_thresh$num.occs), ]	#Remove ghost taxa 
		fitted_func <- fitted_func[!is.na(fitted_func$num.occs), ]	#Remove ghost taxa 
##############
	turn_points <- readRDS('results/turn_points.rds')
	func_points <- readRDS('results/func_points.rds')
	break_points <- readRDS('results/break_points.rds')
	func_points <- readRDS('results/func_points.rds')

	taxa_cats <- assign.taxon(dataset = fitted_thresh[fitted_thresh$num.occs >= 5, ],
		taxon_table = taxa)
	turn_points <- assign.taxon(dataset = turn_points, taxon_table = taxa)
	phylo <- arrange.phylo(timetree = tr, raw_data = thresh.data, taxa_safe = taxa,
		tt_map = map, coefs = fitted_thresh, palette_col = NA)
	taxaXdata <- taxon.dataset(thresh.data)	#Create taxon x survey matrix



#Abstract
	#Number of taxa modelled
		nrow(fitted_thresh)
	#Number of orders containing taxa that were analysed
		sum(!is.na(phylo$numbers3$propTax))
	#Number of functional groups
		nrow(fitted_func)
	#Proactive conservation
		#Value obtained from figures
	#Reactive conservation
		#Value obtained from figures

#Summary data
	#Number of surveys
		length(thresh.data)
	#Number of taxa
		nrow(fitted_thresh)		#all taxa
	#Number of functional groups
		nrow(fitted_func)	
	#Number of modelled taxa
		sum(as.numeric(fitted_thresh$num.occs) >= 5)	#taxa with >= minimum number of occurrences
	#Higher order taxa
		length(unique(taxa_cats$matched.taxa$order))		#Number of orders
		length(unique(taxa_cats$matched.taxa$genus))		#Number of genera
		summary(factor(taxa_cats$dataset$TaxonType))		#Taxonomic categories
		length(which(taxa_cats$matched.taxa$order == 'Coleoptera'))		#Number of beetles
		length(which(taxa_cats$matched.taxa$order == 'Lepidoptera'))	#Number of leps
		length(which(taxa_cats$matched.taxa$family == 'Formicidae'))	#Number of ants
	#Functional groups
	minsize <- print(func.groups[which(func.groups$BodyMass == min(func.groups$BodyMass, na.rm = TRUE)), c('taxon_name', 'BodyMass')])
	maxsize <- print(func.groups[which(func.groups$BodyMass == max(func.groups$BodyMass, na.rm = TRUE)), c('taxon_name', 'BodyMass')])
	log10(maxsize$BodyMass[1]) - log10(minsize$BodyMass[1])		#Orders of magnitude in body size
	func.summary(func_groups = func.groups)



	#Number of taxa instantly impacted
		sum(turn_points$dataset$turn.point == 0, na.rm = TRUE)
		sum(turn_points$dataset$turn.point == 0, na.rm = TRUE) / nrow(turn_points$dataset)
	#Number of functional groups instantly impacted
		sum(func_points$turn.point == 0, na.rm = TRUE)
		sum(func_points$turn.point == 0, na.rm = TRUE) / nrow(func_points)
	#Number negative vs positive taxon responses
		sum(turn_points$dataset$slope < 0 & turn_points$dataset$pval < 0.05)	#number responding negatively
		sum(turn_points$dataset$slope > 0 & turn_points$dataset$pval < 0.05)	#number responding positively
	#Number negative vs positive functional group responses
		sum(func_points$slope < 0 & func_points$pval < 0.05)	#number responding negatively
		sum(func_points$slope > 0 & func_points$pval < 0.05)	#number responding positively
	#Proportion of positive responders
		sum(turn_points$dataset$slope > 0 & turn_points$dataset$pval < 0.05) / nrow(turn_points$dataset)	#taxa
		sum(func_points$slope > 0 & func_points$pval < 0.05) / nrow(func_points)	#functional groups
	

#Figure 1 caption
	#Number of taxa
		nrow(fitted_thresh)		#all taxa
	#Number of functional groups
		nrow(fitted_func)	


#Ecological thresholds
	#Thresholds in taxa turning points
		taxa_turn <- turn_points$dataset$turn.point[!is.na(turn_points$dataset$turn.point)]
		break.points(density(taxa_turn)$x, density(taxa_turn)$y)
	#Thresholds in functional group turning points
		func_turn <- func_points$turn.point[!is.na(func_points$turn.point)]
		break.points(density(func_turn)$x, density(func_turn)$y)
	#Thresholds in taxa peak rate
		taxa_rate <- turn_points$dataset$maxrate[!is.na(turn_points$dataset$maxrate)]
		break.points(density(taxa_rate)$x, density(taxa_rate)$y)
	#Thresholds in functional group peak rate
		func_rate <- func_points$maxrate[!is.na(func_points$maxrate)]
		break.points(density(func_rate)$x, density(func_rate)$y)




##proactive threshold
	#for each of the four break.points, get first 'accelearate' value and take mean + error
#reactive threshold
	#for each of the four break.points, take last 'accelerate' value and take mean + error
#Then apply rounding...
#Use bootstrapping of input data to generate distribution of break.points estimates for getting 95% CI?













	


#Taxonomic categories
	#Number of orders
		sum(!is.na(phylo$numbers3$propTax))
	#Number of orders with impacted taxa
		func <- function(x) sum(!is.na(x))/length(x)
		reps <- by(turn_points$dataset$turn.point, factor(turn_points$dataset$Order), FUN = func)
		sum(reps > 0)				#Number of orders containing impacted taxa
		sum(reps > 0)/length(reps)	#As a proportion
	#Proportion taxa with turning points
		by(turn_points$dataset$turn.point, factor(turn_points$dataset$TaxonType), FUN = func)	#Proportion taxa with turning points
		by(turn_points$dataset$turn.point, factor(turn_points$dataset$TaxonType), FUN = mean, na.rm = TRUE) 	#Mean turning point


#Connection with replicability
	resil <- print(resil.dat(turns = turn_points$dataset, bayes = bayes_results, 
		grouping = 'TaxonType', full_taxa = taxa))
	summary(betareg(resilience ~ bayes_prob, data = resil))


plot(resil$bayes_prob, resil$resilience)






resil_func <- resil.func(func_groups = func.groups, func_points = func_points, turns = turn_points$dataset)







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
	

#Methods
	#Taxon x survey combinations
		sum(!is.na(taxaXdata))
	#Number of data sources
		data_sum <- summarise.data(thresh.data)
		nrow(data_sum)
	#Number of surveys
		sum(data_sum$surveys)
	#Total number of taxa in all surveys combined
		nrow(taxaXdata)
	#Number of taxa identified to species
		used_taxa <- taxa[match(rownames(taxaXdata), taxa$taxon_name), ]
		sum(used_taxa$taxon_level == 'species', na.rm = TRUE)
	#Number identified to morphospecies
		sum(used_taxa$taxon_level == 'morphospecies', na.rm = TRUE)
	#Number of taxa modelled
		nrow(fitted_thresh)
	#Number taxa represented in >1 survey
		modelled <- taxaXdata[match(fitted_thresh$taxon, rownames(taxaXdata)),]	#Subset to just modelled taxa
		sum(apply(X = modelled, MARGIN = 1, FUN = function(x) sum(!is.na(x))) > 1)	#Number of taxa in multiple surveys


#Fig S1 caption
	#Number of orders in all surveys combined
		length(phylo$tr$tip.label)
	#Number of orders containing taxa that were analysed
		sum(!is.na(phylo$numbers3$propTax))


#Methods
	#Number taxa analysed with GLMM
		summary(factor(fitted_thresh$modtype))
		summary(factor(fitted_thresh$modtype)) / (nrow(fitted_thresh) - sum(is.na(fitted_thresh$modtype)))
	#Number functional groups analysed with GLMM
		summary(factor(fitted_func$modtype))
		summary(factor(fitted_func$modtype)) / (nrow(fitted_func) - sum(is.na(fitted_func$modtype)))
	
	

#Table S1
	data_sum
	
