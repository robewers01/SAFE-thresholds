
#Import functions
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

	require(ape)
	require(betareg)
	require(safedata)
	set_safe_dir("C:/Users/rewers/OneDrive - Imperial College London/work/SAFE project - private/safedata",update = FALSE)

#Import data
	thresh.data <- readRDS("data/threshold_taxa_data.rds")	#Site x species matrix for each dataset, containing only taxa that are duplicated across >1 dataset
	lidar.data <- readRDS("data/lidar_percent.rds")			#Lidar data for all sites in full dataset
	func.groups <- readRDS('data/functional_groups.rds')	#Functional groups
	bayes_results <- readRDS("data/bayes.rds")[[1]]			#Results from Bayesian slopes analysis (Replicability analysis)
	tr <- read.tree("data/species.nwk")			#Imported phylogeny from TimeTree (www.timetree.org)
	taxa <- readRDS('data/taxon_table.rds')					#Full list of all taxa
	map <- read.table('data/species_families_order_map.txt', sep = '-')	#Identified one family and example species per order that exists on TimeTree.org
	
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

#Read in pre-calculated versions
	fitted_thresh <- readRDS('results/fitted_thresh.rds')
	fitted_func <- readRDS('results/fitted_func.rds')
	turn_points <- readRDS('results/turn_points.rds')
	func_points <- readRDS('results/func_points.rds')

#Calculate summaries and derived data
	taxa_cats <- assign.taxon(dataset = fitted_thresh[!is.na(fitted_thresh$modtype), ], taxon_table = taxa)
	turn_points <- assign.taxon(dataset = turn_points, taxon_table = taxa)
	phylo <- arrange.phylo(timetree = tr, raw_data = thresh.data, taxa_safe = taxa,
		tt_map = map, coefs = fitted_thresh, palette_col = NA)
	taxaXdata <- taxon.dataset(thresh.data)


#Abstract
	#Number of taxa modelled
		nrow(turn_points$dataset)
	#Number of orders containing taxa that were analysed
		length(unique(taxa_cats$matched.taxa$order))
	#Number of functional groups
		nrow(fitted_func)
	#Conservation thresholds
		estimate.thresholds(turn_points, func_points)
	
#Summary data
	#Number of surveys
		length(thresh.data)
	#Number of taxa
		nrow(fitted_thresh)		#all taxa
	#Number of functional groups
		nrow(fitted_func)	
	#Number of modelled taxa
		sum(!is.na(fitted_thresh$modtype))
	#Higher order taxa
		length(unique(taxa_cats$matched.taxa$order))		#Number of orders
		length(unique(taxa_cats$matched.taxa$genus))		#Number of genera
		summary(factor(taxa_cats$dataset$TaxonType))		#Taxonomic categories
		length(which(taxa_cats$dataset$TaxonType == 'Arachnid' | taxa_cats$dataset$TaxonType == 'Invertebrate' | taxa_cats$dataset$TaxonType == 'Insect'))	#Number of invertebrates
		length(which(taxa_cats$matched.taxa$order == 'Coleoptera'))		#Number of beetles
		length(which(taxa_cats$matched.taxa$order == 'Lepidoptera'))	#Number of leps
		length(which(taxa_cats$matched.taxa$family == 'Formicidae'))	#Number of ants
		length(which(taxa_cats$dataset$TaxonType == 'Arachnid'))		#Number of spiders
	#Survey information
		sum(rowSums(!is.na(taxaXdata)) > 1)				#Number of taxa in >1 survey
		sum(rowSums(!is.na(taxaXdata)) > 1) / nrow(taxaXdata)	#As a proportion
		visits <- repeat.visits(thresh.data)
		sum(visits$mean.visits > 1)					#Number surveys with >1 visit
		sum(visits$mean.visits > 1)/nrow(visits)	#as proportion
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
		nrow(turn_points$dataset)
	#Number of functional groups
		nrow(fitted_func)	


#Ecological thresholds
	estimate.thresholds(turn_points, func_points)
	#Thresholds in taxa peak rate
		taxa_rate <- turn_points$dataset$maxrate[!is.na(turn_points$dataset$maxrate)]
		break.points(density(taxa_rate)$x, density(taxa_rate)$y)
	#Thresholds in functional group peak rate
		func_rate <- func_points$maxrate[!is.na(func_points$maxrate)]
		break.points(density(func_rate)$x, density(func_rate)$y)


#Taxonomic categories
	#Number of orders
		length(unique(taxa_cats$matched.taxa$order))
	#Number of orders with impacted taxa
		func <- function(x) sum(!is.na(x))/length(x)
		reps <- by(turn_points$dataset$turn.point, factor(turn_points$dataset$Order), FUN = func)
		sum(reps > 0)				#Number of orders containing impacted taxa
		sum(reps > 0)/length(reps)	#As a proportion
	#Number of functional groups impacted
		fs <- by(func_points$turn.point, factor(func_points$taxon), FUN = func)
		sum(fs > 0)				#Number of impacted functional groups
		sum(fs > 0)/length(fs)	#As a proportion
	#Proportion taxa with turning points
		by(turn_points$dataset$turn.point, factor(turn_points$dataset$TaxonType), FUN = func)	#Proportion taxa with turning points
		by(turn_points$dataset$turn.point, factor(turn_points$dataset$TaxonType), FUN = mean, na.rm = TRUE) 	#Mean turning point

#Functional resilience
	resil_func <- resil.func(func_groups = func.groups, func_points = func_points, turns = turn_points2$dataset)
	resil_func$funcs[order(resil_func$funcs$resilience), c(17, 20,25)]	#Functional groups ordered from least to most resilient
	sort(by(resil_func$funcs$resilience, resil_func$funcs$category, mean))	#Mean resilience per functional category
	anova(lm(resilience ~ category, data = resil_func$funcs))

#Functional composition - example taxa
	#Body size
	find.egs(Function = 'BodyMass_Fish_Categorised', function_qual = 'low', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Habitat strata generalists
	find.egs(Function = 'StrataGeneralism_Categorised', function_qual = 'high', taxtype = NA, 
		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Trophic
#	find.egs(Function = 'TrophicGeneralism_Categorised', function_qual = 'high', taxtype = NA, 
#		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Diet
#	find.egs(Function = 'DietGeneralism_Categorised', function_qual = 'high', taxtype = NA, 
#		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
	#Parasites
#	find.egs(Function = 'TrophPara', function_qual = 'parasite', taxtype = NA, 
#		turn_points = turn_points, func_points = func_points, func_groups = func.groups)
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
	#Subterranean mammals
	find.egs(Function = 'Movement_Invertebrate', function_qual = 'legless', taxtype = 'mammal', 
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
	summarise.data(thresh.data)
	
#Table S2

	funcs <- rename.funcs(func_groups = func.groups, func_points = func_points)
	funcs_sub <- funcs[ , c('category', 'qualifier', 'TaxType')]
	funcs_sub$quantity <- NA
		funcs_sub$quantity[grep('high', funcs_sub$qualifier)] <- 'high'
		funcs_sub$quantity[grep('medium', funcs_sub$qualifier)] <- 'medium'
		funcs_sub$quantity[grep('low', funcs_sub$qualifier)] <- 'low'
		levels(funcs_sub$quantity) <- c('low', 'medium', 'high')
	funcs_sub$qualifier <- gsub('-', '', funcs_sub$qualifier)
	funcs_sub$qualifier <- gsub('high', '', funcs_sub$qualifier)
	funcs_sub$qualifier <- gsub('medium', '', funcs_sub$qualifier)
	funcs_sub$qualifier <- gsub('low', '', funcs_sub$qualifier)
	
	funcs_sub[order(funcs_sub$category, funcs_sub$qualifier, funcs_sub$quantity, funcs_sub$TaxType), ]
	
	