
#Set working environment
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")


#Read in and summarise data
	turn_points <- readRDS('results/turn_points.rds')
	turn_points <- assign.taxon(dataset = turn_points, taxon_table = taxa)
	func.groups <- readRDS('data/functional_groups.rds')			#Functional groups
	func_points <- readRDS('results/func_points.rds')




#Correlation between sample size and probability of impact
	#Taxonomic groups
		func <- function(x) sum(!is.na(x))/length(x)
		prop_imp <- by(turn_points$dataset$turn.point, factor(turn_points$dataset$TaxonType), FUN = func)	#Proportion taxa with turning points
		num_tax <- by(turn_points$dataset$turn.point, factor(turn_points$dataset$TaxonType), FUN = length)	#Number of taxa per taxonomic group
		cor.test(prop_imp, num_tax)
	
	#Functional groups
		resil_func <- resil.func(func_groups = func.groups, func_points = func_points, turns = turn_points$dataset)
		cor.test(resil_func$susc$prop_imp, resil_func$susc$num.taxa)
		

#Fig 3A vs 3B
	resil_func$susc[rownames(resil_func$susc) == 'StraArb_arboreal', ]
	#which of the 127 taxa had individual models fitted?	
	taxa_func <- taxaXfunc(func_groups = func.groups)
	taxa_arb <- taxa_func$taxon.lists[which(names(taxa_func$taxon.lists) == 'StraArb_arboreal')]	#full list of arboreal taxa

	ind <- match(taxa_arb$StraArb_arboreal, turn_points$dataset$taxon)
	ind <- ind[!is.na(ind)]
	fitted_arb <- turn_points$dataset[ind, ]			#Model results for arboreal taxa
	
	sum(is.na(fitted_arb$turn.point))					#Number arboreal taxa with turning point estimates
	
	
	
