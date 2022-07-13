
#Set working environment
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

#Load packages
	require(paletteer)

#Read in data
	fitted_thresh <- readRDS('results/fitted_thresh.rds')
	turn_points <- readRDS('results/turn_points.rds')

#Summary calculations
	ecdf_out <- cdf(turn_points$turn.point)


##FIGURE 1: 

	par(mai = c(0.8, 0.7, 0.1, 0.3))
	par(oma = c(2, 2, 0, 1.5))
	

	#Panel A: Cumulative distribution function
		plot(ecdf_out$x, ecdf_out$prop, type = "l", lwd = 10 , col = alpha(pal[3], 0.7),
			ylim = c(0,0.5), 
			xlim = c(0, 100), ylab = '', 
			xlab = '', cex.lab = 2.5, cex.axis = 3, xaxt = 'n')
			axis(1, cex.axis = 3, padj = 1, cex.lab = 2.5)
		mtext('Proportion of taxa', side = 2, line = 4, cex = 2)
		mtext('Biomass reduction (%)', 1, line = 0, outer = TRUE, cex = 2)


	#Panel B: rates of change
	
	
	#Function to create matrix of observations and rates of change for all taxa
fitted.matrix <- function(models, predx = c(0:100)){
	#models = list of models to make predictions for
	#predx = set of x-values to make predictions for

	rates <- obs <- matrix(NA, nrow = 101, ncol = length(fitted_thresh$models))
	for(i in 1:length(fitted_thresh$models)){
		target <- fitted_thresh$models[[i]]
		fitted <- fitted.vals(target, predx = predx, mod_coef = NA)$fits
		rates[ , i] <- fitted$d1
		obs[ , i] <- fitted$obs
		}
	rownames(rates) <- rownames(obs) <- predx
	return(obs = obs, rates = rates)
	}
	
	
	
	
	
	
	
	
	
	