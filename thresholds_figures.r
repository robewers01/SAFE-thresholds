
#Set working environment
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

#Load packages
	require(paletteer)
	require(scales)

#Read in data
	fitted_thresh <- readRDS('results/fitted_thresh.rds')
		fitted_thresh <- fitted_thresh[!is.na(fitted_thresh$modtype), ]	#Remove taxa that weren't found for analyses
	turn_points <- readRDS('results/turn_points.rds')

#Summary calculations
	ecdf_out <- cdf(turn_points$turn.point)
	model_fits <- fitted.matrix(models = fitted_thresh)
	turn_points <- assign.taxon(dataset = turn_points$dataset)




##FIGURE 1: 

	par(mai = c(0.8, 0.7, 0.1, 0.3))
	par(oma = c(2, 2, 0, 1.5))
	xvals <- 0:100
	pal <- paletteer_d("ggthemes::excel_Aspect")


	#Panel A: Cumulative distribution function
		plot(ecdf_out$x, ecdf_out$prop, type = "l", lwd = 10 , col = alpha(pal[3], 0.7),
			ylim = c(0,0.5), 
			xlim = c(0, 100), ylab = '', 
			xlab = '', cex.lab = 2.5, cex.axis = 3, xaxt = 'n')
			axis(1, cex.axis = 3, padj = 1, cex.lab = 2.5)
		mtext('Proportion of taxa', side = 2, line = 4, cex = 2)
		mtext('Biomass reduction (%)', 1, line = 0, outer = TRUE, cex = 2)


	#Panel B: rates of change
		plot(x = xvals, y = log10(abs(rowMeans(model_fits$rates, na.rm = TRUE))), 
			#ylim=c(0,0.05),
			pch = 19, col = alpha(pal[4], 0.3), cex = 3,
			xlab = '', ylab = '', cex.lab = 2.5, cex.axis = 3, xaxt = 'n')
			axis(1, cex.axis = 3, padj = 1, cex.lab = 2.5)
		text(x = 5, y = 0.016, labels = c('(B)'), cex = 3)
		mtext('|Mean rate of change|', side = 2, line = 4, cex = 2)
	
	#Panel C: occupancy
		plot(0,0, xlim = c(0,100), ylim = c(0,1), col = 'white', 
			xlab = '', ylab = '', xaxt = 'n', cex.lab = 2.5, cex.axis = 3)
			axis(1, cex.axis = 3, padj = 1, cex.lab = 2.5, )
		for(i in 1:ncol(model_fits$obs)){
			lines(xvals, model_fits$obs[ ,i], col = alpha(pal[2], 0.05))
			}
		text(x = 5, y = 0.97, labels = c('(A)'), cex = 3)
		mtext('Occurrence probability', side = 2, line = 4, cex = 2)
	
		#Add occupancy by taxon
		taxtypes <- unique(turn_points$dataset$TaxonType)
			taxtypes <- taxtypes[!is.na(taxtypes)]
		taxon.names <- names(model_fits$obs)
		taxon.names <- gsub(".+?\\__", '', taxon.names)
		taxcats <- turn_points$matched.taxa$TaxonType[match(taxon.names, turn_points$matched.taxa$taxon_name)]
		for(i in 1:length(taxtypes)){
			#Extract data for that taxon group
			taxobs <- model_fits$obs[ , which(taxcats == taxtypes[i])]
			taxocc <- rowMeans(taxobs, na.rm = TRUE)
			lines(taxocc, lwd = 2, col = i)
			}

		#Add average occupancy
		occ <- rowMeans(model_fits$obs, na.rm = TRUE)
		#Add to plot
		lines(xvals, occ, lwd=10, col = alpha(pal[2], 0.9))
	
	
	#Panel D: slopes histogram
		slopes <- turn_points$dataset$slope
		slopes[slopes < 0] <- 0 - log10(abs(slopes[slopes < 0])+1)
		slopes[slopes > 0] <- log10(slopes[slopes > 0]+1)
		dens <- density(slopes[turn_points$dataset$pval < 0.05], na.rm = TRUE)

		plot(dens, xlim = c(-3, 1), 
			xaxt = 'n', main = "", xlab = '', ylab = '',
			lwd = 5, col = alpha(pal[5], 1), cex.lab = 2.5, cex.axis = 2)
		polygon(dens, col = alpha(pal[5], 0.3), border = NA)
		text(x = 5, y = 0.035, labels = c('(A)'), cex = 2)
		mtext('Density', side = 2, line = 3.4, cex = 2)
		mtext('Slope', 1, line = 3, cex = 2)
		
		axis(1, at=c(-3:1), labels = c(expression('-10'^2), expression('-10'^1),expression('-10'^0),0,expression('10'^0)))
		
	
	
	
	