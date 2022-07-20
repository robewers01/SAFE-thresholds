
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
	fitted_func <- readRDS('results/fitted_func.rds')
		fitted_func <- fitted_func[!is.na(fitted_func$num.occs), ]	#Remove taxa that weren't found for analyses
	turn_points <- readRDS('results/turn_points.rds')
	func_points <- readRDS('results/func_points.rds')

#Summary calculations
	ecdf_out <- cdf(turn_points$turn.point)
	ecdf_func <- cdf(func_points$turn.point)
	model_fits <- fitted.matrix(models = fitted_thresh)
	func_fits <- fitted.matrix(models = fitted_func)
	turn_points <- assign.taxon(dataset = turn_points)




##FIGURE 1: 

png('figures/fig1.png'), width = 800, height = 800)
{
	par(mfrow = c(2,2))
	par(mai = c(0.8, 0.8, 0.3, 0.3))
	par(oma = c(2, 2, 0, 0))
	xvals <- 0:100
	pal <- paletteer_d("ggthemes::excel_Aspect")


	#Panel A: slopes density
	{
		#Taxon data
		slopes <- as.numeric(fitted_thresh$slope)
			slopes <- slopes[!is.na(slopes) & fitted_thresh$pval < 0.05]
		slopes[slopes < 0] <- 0 - log10(abs(slopes[slopes < 0])+1)
		slopes[slopes > 0] <- log10(slopes[slopes > 0]+1)
		dens <- density(slopes[fitted_thresh$pval < 0.05], na.rm = TRUE)
		plot(dens, xlim = c(-2, 1), ylim = c(0,40), 
			xaxt = 'n', main = "", xlab = '', ylab = '',
			lwd = 2, col = alpha(pal[4], 1), cex.lab = 2.5, cex.axis = 2)
		polygon(dens, col = alpha(pal[4], 0.3), angle = -45, border = pal[4], lwd = 2)
		
		#Functional groups data
		slopes <- as.numeric(fitted_func$slope)
			slopes <- slopes[!is.na(slopes) & fitted_func$pval < 0.05]
		slopes[slopes < 0] <- 0 - log10(abs(slopes[slopes < 0])+1)
		slopes[slopes > 0] <- log10(slopes[slopes > 0]+1)
		dens <- density(slopes, na.rm = TRUE)
		polygon(dens, col = alpha(pal[1], 0.3), angle = -45, border = pal[1], lwd = 2)
			
		text(x = 5, y = 0.035, labels = c('(A)'), cex = 2)
		mtext('Density', side = 2, line = 3.4, cex = 2)
		mtext('Slope', 1, line = 3, cex = 2)
		axis(1, at=c(-3:1), labels = c(expression('-10'^2), expression('-10'^1),expression('-10'^0),0,expression('10'^0)),
			cex.lab = 2.5, cex.axis = 2)

		text(x = 0.85, y = 38, labels = c('(A)'), cex = 3)
		legend('topleft', legend = c('taxa', 'groups'), pch = 22, 
			pt.bg = c(alpha(pal[4], 0.3),alpha(pal[1],0.3)), col = c(pal[4], pal[1]), 
			text.col=c(pal[4], pal[1]),
			cex = 2, pt.cex = 4, lty = 0, bty = 'n', lwd = 2)
	}

	#Panel B: mean occurrence
	{
		plot(0,0, xlim = c(0,100), ylim = c(0,1), col = 'white', 
			xlab = '', ylab = '',  cex.lab = 2.5, cex.axis = 2)

		#Taxa
		taxon_occ <- apply(model_fits$obs, MARGIN = 1, FUN = mean, na.rm = TRUE)
		lines(xvals, taxon_occ, lwd=10, col = pal[4])
		confints <- apply(model_fits$obs, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
		polygon(x = c(xvals, rev(xvals)), y = c(confints[1,], rev(confints[2,])),
			col = alpha(pal[4], 0.3))

		#Functional groups
		func_occ <- apply(func_fits$obs, MARGIN = 1, FUN = mean, na.rm = TRUE)
		lines(xvals, func_occ, lwd=10, col = pal[1])
		confints <- apply(func_fits$obs, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
		polygon(x = c(xvals, rev(xvals)), y = c(confints[1,], rev(confints[2,])),
			col = alpha(pal[1], 0.3))

		#Tthresholds
		plot.breaks(x = xvals, y = taxon_occ, pch = 21, cex = 5, lwd = 2 )
		plot.breaks(x = xvals, y = func_occ, pch = 21, cex = 5, lwd = 2 )

		text(x = 95, y = 0.97, labels = c('(B)'), cex = 3)
		mtext('Occurrence probability', side = 2, line = 4, cex = 2)
		mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	}


	#Panel C: Cumulative distribution function
	{
		#Taxa
		plot(ecdf_out$x, ecdf_out$prop, type = "l", lwd = 10 , col = pal[4],
			ylim = c(0,0.6), 
			xlim = c(0, 100), ylab = '', 
			xlab = '', cex.lab = 2.5, cex.axis = 2)
		#Functional groups
		lines(ecdf_func$x, ecdf_func$prop, lwd = 10 , col = pal[1])
		
		#Thresholds
		plot.breaks(x = ecdf_out$x, y = ecdf_out$prop, pch = 21, cex = 5, lwd = 2 )
		plot.breaks(ecdf_func$x, ecdf_func$prop, pch = 21, cex = 5, lwd = 2 )
		
		legend('bottomleft', legend = c('taxa', 'groups', 'accel.', 'decel.'), 
			pch = c(22,22, 21, 21),
			pt.bg = c(alpha(pal[4], 0.3), alpha(pal[1],0.3), 'white', alpha(pal[3], 0.8)),
			col = c(pal[4], pal[1], 1, 1), 
			text.col=c(pal[4], pal[1], 1, pal[3]),
			cex = 2, pt.cex = 4, lty = 0, bty = 'n', lwd = 2)
		
		text(x = 95, y = 0.58, labels = c('(B)'), cex = 3)
		mtext('Proportion of taxa', side = 2, line = 4, cex = 2)
		mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	}

	#Panel D: rates of change
		plot(0,0, xlim = c(0,100), ylim = c(0,0.03), col = 'white', 
			xlab = '', ylab = '',  cex.lab = 2.5, cex.axis = 2)

		#Taxa
		taxon_occ <- apply(abs(model_fits$rate), MARGIN = 1, FUN = mean, na.rm = TRUE)
		lines(xvals, taxon_occ, lwd=10, col = pal[4])
		confints <- apply(abs(model_fits$rate), MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
		polygon(x = c(xvals, rev(xvals)), y = c(confints[1,], rev(confints[2,])),
			col = alpha(pal[4], 0.3))

		#Functional groups
		func_occ <- apply(abs(func_fits$rate), MARGIN = 1, FUN = mean, na.rm = TRUE)
		lines(xvals, func_occ, lwd=10, col = pal[1])
		confints <- apply(abs(func_fits$rate), MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
		polygon(x = c(xvals, rev(xvals)), y = c(confints[1,], rev(confints[2,])),
			col = alpha(pal[1], 0.3))

		text(x = 95, y = 0.029, labels = c('(D)'), cex = 3)
		mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
		mtext('|Mean rate of change|', side = 2, line = 4, cex = 2)

		#Thresholds
		try(plot.breaks(x = xvals, y = func_occ, pch = 21, cex = 5, lwd = 2 ), silent = TRUE)
		try(plot.breaks(x = xvals, y = taxon_occ, pch = 21, cex = 5, lwd = 2 ), silent = TRUE)

	
	}
	dev.off()
		
	
	
	
	