
#Set working environment
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

#Load packages
	require(paletteer)
	require(scales)

#Read in data
	func.groups <- readRDS('data/functional_groups.rds')			#Functional groups
	fitted_thresh <- readRDS('results/fitted_thresh.rds')
		fitted_thresh <- fitted_thresh[!is.na(fitted_thresh$modtype), ]	#Remove taxa that weren't found for analyses
	fitted_func <- readRDS('results/fitted_func.rds')
		fitted_func <- fitted_func[!is.na(fitted_func$num.occs), ]	#Remove taxa that weren't found for analyses
	turn_points <- readRDS('results/turn_points.rds')
	func_points <- readRDS('results/func_points.rds')
	break_points <- readRDS('results/break_points.rds')

#Summary calculations
	ecdf_out <- cdf(turn_points$turn.point)
	ecdf_func <- cdf(func_points$turn.point)
	model_fits <- fitted.matrix(models = fitted_thresh)
	func_fits <- fitted.matrix(models = fitted_func)
	turn_points <- assign.taxon(dataset = turn_points)


##FIGURE 1: MAIN RESULTS

png('figures/fig1.png', width = 800, height = 800)
{
	par(mfrow = c(2,2))
	par(mai = c(0.8, 1.0, 0.3, 0.3))
	par(oma = c(3, 0, 0, 0))
	xvals <- 0:100
	pal <- paletteer_d("ggthemes::excel_Aspect")

	#Panel A: Cumulative distribution function
	{
		#Taxa
		plot(ecdf_out$x, ecdf_out$prop, type = "l", lwd = 5 , col = pal[4],
			ylim = c(0,0.8), 
			xlim = c(0, 100), ylab = '', 
			xlab = '', cex.lab = 2.5, cex.axis = 2)
		#Functional groups
		lines(ecdf_func$x, ecdf_func$prop, lwd = 5 , col = pal[1])
		
		#Thresholds
	#	plot.breaks(x = ecdf_out$x, y = ecdf_out$prop, pch = 21, cex = 5, lwd = 2 )
	#	plot.breaks(ecdf_func$x, ecdf_func$prop, pch = 21, cex = 5, lwd = 2 )
		
		legend('bottomleft', legend = c('taxa', 'functional groups'), 
			pch = c(22,22),
			pt.bg = c(alpha(pal[4], 0.3), alpha(pal[1],0.3)),
			col = c(pal[4], pal[1]), 
			text.col=c(pal[4], pal[1]),
			cex = 2, pt.cex = 4, lty = 0, bty = 'n', lwd = 2,
			x.intersp = 0.1)

#		legend('bottomleft', legend = c('taxa', 'groups', 'accel.', 'decel.'), 
#			pch = c(22,22, 21, 21),
#			pt.bg = c(alpha(pal[4], 0.3), alpha(pal[1],0.3), 'white', alpha(pal[3], 0.8)),
#			col = c(pal[4], pal[1], 1, 1), 
#			text.col=c(pal[4], pal[1], 1, pal[3]),
#			cex = 2, pt.cex = 4, lty = 0, bty = 'n', lwd = 2)
		
		text(x = 0, y = 0.78, labels = c('(A) Cumulative impact'), cex = 2, pos = 4)
		mtext('Proportion of taxa', side = 2, line = 4, cex = 2)
#		mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	}

	#Panel B: mean occurrence
	{
		plot(0,0, xlim = c(0,100), ylim = c(0,1), col = 'white', 
			xlab = '', ylab = '',  cex.lab = 2.5, cex.axis = 2)

		#Individual fits
		for(i in 1:ncol(model_fits$obs)){
			lines(x = xvals, y = model_fits$obs[,i], col = alpha('grey', 0.2))
			}
		for(i in 1:ncol(func_fits$obs)){
			lines(x = xvals, y = func_fits$obs[,i], col = alpha('grey', 0.2))
			}
		
		#Taxa
		taxon_occ <- apply(model_fits$obs, MARGIN = 1, FUN = mean, na.rm = TRUE)
		lines(xvals, taxon_occ, lwd=5, col = pal[4])
	#	confints <- apply(model_fits$obs, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
	#	polygon(x = c(xvals, rev(xvals)), y = c(confints[1,], rev(confints[2,])),
	#		col = alpha(pal[4], 0.3))

		#Functional groups
		func_occ <- apply(func_fits$obs, MARGIN = 1, FUN = mean, na.rm = TRUE)
		lines(xvals, func_occ, lwd=5, col = pal[1])
	#	confints <- apply(func_fits$obs, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
	#	polygon(x = c(xvals, rev(xvals)), y = c(confints[1,], rev(confints[2,])),
	#		col = alpha(pal[1], 0.3))

		#Thresholds
	#	plot.breaks(x = xvals, y = taxon_occ, pch = 21, cex = 5, lwd = 2 )
	#	plot.breaks(x = xvals, y = func_occ, pch = 21, cex = 5, lwd = 2 )

		text(x = 0, y = 0.97, labels = c('(B) Mean occurrence'), cex = 2, pos = 4)
		mtext('Occurrence probability', side = 2, line = 4, cex = 2)
#		mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	}

	#Panel C: turning points density
	{
		#Taxa
		taxa_turn <- turn_points$turn.point[!is.na(turn_points$turn.point)]
		plot(density(taxa_turn), xlim = c(0,100), main = "", xlab = '', ylab = '',
			ylim = c(0, 0.03), 
			lwd = 5, col = alpha(pal[4], 1), cex.lab = 2.5, cex.axis = 2)
		polygon(density(taxa_turn), col = alpha(pal[4], 0.3), angle = -45, 
			border = pal[4], lwd = 5, density = 15)
		plot.breaks(density(taxa_turn)$x, density(taxa_turn)$y, add_plot = TRUE,
			pch = 21, cex = 4, col = pal[4], lwd = 2, decel.col = pal[4])
		#Functional groups
		func_turn <- func_points$turn.point[!is.na(func_points$turn.point)]
		polygon(density(func_turn), col = alpha(pal[1], 0.3), angle = 45, 
			border = pal[1], lwd = 5, density = 15)
		plot.breaks(density(func_turn)$x, density(func_turn)$y, add_plot = TRUE,
			pch = 21, cex = 4, col = pal[1], lwd = 2, decel.col = pal[1])
	
		legend('right', legend = c('accel.', 'decel.'), 
			pch = c(21, 21), pt.bg = c('white', alpha(1, 0.5)),
			cex = 2, pt.cex = 4, lty = 0, bty = 'n', lwd = 2,
			x.intersp = 0.1)

		text(x = 0, y = 0.029, labels = c('(C) Turning points'), cex = 2, pos = 4)
		mtext('Density', side = 2, line = 4, cex = 2)
#		mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	}
	
	#Panel D: rates of change density
	{
	#par(mai = c(0.8, 0.8, 0.3, 0.3))
	#par(oma = c(2, 2, 0, 0))
		#Taxa
		taxa_rate <- turn_points$maxrate[!is.na(turn_points$maxrate)]
		plot(density(taxa_rate), xlim = c(0,100), main = "", xlab = '', ylab = '',
			ylim = c(0, 0.03),
			lwd = 5, col = alpha(pal[4], 1), cex.lab = 2.5, cex.axis = 2)
		polygon(density(taxa_rate), col = alpha(pal[4], 0.3), angle = -45, 
			border = pal[4], lwd = 5, density = 15)
		plot.breaks(density(taxa_rate)$x, density(taxa_rate)$y, add_plot = TRUE,
			pch = 21, cex = 4, lwd = 2, col = pal[4], decel.col = pal[4])
		#Functional groups
		func_rate <- func_points$maxrate[!is.na(func_points$maxrate)]
		polygon(density(func_rate), col = alpha(pal[1], 0.3), angle = 45, 
			border = pal[1], lwd = 5, density = 15)
		plot.breaks(density(func_rate)$x, density(func_rate)$y, add_plot = TRUE,
			pch = 21, cex = 4, lwd = 2, col = pal[1], decel.col = pal[1])
	
		text(x = 0, y = 0.029, labels = c('(D) Peak change points'), cex = 2, pos = 4)
		mtext('Density', side = 2, line = 4, cex = 2)
#		mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	}
	
	mtext('Biomass reduction (%)', 1, line = 0, outer = TRUE, cex = 2)

	}
	dev.off()
		
	

##FIGURE 3 - FUNCTIONAL THRESHOLDS
png('figures/fig3.png', width = 320, height = 800)
{
	par(mai = c(0.8, 0.1, 0.1, 0.1))
	par(oma = c(0, 0, 0, 0))
	#pal <- paletteer_d("ggthemes::excel_Aspect")
	pal <- paletteer_d("Redmonder::qPBI")
	
	#Get and filter data
	funcs <- rename.funcs(func_groups = func.groups, func_points = func_points)
	funcs <- funcs[!is.na(funcs$turn.point) ,]
	funcs <- funcs[!is.na(funcs$category), ]
	funcs <- funcs[!duplicated(funcs), ]
	#Exclude groups with impacts that only start/end at the ends of the gradient
	funcs <- funcs[-which(funcs$turn.point == funcs$maxrate), ]
	
	#Sort data
	funcs$category <- factor(funcs$category, levels = c('Red List status', 'Habitat strata', 
		'Plant', 'Physiology', 'Development', 'Sociality', 'Movement', 'Diet', 'Trophic',  'Body mass'))
	funcs$qualifier <- factor(funcs$qualifier, levels = rev(c('high', 'medium', 'low', 
		'generalism-high','generalism-low',
		'aerial', 'arboreal', 'understory', 'terrestrial', 'subterranean', 'aquatic', 'constant', 'variable',
		'endotherm', 'ectotherm', 'direct', 'indirect', 'social', 'pair', 'solitary', 
		'winged', 'legless', 'threatened', 'not threatened', 
		'parasitoid', 'parasite', 'predator', 'herbivore', 'producer', 'saprophage', 
		'hematophage', 'vertivore', 'invertivore', 'florivore', 'folivore', 'frugivore', 
		'granivore', 'nectarivore', 'palynivore', 'rhizophage', 'xylophage', 'mycophage', 
		'necrophage', 'detritivore',
		rev(funcs$qualifier[funcs$TaxType == 'plant']))))
	
	funcs <- funcs[order(funcs$category, funcs$qualifier, funcs$TaxType, funcs$turn.point), ]
	
	#Add colours
	funcs$col <- pal[as.numeric(factor(funcs$category))]
	#Add line types
	funcs$lty <- 1
	funcs$lty[funcs$slope < 0] <- 2
	#symbol types
	types <- c('all taxa', 'invertebrate', 'amphibian', 'bird', 'fish', 'mammal', 'plant')
	pchtypes <- c(21:25,8,12)
	pchmatch <- data.frame(cbind(types, pchtypes))
	funcs$pch <- as.numeric(pchmatch$pchtypes[match(funcs$TaxType, pchmatch$types)])

	#plot figure
	plot(0,0, xlim = c(-85,180), ylim = c(1, nrow(funcs) + 3), col = 'white', 
		xlab = '', ylab = '',  yaxt = 'n', xaxt = 'n', cex.lab = 2.5, cex.axis = 2)
	axis(1, at = c(0, 50, 100), cex.axis = 2)

	#Add categories
	cats <- unique(funcs$category)
		cats <- cats[!is.na(cats)]
	for(k in 1:length(cats)){
		ymin <- min(which(funcs$category == cats[k])) - 0.5
		ymax <- max(which(funcs$category == cats[k])) + 0.5
		ymid <- mean(c(ymin, ymax))
		polycol <- funcs$col[which(funcs$category == cats[k])[1]]
		polygon(x = c(-100, -100, 200, 200), y = c(ymin, ymax, ymax, ymin),
			col = alpha(polycol, 0.1), density=100, border = NA)
		text(x = 145, y = ymid, cats[k], cex = 1, adj = c(0.5))
		}
	#Add data
	for(i in 1:nrow(funcs)){
		points(x = funcs$turn.point[i], y = i, pch = funcs$pch[i], col = funcs$col[i], lty = funcs$lty[i])
		if(!is.na(funcs$maxrate[i]) & !(funcs$maxrate[i] == 100 & funcs$turn.point[i] == 100)){
			points(x = funcs$maxrate[i], y = i, pch = funcs$pch[i], col = funcs$col[i], lty = funcs$lty[i])
			lines(x = c(funcs$turn.point[i], funcs$maxrate[i]), y = rep(i,2), col = funcs$col[i], lty = funcs$lty[i])
			}
		text(x = -5, y = i, funcs$qualifier [i], cex = 0.8, adj = c(1, 0.5))
		}

	legend('top', legend = c('increasing', 'decreasing', pchmatch$types), 
		pch = c(NA, NA, as.numeric(pchmatch$pchtypes)),
		lty = c(1, 2, rep(0, nrow(pchmatch))),
		lwd = c(2,2, rep(1, nrow(pchmatch))),
		bty = 'n', pt.cex = 1.7, ncol = 3)
	abline(nrow(funcs) +0.5, 0, lwd = 1)
	mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	}
	dev.off()


















###########################
##NOT USED###########

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






	#Panel D: rates of change
	{
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
	



##FIGURE - TAXON THRESHOLDS
	
	taxa <- turn_points$dataset[!is.na(turn_points$dataset$turn.point) ,]
	taxa <- taxa[order(taxa$TaxonType, taxa$turn.point, taxa$maxrate),]
	plot(0,0, xlim = c(0,100), ylim = c(1, nrow(taxa)), col = 'white', 
		xlab = '', ylab = '',  yaxt = 'n', cex.lab = 2.5, cex.axis = 2)

	for(i in 1:nrow(taxa)){
		points(x = taxa$turn.point[i], y = i, pch = 21, col = as.numeric(taxa$TaxonType[i]))
		if(!is.na(taxa$maxrate[i])){
			lines(x = c(taxa$turn.point[i], taxa$maxrate[i]), y = rep(i,2), col = as.numeric(taxa$TaxonType[i]))
			points(x = taxa$maxrate[i], y = i, pch = 19, col = as.numeric(taxa$TaxonType[i]))
			}
		}
	
