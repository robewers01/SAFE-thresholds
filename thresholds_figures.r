
#Set working environment
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

#Load packages
	require(ape)
	require(grid)
	require(jpeg)
	require(lme4)
	require(paletteer)
	require(png)
	require(scales)

#Read in data
	thresh_data <- readRDS('data/threshold_taxa_data.rds')
	func.groups <- readRDS('data/functional_groups.rds')			#Functional groups
	fitted_thresh <- readRDS('results/fitted_thresh.rds')
		fitted_thresh <- fitted_thresh[!is.na(fitted_thresh$modtype), ]	#Remove taxa that weren't found for analyses
	fitted_func <- readRDS('results/fitted_func.rds')
		fitted_func <- fitted_func[!is.na(fitted_func$num.occs), ]	#Remove taxa that weren't found for analyses
	turn_points <- readRDS('results/turn_points.rds')
	func_points <- readRDS('results/func_points.rds')
	break_points <- readRDS('results/break_points.rds')
	taxa <- readRDS('data/taxon_table.rds')					#Full list of all taxa
	map <- read.table('data/species_families_order_map.txt', sep = '-')	#Identified one family and example species per order that exists on TimeTree.org
	tr <- read.tree("data/species.nwk")			#Imported phylogeny from TimeTree (www.timetree.org)
		bayes_results <- readRDS("data/bayes.rds")[[1]]			#Results from Bayesian slopes analysis (Replicability analysis)

#Summary calculations
#	full_funcs <- expand.funcs(func.groups)
	ecdf_out <- cdf(turn_points$turn.point)
	ecdf_func <- cdf(func_points$turn.point)
	model_fits <- fitted.matrix(models = fitted_thresh)
	func_fits <- fitted.matrix(models = fitted_func)
	turn_points <- assign.taxon(dataset = turn_points, taxon_table = taxa)
	resil <- resil.dat(turns = turn_points$dataset, bayes = bayes_results)


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
		
		legend('bottomleft', legend = c('taxa', 'functional groups'), 
			pch = c(22,22),
			pt.bg = c(alpha(pal[4], 0.3), alpha(pal[1],0.3)),
			col = c(pal[4], pal[1]), 
			text.col=c(pal[4], pal[1]),
			cex = 2, pt.cex = 4, lty = 0, bty = 'n', lwd = 2,
			x.intersp = 0.1)

		text(x = 0, y = 0.78, labels = c('(A) Cumulative impact'), cex = 2, pos = 4)
		mtext('Proportion of taxa', side = 2, line = 4, cex = 2)
	}

	#Panel B: mean occurrence
	{
		plot(0,0, xlim = c(0,100), ylim = c(0,1), col = 'white', 
			xlab = '', ylab = '',  cex.lab = 2.5, cex.axis = 2)

		#Individual fits
		for(i in 1:ncol(model_fits$obs)){
			lines(x = xvals, y = model_fits$obs[,i], col = alpha(pal[4], 0.1))
			}
		for(i in 1:ncol(func_fits$obs)){
			lines(x = xvals, y = func_fits$obs[,i], col = alpha(pal[1], 0.1))
			}
		
		#Taxa
		taxon_occ <- apply(model_fits$obs, MARGIN = 1, FUN = mean, na.rm = TRUE)
		lines(xvals, taxon_occ, lwd=5, col = pal[4])

		#Functional groups
		func_occ <- apply(func_fits$obs, MARGIN = 1, FUN = mean, na.rm = TRUE)
		lines(xvals, func_occ, lwd=5, col = pal[1])

		text(x = 0, y = 0.97, labels = c('(B) Mean occurrence'), cex = 2, pos = 4)
		mtext('Occurrence probability', side = 2, line = 4, cex = 2)
	}

	#Panel C: turning points density
	{
		#Taxa
		taxa_turn <- turn_points$dataset$turn.point[!is.na(turn_points$dataset$turn.point)]
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
	}
	
	#Panel D: rates of change density
	{
		#Taxa
		taxa_rate <- turn_points$dataset$maxrate[!is.na(turn_points$dataset$maxrate)]
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
	}
	
	mtext('Biomass reduction (%)', 1, line = 0, outer = TRUE, cex = 2)

	}
	dev.off()
		


##FIGURE 2 - TAXONOMIC SENSITIVITY
png('figures/fig2.png', width = 800, height = 400)
{

	par(mfrow = c(1,2))
	par(oma = c(0, 0, 0, 0))
	par(mai = c(1, 1.0, 0.3, 0.3))
	pal <- paletteer_d("ggthemes::excel_Aspect")

	plant<-readPNG("data/leaf.png")
	insect<-readPNG("data/invert.png")
	mammal<-readPNG("data/mammal.png")
	bird<-readPNG("data/aves.png")
	reptile<-readPNG("data/reptile.png")
	amphib<-readPNG("data/frog.png")
	fish<-readPNG("data/fish.png")

	#Panel A: Proportion taxa with turning points
	{
		turn_points$dataset$TaxonType <- factor(turn_points$dataset$TaxonType, levels = c('plant', 'invertebrate', 
			'mammal', 'bird', 'reptile', 'amphibian', 'fish'))
		func <- function(x) sum(!is.na(x))/length(x)
		prop_imp <- by(turn_points$dataset$turn.point, factor(turn_points$dataset$TaxonType), FUN = func)	#Proportion taxa with turning points
		mean.turn <- by(turn_points$dataset$turn.point, factor(turn_points$dataset$TaxonType), FUN = mean, na.rm = TRUE) 	#Mean turning point
		
		#Plot figure
		plot(mean.turn, prop_imp, cex.axis = 2,
			xlim = c(0,40), ylim = c(0.3,0.9), xlab = '', ylab = '')
		
		x.scale <- function(x) 0.11 + 0.342*(x/40)
		y.scale <- function(x) 0.19 + 0.73* ((x-0.3)/(0.9-0.3))
			
		#Plot images
		grid.raster(plant, x=x.scale(mean.turn[1]), y=y.scale(prop_imp[1]), width = 0.05)
		grid.raster(insect, x=x.scale(mean.turn[2]), y=y.scale(prop_imp[2]), width = 0.05)
		grid.raster(mammal, x=x.scale(mean.turn[3]), y=y.scale(prop_imp[3]), width = 0.05)
		grid.raster(bird, x=x.scale(mean.turn[4]), y=y.scale(prop_imp[4]), width = 0.05)
		grid.raster(reptile, x=x.scale(mean.turn[5]), y=y.scale(prop_imp[5]), width = 0.05)
		grid.raster(amphib, x=x.scale(mean.turn[6]), y=y.scale(prop_imp[6]), width = 0.05)
		grid.raster(fish, x=x.scale(mean.turn[7]), y=y.scale(prop_imp[7]), width = 0.05)

		#Add arrows
		arrows(x = 35, y = 0.9, x1 = 5, y1 = 0.9, lwd = 2, angle = 20, length = 0.15, col = pal[6])
		text(x = 20, y = 0.88, labels = 'increasing sensitivity', adj = 0.5, cex = 1.5, col = pal[6])
		arrows(x = 40, y = 0.5, x1 = 40, y1 = 0.9, lwd = 2, angle = 20, length = 0.15, col = pal[2])
		text(x = 38, y = 0.7, labels = 'increasing susceptibility', adj = 0.5, cex = 1.5, srt = 270, col = pal[2])
		arrows(x=7, y = 0.71, x1 = 30, y1 = 0.45, lwd=3, angle = 20, length = 0.15, col = pal[4])
		text(x = 17, y = 0.57, labels = 'increasing resilience', adj = 0.5, cex = 1.5, srt = -38, col = pal[4])

		text(x = -1, y = 0.88, labels = c('(A)'), cex = 2, pos = 4)
		mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
		mtext('Proportion taxa impacted', 2, line = 3, cex = 2)
	}
	
	#Panel B: Connection with replicability
	{
		plot(resil$bayes_prob, resil$resilience, cex.axis = 2,
			 xlim = c(0, 0.15), ylim = c(0.2, 1), xlab = '', ylab = '')
		
		#Add fitted line
		resil_mod <- lm(resilience ~ bayes_prob, data = resil)
		vals <- coef(resil_mod)
#		abline(vals[1], vals[2], lwd = 2, lty = 2, col = 'black')
		#Add confidence intervals
		preds <- as.data.frame(predict(resil_mod, newdata = data.frame(bayes_prob = seq(-0.05, 0.2, 0.01)), interval = 'confidence'))
		preds$x <- seq(-0.05, 0.2, 0.01)
#		polygon(x = c(preds$x, rev(preds$x)), y = c(preds$lwr, rev(preds$upr)),
#			col = alpha(1, 0.1), density = 200)
	
		#Plot images
#		x.scale <- 	function(x) 0.61 + (0.95-0.61)*(((x-0.85) + (1/0.15)*(x-0.85))*0.85)
		x.scale <- 	function(x) 0.615 + (0.95-0.615)*(x/0.15)
		y.scale <- 	function(x) 0.2 + 0.72* ((x-0.2)/(1 - 0.2))

		grid.raster(plant, x=x.scale(resil$bayes_prob[1]), y=y.scale(resil$resilience[1]), width=0.05)
		grid.raster(insect, x=x.scale(resil$bayes_prob[2]), y=y.scale(resil$resilience[2]), width=0.05)
		grid.raster(mammal, x=x.scale(resil$bayes_prob[3]), y=y.scale(resil$resilience[3]), width=0.05)
		grid.raster(bird, x=x.scale(resil$bayes_prob[4]), y=y.scale(resil$resilience[4]), width=0.05)
		grid.raster(reptile, x=x.scale(resil$bayes_prob[5]), y=y.scale(resil$resilience[5]), width=0.05)
		grid.raster(amphib, x=x.scale(resil$bayes_prob[6]), y=y.scale(resil$resilience[6]), width=0.05)
		grid.raster(fish, x=x.scale(resil$bayes_prob[7]), y=y.scale(resil$resilience[7]), width=0.05)

		box()
		text(x = -0.005, y = 0.97, labels = c('(B)'), cex = 2, pos = 4)
		mtext('Context dependence', 1, line = 3, cex = 2)
		mtext('Resilience', 2, line = 3, cex = 2)
	}	
	dev.off()
	}
	


##FIGURE 3 - FUNCTIONAL THRESHOLDS
png('figures/fig3.png', width = 320, height = 800)
{
	par(mai = c(0.8, 0.1, 0.1, 0.1))
	par(oma = c(0, 0, 0, 0))
	#pal <- paletteer_d("ggthemes::excel_Aspect")
	pal <- paletteer_d("Redmonder::qPBI")
	
	#Get and filter data

	funcs <- arrange.funcplot(func_groups = func.groups, func_points = func_points, pal = pal)

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
		#	points(x = funcs$maxrate[i], y = i, pch = funcs$pch[i], col = funcs$col[i], lty = funcs$lty[i])
			arrows(x0 = funcs$turn.point[i], y0 = i, x1 = funcs$maxrate[i], y1 = i, 
				col = funcs$col[i], lty = funcs$lty[i],
				angle = 20, length =0.07)
			}
		text(x = -5, y = i, funcs$qualifier [i], cex = 0.8, adj = c(1, 0.5))
		}





#Add resilience data
	resil_func <- resil.func(func_groups = func.groups, func_points = func_points, turns = turn_points$dataset)


	
	keeps <- match(rownames(resil_func), funcs$taxon)
	resil_keeps <- resil_func[!is.na(keeps), ]
	funcs$resilience <- NA
	funcs$resilience[match(rownames(resil_keeps), funcs$taxon)] <- resil_keeps$resilience
	
	barplot(funcs$resilience * 100, horiz = TRUE, add = TRUE, width = 0.85, col = alpha('grey', 0.2))



require(emmeans)
require(lmtest)


	resil_func <- resil.func(func_groups = func.groups, func_points = func_points, turns = turn_points$dataset)
	funcs <- arrange.funcplot(func_groups = func.groups, func_points = func_points)
	resil_func <- merge(resil_func, funcs)
	resil_func$category <- factor(resil_func$category, levels = c('Red List status', 
		'Body mass', 'Physiology', 'Development', 'Movement', 'Sociality', 
		'Trophic', 'Diet', 'Habitat strata', 'Plant'))
	cats <- summary(resil_func$category)
	keeps <- names(cats)[cats >= 5]			#Only retain categories with >=5 estimates
	resil_dat <- resil_func[which(!is.na(match(resil_func$category, keeps))), ]
	resil_dat$category <- factor(resil_dat$category, levels = keeps)	
	
	#Convert zeros and ones for use in betareg
	resil_dat[, c(2:4)][resil_dat[, c(2:4)] == 0] <- 0.0001
	resil_dat[, c(2:4)][resil_dat[, c(2:4)] == 1] <- 0.9999

	mod_resil <- betareg(resilience ~ category, weights = num.taxa, data = resil_dat)
	lrtest(mod_resil)
	contrast(emmeans(mod_resil, specs = ~category), method = 'pairwise')
	
	mod_sens <- betareg(sensitivity ~ category, weights = num.taxa, data = resil_dat)
	lrtest(mod_sens)
	mod_susc <- betareg(susceptibility ~ category, weights = num.taxa, data = resil_dat)
	lrtest(mod_susc)




	res <- by(resil_dat$resilience, resil_dat$category, mean)
	sens <- by(resil_dat$sensitivity, resil_dat$category, quantile, probs = c(0.025, 0.5, 0.975))
	susc <- by(resil_dat$susceptibility, resil_dat$category, quantile, probs = c(0.025, 0.5, 0.975))

	par(oma = c(0, 0, 0, 0))
	par(mai = c(1, 1.0, 0.3, 0.3))
	plot(0,0, xlim = c(0.2, 1), ylim = c(0.2, 1), col = 'white', 
		xlab = '', ylab = '', cex.lab = 2.5, cex.axis = 2)

	for(i in 1:length(sens)){
		xvals <- sens[i][[1]]
		yvals <- susc[i][[1]]
		arrows(xvals[1], yvals[2], xvals[3], yvals[2], length = 0.05, lwd = 2, angle = 90, code = 3, col = 'grey')
		arrows(xvals[2], yvals[1], xvals[2], yvals[3], length = 0.05, lwd = 2, angle = 90, code = 3, col = 'grey')
		points(xvals[2], yvals[2], pch = 19, cex = 5*res[i])
		}
	#Plot text on top
	for(i in 1:length(sens)){
		xvals <- sens[i][[1]]
		yvals <- susc[i][[1]]
		text(x = xvals[2]+0.03, y = yvals[2]+0.03, adj = 0, labels = names(sens)[i], cex = 1.2)
		}
	
	legend('bottomright', pch = 19, pt.cex = c(0.4, 0.6, 0.8)*5, legend = c(0.4, 0.6, 0.8),
		cex = 2, lty = 0, bty = 'n', x.intersp = 0.8)
	mtext('Sensitivity', 1, line = 3, cex = 2)
	mtext('Susceptibility', 2, line = 3, cex = 2)

	#still to do...
		#colours - use colours from categories in Fig. 3
		#add bacground semi-transparent box to text so that it shows through better?









	legend('top', legend = c('increasing', 'decreasing', pchmatch$types), 
		pch = c(NA, NA, as.numeric(pchmatch$pchtypes)),
		lty = c(1, 2, rep(0, nrow(pchmatch))),
		lwd = c(2,2, rep(1, nrow(pchmatch))),
		bty = 'n', pt.cex = 1.7, ncol = 3)
	abline(nrow(funcs) +0.5, 0, lwd = 1)
	mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	}
	dev.off()



##FIGURE S1 - PHYLOGENY
png(paste('figures/figS1.png', sep = ""), width = 500, height = 500)
{
	phylo <- arrange.phylo(timetree = tr, raw_data = thresh_data, taxa_safe = taxa,
		tt_map = map, coefs = fitted_thresh, palette_col = paletteer_d("ggsci::default_igv"))

	par(oma = c(1,0,0,0), mar=c(3, 0, 0, 0))
	plot.phylo(phylo$tr, show.tip.label = FALSE, edge.color = 'white',
		edge.width = 1, font = 0.8, x.lim = c(0,2300))
	phydataplot(phylo$numbers2$numTax * 200, phylo$tr, offset = 30, col = alpha(phylo$pal[phylo$col_index], 0.3), legend = 'none')	#Number of taxa
	phydataplot(phylo$numbers3$propTax * 200, phylo$tr, offset = 30, col = alpha(phylo$pal[phylo$col_index], 0.3), legend = 'none')	#Proportion of taxa modelled
	
	polygon(x = c(1520,2300, 2300, 1520), y = c(0.5, 0.5, 5.5, 5.5), col = alpha(phylo$pal[1], 0.1), border = NA)
		text(x = 2280, y = 2.5, labels = 'Actinopterygii', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(5.5, 5.5, 6.5, 6.5), col = alpha(phylo$pal[2], 0.1), border = NA)
		text(x = 2280, y = 5, labels = 'Amphibia', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(6.5, 6.5, 7.5, 7.5), col = alpha(phylo$pal[3], 0.1), border = NA)
		text(x = 2280, y = 8, labels = 'Reptilia', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(7.5, 7.5, 21.5, 21.5), col = alpha(phylo$pal[4], 0.1), border = NA)
		text(x = 2280, y = 14.5, labels = 'Aves', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(21.5, 21.5, 30.5, 30.5), col = alpha(phylo$pal[5], 0.1), border = NA)
		text(x = 2280, y = 26, labels = 'Mammalia', adj = 1, cex = 0.8)
#	polygon(x = c(1520,2300, 2300, 1520), y = c(30.5, 30.5, 31.5, 31.5), col = alpha(phylo$pal[6], 0.1), border = NA)
#		text(x = 2280, y = 30.5, labels = 'Diplopoda', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(30.5, 30.5, 32.5, 32.5), col = alpha(phylo$pal[6], 0.1), border = NA)
		text(x = 2280, y = 32.5, labels = 'Malacostraca', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(32.5, 32.5, 51.5, 51.5), col = alpha(phylo$pal[7], 0.1), border = NA)
		text(x = 2280, y = 42, labels = 'Insecta', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(51.5, 51.5, 55.5, 55.5), col = alpha(phylo$pal[8], 0.1), border = NA)
		text(x = 2280, y = 53.5, labels = 'Arachnida', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(55.5, 55.5, 57.5, 57.5), col = alpha(phylo$pal[9], 0.1), border = NA)
		text(x = 2280, y = 56.5, labels = 'Clitellata', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(57.5, 57.5, 59.5, 59.5), col = alpha(phylo$pal[10], 0.1), border = NA)
		text(x = 2280, y = 58.5, labels = 'Lycopodiopsida', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(59.5, 59.5, 64.5, 64.5), col = alpha(phylo$pal[12], 0.1), border = NA)
		text(x = 2280, y = 62, labels = 'Polypodiopsida', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(64.5, 64.5, 68.5, 68.5), col = alpha(phylo$pal[13], 0.1), border = NA)
		text(x = 2280, y = 68.5, labels = 'Magnoliopsida', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(68.5, 68.5, 77.5, 77.5), col = alpha(phylo$pal[14], 0.1), border = NA)
		text(x = 2280, y = 73, labels = 'Liliopsida', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(77.5, 77.5, 102.5, 102.5), col = alpha(phylo$pal[13], 0.1), border = NA)
		text(x = 2280, y = 92, labels = 'Magnoliopsida', adj = 1, cex = 0.8)
	polygon(x = c(1520,2300, 2300, 1520), y = c(102.5, 102.5, 103.5, 103.5), col = alpha(phylo$pal[15], 0.1), border = NA)
		text(x = 2280, y = 103, labels = 'Gnetopsida', adj = 1, cex = 0.8)


	#Polygon to hide the axis I can't suppress...
	par(fig=c(0, 1, 0, 1), oma=c(0, 2.5, 2, 0), mar=c(0, 0, 0, 0), new=TRUE)
	plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
	legend('bottomright', legend = c('a','somethingreallyreallylongandboring'), text.col = 'white', bg = 'white', box.col = 'white')

	par(oma = c(1,0,0,0), mar=c(3, 0, 0, 0), new = TRUE)
	plot.phylo(tr, show.tip.label = FALSE, edge.color = 'grey',
		edge.width = 2, font = 0.8, x.lim = c(0,2300))

	axis(1, at = seq(1530, 2130, 200), labels = c(parse(text = '10^0'), parse(text = '10^1'), parse(text = '10^2'), parse(text = '10^3')))
	mtext('Number of taxa              ', 1, outer = TRUE, line = -0.8, at = c(0, 1), adj = 1, cex = 1.3)

	dev.off()
}



##FIGURE S2 - ANALYSIS METHODS
png('figures/figS2.png', width = 800, height = 800)
{
	par(mfrow = c(2,2))
	par(mai = c(0.8, 0.7, 0.1, 0.3))
	par(oma = c(2, 2, 0, 1.5))
	pal <- paletteer_d("ggthemes::excel_Aspect")

	#Panel A - single species data
	{
		sppCount = 1
		rangeCount = 1
		fullrange <- 0:100
		predx <- -300:400
		ranges <- data.frame(lower = c(15,45, 7), upper = c(82,98, 63))
		as <- c(5.44757, 1.3198549  , -2.173864, 2, 1.5)
		bs <- c(-0.1389627, -0.0468397, 0.0153978, -0.045, -0.073)
		params <- data.frame(a = as, b = bs)
		#Open empty plot window
		plot(0,0, col = 'white', xlim = c(0,100), ylim = c(-0.02,1.02),
			xlab = "", ylab = "", cex.axis = 2)

		#Generate random data
		set.seed(71842)
		presence <- rnorm(10, mean = 28, sd = 8)
			points(presence, rep(1, 10), pch = 19, lwd = 2, cex = 3, col = alpha(pal[3], 0.5))
		absence <- runif(30, ranges[rangeCount, 1], ranges[rangeCount, 2])
			points(absence, rep(0, 30), pch = 21, lwd = 2, cex = 3, col = alpha(pal[3], 0.5))
		#Fit model
		x.data <- c(presence, absence)
		mod <- glm(c(rep(1, 10), rep(0, 30)) ~ x.data, family = 'binomial', weights = c(rep(1, 10), rep(0.3, 30)))
		a <- coef(mod)[1]
		b <- coef(mod)[2]
		params[rangeCount, ] <- coef(mod)
		#Add fitted lines
		x <- seq(ranges[rangeCount,1], ranges[rangeCount,2])
		y <- exp(a + b*fullrange)/(1 + exp(a + b*fullrange))
		lines(fullrange, y, lwd = 5, col = alpha(pal[3], 1))
			
		text(x = 97, y = 0.99, labels = c('(A)'), cex = 2)
		}
			
			
	#Panel B - taxa among datasets
	{
		#Open empty plot window
		plot(0,0, col = 'white', xlim = c(0,100), ylim = c(-0.02,1.02),
			xlab = "", ylab = "", cex.axis = 2)

		#Predicted values for derivatives
		mod_exprA <- expression((exp(a + b*fullrange)/(1 + exp(a + b*fullrange))))
		x_p <- D(mod_exprA, 'fullrange')
		x_pp <- D(x_p, 'fullrange')

		#Generate data for use in combined model
		mod1 <- glm(c(rep(1, 10), rep(0, 30)) ~ x.data, family = 'binomial')
		mod1 <- list(a = coef(mod1)[1], b = coef(mod1)[2])
		obsmod1 <- with(mod1, eval(mod_exprA))
		lines(fullrange, obsmod1, lwd = 2, col = alpha(pal[3], 0.4))

		set.seed(7931)
		presence2 <- rnorm(10, mean = 35, sd = 5)
		absence2 <- runif(30, 20, 100)
		x.data2 <- c(presence2, absence2)
		mod2 <- glm(c(rep(1, 10), rep(0, 30)) ~ x.data2, family = 'binomial')
		mod2 <- list(a = coef(mod2)[1], b = coef(mod2)[2])
		obsmod2 <- with(mod2, eval(mod_exprA))
		lines(fullrange, obsmod2, lwd = 2, col = alpha(pal[3], 0.4))

		set.seed(9270)
		presence3 <- rnorm(10, mean = 27, sd = 3)
		absence3 <- runif(30, 0, 100)
		x.data3 <- c(presence3, absence3)
		mod3 <- glm(c(rep(1, 10), rep(0, 30)) ~ x.data3, family = 'binomial')
		mod3 <- list(a = coef(mod3)[1], b = coef(mod3)[2])
		obsmod3 <- with(mod3, eval(mod_exprA))
		lines(fullrange, obsmod3, lwd = 2, col = alpha(pal[3], 0.4))
	
		dat <- data.frame(
			obs = rep(c(rep(1, 10), rep(0, 30)), 3),
			agb = c(x.data, x.data2, x.data3),
			survey = c(rep('dat1', 40), rep('dat2', 40), rep('dat3', 40))
			)
	
		#Fit and plot mixed model
		model <- glmer(obs ~ agb + (1 + agb|factor(dat$survey)), family = binomial, data = dat)
		mod <- list(a = fixef(model)[1], b = fixef(model)[2])
		obsmod <- with(mod, eval(mod_exprA))
		lines(fullrange, obsmod, lwd = 5, col = alpha(pal[3], 1))

		text(x = 97, y = 0.99, labels = c('(B)'), cex = 2)
		}

	#Panel C - turning points
	{
		sppCount = 1
		rangeCount = 1
		#Find turning points
		cf <- list(a = params[rangeCount,1], b = params[rangeCount,2])
		mod_expr <- expression((exp(a + b*fullrange)/(1 + exp(a + b*fullrange))))
		#Calculate expressions for derivatives of the model
		x_p <- D(mod_expr, 'fullrange')
		x_pp <- D(x_p, 'fullrange')
		#Predicted values for derivatives
		obs.vals <- with(cf, eval(mod_expr))
		d1.vals <- with(cf, eval(x_p))
		d2.vals <- with(cf, eval(x_pp))
		#Make list for use in root.finder function
		fits <- data.frame (agb = fullrange, obs = obs.vals, d1 = d1.vals, d2 = d2.vals)
		fitvals <- list(fits = fits)

		#Find turning points
		turns <- root.finder(fitted_vals = fitvals)
		mod_expr2 <- expression((exp(a + b*predx)/(1 + exp(a + b*predx))))
		predx <- turns$tps
		tps.d2 <- with(cf, eval(mod_expr2))
		predx <- turns$maxrate
		tps.d1 <- with(cf, eval(mod_expr2))


		#Open empty plot window
		plot(0,0, col = 'white', xlim = c(0,100), ylim = c(-0.02,1.02),
			xlab = "", ylab = "", cex.axis = 2)
		#Add fitted lines
		x <- seq(ranges[rangeCount,1], ranges[rangeCount,2])
		y <- exp(a + b*fullrange)/(1 + exp(a + b*fullrange))
		lines(fullrange, y, lwd = 5, col = alpha(pal[3], 1))
		#Plot 1st derivatives
		d1 <- (0.9 - 0.1) * (d1.vals - min(d1.vals)) / (max(d1.vals) - min(d1.vals)) + 0.1		#Rescale to fit in plot window
		lines(d1, col = alpha(pal[3], 1))
		#Plot 2nd derivatives
		d2 <- (0.8 - 0.2) * (d2.vals - min(d2.vals)) / (max(d2.vals) - min(d2.vals)) + 0.2		#Rescale to fit in plot window
		lines(d2, lty = 3, col = alpha(pal[3], 1))
		#Add turning points
		points(turns$tps, tps.d2, pch = 17, lwd = 2, cex = 5, col = alpha(pal[3], 0.8))
		points(turns$maxrate, tps.d1, pch = 25, lwd = 2, cex = 5, col = alpha(pal[3], 0.8))
			
		text(x = turns$tps+5, y = tps.d2+0.05, labels = c('turning point'), cex = 1.5, adj = 0)
		text(x = turns$maxrate+5, y = tps.d1-0.05, labels = c('peak change point'), cex = 1.5, adj = 0)
		text(x = 97, y = 0.99, labels = c('(C)'), cex = 2)
		}
		
		
	#Panel D - turning point rules
	{
		#Open empty plot window
		plot(0,0, col = 'white', xlim = c(0,100), ylim = c(-0.02,1.02),
			xlab = "", ylab = "", cex.axis = 2)
		#Add data range
		predx <- 0:100
		
		#Add fitted lines
		
		cf1 <- cf
		cf2 <- as.list(params[2,])
		cf3 <- as.list(params[3,])
			names(cf1) <- names(cf2) <- names(cf3) <- c("a","b")
		#Predicted values for derivatives
		mod_exprA <- expression((exp(a + b*predx)/(1 + exp(a + b*predx))))
		x_p <- D(mod_exprA, 'predx')
		x_pp <- D(x_p, 'predx')
		obs2 <- with(cf2, eval(mod_exprA))
		obs3 <- with(cf3, eval(mod_exprA))
		obs2.d2 <- with(cf2, eval(x_pp))
		obs3.d2 <- with(cf3, eval(x_pp))
		
		#Add lines to plot
		lines(fullrange, y, lwd = 5, col = alpha(pal[3], 1))
		lines(fullrange, obs2, lwd = 5, col = alpha(pal[3], 1))
		lines(fullrange, obs3, lwd = 5, col = alpha(pal[3], 1))
	
		#Turning points
		mod_expr2 <- expression((exp(a + b*predx)/(1 + exp(a + b*predx))))
		predx <- 0
		tps2 <- with(cf2, eval(mod_expr2))
		predx <- 100
		tps3 <- with(cf3, eval(mod_expr2))

		#Add turning points
		points(turns$tps, tps.d2, pch = 17, lwd = 2, cex = 5, col = alpha(pal[3], 0.8))
		points(0, tps2, pch = 24, lwd = 2, cex = 5, col = alpha(pal[3], 0.8))
		points(100, tps3, pch = 24, lwd = 2, cex = 5, col = alpha(pal[3], 0.8))

		#Add text to indicate special cases
		text(x = 6, y = tps2+0.05, labels = c('1'), cex = 1.5)	#Outside 0:100 range so gets set to zero with rangeweight applied
		text(x = 94, y = tps3+0.05, labels = c('2'), cex = 1.5)	#Outside data bounds but within 0:100 range so gets rangeweight applied

		text(x = 97, y = 0.99, labels = c('(D)'), cex = 2)
		}
		
	
	mtext('Biomass reduction (%)', 1, line = 0, outer = TRUE, cex = 2)
	mtext('Occurrence probability', 2, line = 0, outer = TRUE, cex = 2)

	dev.off()
}








###########################
##NOT USED###########



##FIGURE - VIOLIN PLOTS OF TURNING POINTS PER TAXA
{	

	par(oma = c(0, 0, 0, 0))
	par(mai = c(0.8, 1.0, 0.3, 0.3))
	pal <- paletteer_d("ggthemr::solarized")
	turn_points$dataset$TaxonType <- factor(turn_points$dataset$TaxonType, levels = c('plant', 'invertebrate', 
		'mammal', 'bird', 'reptile', 'amphibian', 'fish'))
		
	#Plot graph
	vioplot(turn.point ~ factor(TaxonType), data = turn_points$dataset,
		ylim = c(0,100), ylab = '', xlab = '', yaxt = 'n', cex.axis = 2,
		col = alpha(pal[1:(length(taxcats))], 0.3), linecol = pal[1:(col.index - 1)],
		horizontal = TRUE)
	axis(1, at = seq(0, 100, 20), cex.axis = 2)
	stripchart(turn.point ~ factor(TaxonType), data = turn_points$dataset, 
		method = 'jitter', jitter = 0.25, vertical = FALSE, add = TRUE, pch = 19, cex = 1.5,
		col = alpha(pal[1:(length(taxcats))], 0.4))


	#Add image labels
	plant<-readPNG("data/leaf.png")
	grid.raster(plant, x=0.07, y=0.2, width=0.1)

	insect<-readPNG("data/invert.png")
	grid.raster(insect, x=0.07, y=0.32, width=0.1)

	mammal<-readPNG("data/mammal.png")
	grid.raster(mammal, x=0.07, y=0.44, width=0.1)

	bird<-readPNG("data/aves.png")
	grid.raster(bird, x=0.07, y=0.54, width=0.1)

	reptile<-readPNG("data/reptile.png")
	grid.raster(reptile, x=0.07, y=0.65, width=0.1)

	amphib<-readPNG("data/frog.png")
	grid.raster(amphib, x=0.07, y=0.77, width=0.1)
	
	fish<-readPNG("data/fish.png")
	grid.raster(fish, x=0.07, y=0.88, width=0.1)

	mtext('Biomass reduction (%)', 1, line = 3, cex = 2)
	
	dev.off()
	}



##FIGURE SXX - FUNCTIONAL RESPONSE PATTERNS


	par(mai = c(0.8, 0.1, 0.1, 0.1))
	par(oma = c(0, 0, 0, 0))
	#pal <- paletteer_d("ggthemes::excel_Aspect")
	pal <- paletteer_d("Redmonder::qPBI")
	xvals <- 0:100
	
	#Get and filter data
	funcs <- arrange.funcplot(func_groups = func.groups, func_points = func_points)
	cats <- unique(funcs$category)
		cats <- cats[!is.na(cats)]

	#Panel A: Body mass
	plot(0,0, xlim = c(-0,100), ylim = c(0, 1), col = 'white', 
		xlab = '', ylab = '',  yaxt = 'n', xaxt = 'n', cex.lab = 2.5, cex.axis = 2)
	axis(1, at = c(0, 50, 100), cex.axis = 2)

	subdata <- funcs[funcs$category == 'Body mass' & funcs$pval < 0.05, ]	#Significant relationshps within that category
	
	for(i in 1:nrow(subdata)){
		fitvals <- func_fits$obs[ , which(names(func_fits$obs) == subdata$taxon[i])[1]]
		lines(x = xvals, y = func_fits$obs[,i], col = subdata$col[i], lwd = 2)
		}



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
	
