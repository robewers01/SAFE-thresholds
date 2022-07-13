
require(lme4)
require(pastecs)


#Function to align site x species matrix with lidar data
#Returns single dataframe for analysis
align.data <- function(comm.out, lidar, predictor, min.obs = 5){
	#comm.out = comm.out section from taxa.clean function
	#lidar = dataframe containing lidar estimates at point locations
	#predictor = dependent variable(s) for analyses
	#min.obs = minimum number of occurrences needed to model the taxon

	tax <- comm.out[  4:ncol(comm.out)]
	comm.out[, 4:ncol(comm.out)] <- tax
	
	#Add predictor data
	preds.include <- match(predictor, names(lidar))						#Select just the lidar variables used in analysis
	preds <- as.data.frame(lidar[match(comm.out$site, lidar$site), preds.include])		#Extract lidar data for all sites in the analysis
		names(preds) <- predictor
		rownames(preds) <- 1:nrow(preds)
	comm.out <- cbind(preds, comm.out)									#Combine into single dataframe
	names(comm.out)[1:length(preds.include)] <- predictor
	comm.out <- comm.out[complete.cases(comm.out[ , 1:length(predictor)]), ]	#Remove cases with incomplete lidar data

	#Exclude taxa with < min.obs occurrences
	tax.start <- which(names(comm.out) == 'date') + 1		#Index of column contining data for first taxon in dataset
	if(tax.start == ncol(comm.out)){
		if(sum(comm.out[ , tax.start], na.rm = TRUE) >= min.obs){
			taxa.include <- tax.start
			}else{
				taxa.include <- NULL
				}
		}else{
		taxa.include <- which(colSums(comm.out[, tax.start:ncol(comm.out)], na.rm = TRUE) >= min.obs) + (tax.start-1)
		}
	if(length(taxa.include) >= 1){
		comm.out <- comm.out[ , c(1:(tax.start-1), taxa.include)]
		}else{
			comm.out <- comm.out[ , 1:(tax.start-1)]
			}

	return(comm.out)
	}


###################################################################################
###################################################################################


#Function to extract and align data from multiple surveys per taxon
align.data.multi <- function(taxon, taxa_survey, full_data, lidar, predictor, min.obs = 1){
	#taxon = taxon name to extract
	#taxa_survey = output from call to taxon.dataset
	#full_data = full dataset containing all data to model
	#lidar = dataframe containing lidar estimates at point locations
	#predictor = dependent variable(s) for analyses
	#min.obs = minimum number of occurrences needed to model the taxon


	surveys <- names(taxa_survey)[which(!is.na(taxa_survey[i, ]))]	#Which surveys contain that taxon
	data.out <- NULL
	for(k in 1:length(surveys)){		#For each of those surveys
		target <- full_data[[match(surveys[k], names(full_data))]]
		comm <- align.data(target$comm.out, lidar, predictor, min.obs = min.obs)	#Add lidar data
		taxonind <- match(rownames(taxa_survey)[i], names(comm))
		if(!is.na(taxonind)){
			commsub <- comm[ , c(1:(length(predictor)+ 3),taxonind)]
			commsub$survey <- surveys[k]
			data.out <- rbind(data.out, commsub)
			}
		}
	 return(data.out)
	 }

###################################################################################
###################################################################################

		
#Function to fit univariate glm models
#Returns best model and name of the best univariate predictor
fit.glm <- function(comm, predictor, taxon_ind = ncol(target)){
	#comm = output from align.data
	#taxon_ind = column index for taxon within comm to be analysed
	#predictor = list of predictor variable(s) to be analysed

	#Fit GLMs and extract AIC values
	mod0 <- try(glm(comm[,taxon_ind] ~ 1, family = binomial, ), FALSE)	#Fit null model
	aic <- numeric()
	for(scale in 1:length(predictor)){						#For each predictor...
		assign('pred', comm[, which(names(comm) == predictor[scale])])
		model <- try(glm(comm[,taxon_ind] ~ pred, family = binomial, ), FALSE)
		assign(paste('mod', scale, sep = ""), model)
		aic <- c(aic, AIC(model))
		rm(model)
		}
	#Select best model
	if(sum(diff(aic)) == 0){		#Where all models give same AIC or there's only a single model
		best <- mod1		#Select model with smallest spatial scale
		bestpred <- predictor[1]
		}else{
			best <- get(paste('mod', (which(aic == min(aic))), sep = ""))
			bestpred <- predictor[which(aic == min(aic))]
			}
			
	pval <- anova(best, mod0, test = 'Chi')$Pr[2]		#Likelihood ratio test
	
	return(list(bestmod = best, bestpred = bestpred, pval = pval))
	}


###################################################################################
###################################################################################


#Function to extract summary statistics from glm
glm.extract <- function(glms_output, taxa_data, comm = comm, taxon, coefs = NA){
	#glms_output = output from call to fit.glm
	#taxa_data = dataset list object containing data used to fit model
	#comm = comm data used in analysis; output from align.data
	#taxon = taxon name
	#coefs = previously extracted stats from glm.extract

	if(!is.data.frame(coefs)){
		coefs <- data.frame(matrix(ncol=17, nrow=0))		#Empty dataframe to store summary statistics
		colnames(coefs) <- c('dataset', 'data.group', 'year', 'taxon', 'modtype', 'r2', 'num.obs', 
			'num.occs', 'best.pred', 'min.x', 'max.x', 'intercept', 'slope', 'interceptSE', 'slopeSE', 'pval', 'method')
		}
	
	dataset <- taxa_data$data.name
	data.group <- paste(taxa_data$data.group, taxa_data$method, sep = '')
	data.group <- gsub(' ', '', data.group)
	year <- taxa_data$sample.year
	taxon <- taxon
	modtype <- 'glm'
	r2 <- with(summary(glms_output$bestmod), 1 - deviance/null.deviance)
	num.obs <- nrow(comm)
	num.occs <- sum(comm[, match(taxon, names(comm))], na.rm = TRUE)
	best.pred <- glms_output$bestpred
	min.x <- min(comm[, match(best.pred, names(comm))], na.rm = TRUE)
	max.x <- max(comm[, match(best.pred, names(comm))], na.rm = TRUE)
	intercept <- summary(glms_output$bestmod)$coef[1]
	slope <- summary(glms_output$bestmod)$coef[2]
	interceptSE <- summary(glms_output$bestmod)$coef[3]
	slopeSE <- summary(glms_output$bestmod)$coef[4]
	pval <- glms_output$pval
	method <- taxa_data$method
	
	stats <- c(dataset, data.group, year, taxon, modtype, r2, num.obs, num.occs, best.pred,
		min.x, max.x, intercept, slope, interceptSE, slopeSE, pval, method)
	coefs[nrow(coefs)+1, ] <- stats
	
	return(coefs)
	}

###################################################################################
###################################################################################


#Function to fit model-averaged binomial GLM and estimate turning points in predicted values
fit.models <- function(taxa_data, lidar, predictor, min.obs = 5){
	#taxa_data = output from taxa.clean function
	#lidar = dataframe containing lidar estimates at point locations
	#predictor = dependent variable; one of 'agb', 'tch'
	#min.obs = minimum number of presences required to analyse taxa


	coefs <- NA
	models <- list()
	for(d in 1:length(taxa_data)){		#For each dataset...
		print(paste('OPENING DATASET', d, 'OUT OF', length(taxa_data)))
		analysis_data <- taxa_data[[d]]
		comm <- align.data(comm.out = analysis_data$comm.out, lidar = lidar, 
			predictor = predictor, min.obs = min.obs)

		if(nrow(comm) >= 5 && ncol(comm) > length(predictor) + 3){		#Only continue if there are at least one taxon AND five sites in the dataset
			#Taxon-by-taxon analysis
			for(species in (4 + length(predictor)): ncol(comm)){		#For each taxon in the dataset....
				target <- comm[ , c(1:(length(predictor) + 3), species)]		#Extract taxon

				print(paste("Fitting models to taxa ", species-(3 + length(predictor)), " out of ", ncol(comm)-(3 + length(predictor)), " in the dataset"))

				glms <- fit.glm(comm = target, taxon_ind = ncol(target), predictor = predictor)		#Fit series of univariate glms and return best model
				coefs <- glm.extract(glms_output = glms, taxa_data = analysis_data, comm = comm,
					taxon = names(comm)[species], coefs = coefs)
				models[[nrow(coefs)]] <- glms								#Record the model itself
				rm(glms)
				}
			}
		}
	
	#Remove any models that weren't fully fitted or tested
	if(!is.na(coefs)[1]){
		toretain <- !is.na(coefs$pval)
		coefs <- coefs[toretain, ]
		models <- models[toretain]
		}
	names(models) <- paste(coefs$dataset, coefs$taxon, sep = "__")
		
	return(list(coefs = coefs, models = models))
	}

###################################################################################
###################################################################################

#Function to calculate predicted values and derivatives
fitted.vals <- function(model, predx = seq((-50), 150, 0.1), mod_coef){
	#model = fitted model for which values should be predicted
	#predx = x values for which to return fitted values
	#mod_coef = summary model information (passed for follow-up functions)

	cf <- as.list(coef(model$bestmod))
	names(cf) <- c("a","b")
	mod_expr <- expression((exp(a + b*predx)/(1 + exp(a + b*predx))))

	#Calculate expressions for derivatives of the model
	x_p <- D(mod_expr, 'predx')
	x_pp <- D(x_p, 'predx')
	x_ppp <- D(x_pp, 'predx')

	#Predicted values for derivatives
	obs <- with(cf, eval(mod_expr))
	d1.vals <- with(cf, eval(x_p))
	d2.vals <- with(cf, eval(x_pp))
	d3.vals <- with(cf, eval(x_ppp))
	
	#Combine for output
	out <- data.frame(agb = predx, obs = obs, d1 = d1.vals, d2 = d2.vals)
	
	return(list(fits = out, coefs = mod_coef))
	}

###################################################################################
###################################################################################


	
#Function to find first turning point
root.finder <- function(fitted_vals){
	#fitted_vals = output from fitted.vals
	
	if (!require(pastecs)) install.packages("pastecs") && require(pastecs)   ## Check if required packages are installed
	 
	x <- fitted_vals$fits$agb
	y <- round(fitted_vals$fits$d2,8)
	keep <- which(y != Inf & !is.na(y))
	x <- x[keep]
	y <- y[keep]
	
	past <- NA
	try(past <- turnpoints(y)$tp, silent = TRUE)
	try(peak <- turnpoints(y)$firstispeak, silent = TRUE)
	if(!is.na(past[1])){
		first <- x[past]
		tps = min(first)
		}else{
			tps <- NA
			peak <- NA
			}
		
	return(list(tps = tps, peak = peak, coefs = fitted_vals$coefs))		
	}

###################################################################################
###################################################################################


#Function to extract turning point data
extract.turns <- function(turns_estimate){
	#turns_estimate = output from root.finder

	obs.x.range <- as.numeric(c(turns_estimate$coefs[, 10:11]))
	turnings <- turns_estimate$tps
	
	#If no turning point was detected:
	if(is.na(turnings) & as.numeric(turns_estimate$coefs$pval) < 0.05){			
		turnings <- 0	#Set turning points to be the min (because fitted model is significant)
		}

	#Check for turning points that fall outside range of valid AGB values
	if(as.numeric(turns_estimate$coefs$slope) > 0){
		if(!is.na(turns_estimate$peak) & turns_estimate$peak == FALSE) turnings <- 0		#Turning point returned has missed initial acceleration point
		}
	if(as.numeric(turns_estimate$coefs$slope) < 0){ 
		if(!is.na(turns_estimate$peak) & turns_estimate$peak == TRUE) turnings <- 0		#Turning point returned has missed initial deceleration point
		}

	#Constrain to fit within valid range of agb values
	if(turnings > 100) turnings <- 100
	if(turnings < 0) turnings <- 0

	#Generate weights for turning point estimates
	if(turnings >= min(obs.x.range) & turnings <= max(obs.x.range)) turn.weight <- 1
	if(turnings < min(obs.x.range)) turn.weight <- abs(1/(turnings-min(obs.x.range)))
	if(turnings > max(obs.x.range)) turn.weight <- abs(1/(turnings-max(obs.x.range)))

	return(list(turn.point = turnings, turn.weight = turn.weight))
	}
	
###################################################################################
###################################################################################



#Function to estimate turning points for fitted models
turns <- function(fitted_mod){
	#fitted_mod = output from fit.models

	turnpoints <- matrix(NA, nrow = length(fitted_mod$models), ncol = 2)
	for(i in 1:length(fitted_mod$models)){			#For each fitted model....
		print(paste('getting turning points for model', i, 'of', length(fitted_mod$models)))
		target = fitted_mod$models[[i]]
		#Check if fitted model was significant or not
		if(as.numeric(fitted_mod$coefs$pval[i]) < 0.05){	#If it was a significant model
			#Find turning points
			fits <- fitted.vals(target, mod_coef = fitted_mod$coefs[i , ])		#Estimate fitted values and derivatives
			turnings <- root.finder(fitted_vals = fits)
			turnsdat <- unlist(extract.turns(turns_estimate = turnings))
			turnpoints[i ,] <- turnsdat
			}
		}
	out <- data.frame(cbind(fitted_mod$coefs, turnpoints))
	names(out)[(ncol(out)-1):ncol(out)] <- c('turn.point', 'range.weight')
	out[ , c(3,6:8,10:16,18:19)] <- sapply(out[, c(3,6:8,10:16,18:19)], as.numeric)
	
	return(out)
	}

###################################################################################
###################################################################################



#Function to summarise datasets
summarise.data <- function(full.data){
	#full.data = raw data
	
	datasets <- names(full.data)
	datagroup <- NULL
	for(i in 1:length(datasets)){
		datagroup[i] <- full.data[[i]]$data.group
		}
	groups <- unique(datagroup)
	surveys <- methods <- periods <- taxa <- numeric()
	for(k in 1:length(groups)){
		targets <- which(datagroup == groups[k])
		surveys[k] <- length(targets)				#Number of surveys
		surv_meth <- surv_year <- surv_taxa <- NULL
		for(m in 1:length(targets)){
			#Extract methods
			surv_meth <- c(surv_meth, full.data[[targets[m]]]$method)
			#Extract year
			surv_year <- c(surv_year, full.data[[targets[m]]]$sample.year)
			#Extract taxa
			surv_taxa <- c(surv_taxa, names(full.data[[targets[m]]]$comm.out)[4:length(names(full.data[[targets[m]]]$comm.out))])
			}
		taxa[k] <- length(unique(surv_taxa))
		methods[k] <- length(unique(surv_meth))
		periods[k] <- length(unique(surv_year))
		}
	
	out <- data.frame(groups, taxa, methods, periods, surveys)
	out <- out[order(groups), ]
	
	return(out)
	
	}


##########################################################################################################
##########################################################################################################


#Function to calculate number of taxa within datasets
summarise.taxa <- function(full_data, min.obs = 5){
	#full_data = dataset to summarise
	#min.obs = minimum number of occurrences needed to model a taxon

	rawtaxa <- mintaxa <- NULL
	for(i in 1:length(full_data)){
		target <- full_data[[i]]$comm.out
		rawtaxa <- c(rawtaxa, names(target)[4:ncol(target)])
		if(ncol(target) > 4){
			taxind <- which(colSums(target[ , 4:ncol(target)]) >= min.obs)
			}else{
				taxind <- which(sum(target[ , 4]) >= min.obs)
				}
		mintaxa <- c(mintaxa, names(target)[taxind + 3])
		}
	all.taxa <- unique(rawtaxa)
	modelled.taxa <- unique(mintaxa)
	num.all.taxa <- length(all.taxa)
	num.modelled.taxa <- length(modelled.taxa)

	out <- list(all.taxa = all.taxa, modelled.taxa = modelled.taxa, 
		num.all.taxa = num.all.taxa, num.modelled.taxa = num.modelled.taxa)
	return(out)
	}

##########################################################################################################
##########################################################################################################


#Function to assign taxa to broad taxonomic categories and level of taxonomic precision
assign.taxon <- function(dataset){
		#dataset = dataframe containing 'taxon' field
	
	if (!require(safedata)) install.packages("safedata") && require(safedata)   ## Check if required packages are installed

	#Import full taxa list from safedata
	taxa <- jsonlite::fromJSON('https://www.safeproject.net/api/taxa')
	taxa$worksheet_name <- taxa$taxon_name
	taxa <- safedata:::taxon_index_to_taxon_table(taxa)

	#Extract taxa that match directly
	matched.ind <- match(unique(dataset$taxon), taxa$taxon_name)
	matched.taxa <- taxa[matched.ind[!is.na(matched.ind)] , ]
	
	#Align aggregated taxa that don't match directly
	if(length(grep('Agg', dataset$taxon)) > 0){
		aggs <- unique(dataset$taxon)[is.na(matched.ind)]		#Vector of aggregated taxa names
		agg.levels <- c('genus')
		for(i in 1:length(agg.levels)){
			targets <- aggs[grep(paste('Agg_', agg.levels[i], sep = ''), aggs)]
			targ.names <- gsub(pattern = paste('Agg_', agg.levels[i], '_', sep = ''), replacement = '', x = targets)
			targ.ind <- match(targ.names, taxa[, match(agg.levels[i], names(taxa))])
			targ.sub <- taxa[targ.ind, ]
			targ.sub$taxon_name <- targets
			targ.sub$taxon_level <- agg.levels[i]
			matched.taxa <- rbind(matched.taxa, targ.sub)
			}
		}
		
	#Assign taxonomic categories
	TaxonType <- rep(NA, nrow(matched.taxa))
	TaxonType[matched.taxa$kingdom == 'Plantae'] <- 'plant'
	TaxonType[matched.taxa$class == 'Mammalia'] <- 'mammal'
	TaxonType[matched.taxa$class == 'Aves'] <- 'bird'
	TaxonType[matched.taxa$class == 'Reptilia'] <- 'reptile'
	TaxonType[matched.taxa$class == 'Amphibia'] <- 'amphibian'
	TaxonType[matched.taxa$class == 'Actinopterygii'] <- 'fish'
	TaxonType[matched.taxa$phylum == 'Platyhelminthes' | matched.taxa$phylum == 'Nematoda' | matched.taxa$phylum == 'Annelida' | matched.taxa$phylum == 'Arthropoda' | matched.taxa$phylum == 'Mollusca' | matched.taxa$phylum == 'Nematomorpha'] <- 'invertebrate'
	TaxonType[matched.taxa$kingdom == 'Chromista'] <- 'microbe'
	TaxonType[matched.taxa$kingdom == 'Fungi'] <- 'fungi'
	matched.taxa$TaxonType <- TaxonType
	
	#Assign level of taxonomic resolution
	TaxonAgg <- matched.taxa$taxon_level
	TaxonAgg[grep('Agg_order', matched.taxa$taxon_name)] <- 'order'
	#Tidy up classifications
	TaxonAgg[TaxonAgg == 'morphospecies'] <- 'species'
	TaxonAgg[TaxonAgg == 'class'] <- NA
	TaxonAgg[TaxonAgg == 'subclass'] <- NA
	TaxonAgg[TaxonAgg == 'infraorder'] <- NA
	TaxonAgg[TaxonAgg == 'phylum'] <- NA
	matched.taxa$TaxonAgg <- TaxonAgg

	
	#Append information to dataset
	dataset$TaxonType <- as.factor(matched.taxa$TaxonType[match(dataset$taxon, matched.taxa$taxon_name)])
	dataset$TaxonAgg <- as.factor(matched.taxa$TaxonAgg[match(dataset$taxon, matched.taxa$taxon_name)])
	
	out <- list(dataset = dataset, matched.taxa = matched.taxa)
	
	return(out)
	}

###################################################################################
###################################################################################


#Function to return cumulative distribution
cdf <- function(data, x_range = c(0:100)){
	#data = vector of data to calculate cdf for
	#x_range = range of x-values to calculate cdf over
	
	sums <- numeric()
	for(i in min(x_range):max(x_range)){
		sums <- c(sums, sum(data <= i, na.rm = TRUE))
		}
	props <- sums/ length(data)
	
	out <- data.frame(x = x_range, cdf = sums, prop = props)
	return(out)
	}

###################################################################################
###################################################################################


	#Function to create matrix of observations and rates of change for all taxa
fitted.matrix <- function(models, predx = c(0:100)){
	#models = list of models to make predictions for
	#predx = set of x-values to make predictions for

	rates <- obs <- data.frame(matrix(NA, nrow = 101, ncol = length(models)))
	for(i in 1:length(models)){
		target <- models[[i]]
		if(target$pval < 0.05){
			fitted <- fitted.vals(target, predx = predx, mod_coef = NA)$fits
			rates[ , i] <- fitted$d1
			obs[ , i] <- fitted$obs
			}
		}
	names(rates) <- names(obs) <- names(models)
	rownames(rates) <- rownames(obs) <- predx
	return(list(obs = obs, rates = rates))
	}
	
###################################################################################
###################################################################################


#Function to create taxon x dataset list
taxon.dataset <- function(full_data){
	#full_data = output from taxa.clean

	taxa.dataset <- NA
	for(i in 1:length(full_data)){
		comm.out <- full_data[[i]]$comm.out
		if(!is.data.frame(taxa.dataset)){	
			if(ncol(comm.out) > 4){
				taxa.dataset <- as.data.frame(matrix(colSums(comm.out[, 4:ncol(comm.out)], na.rm = TRUE),
						nrow = length(4:ncol(comm.out)), ncol = 1))
				}else{
					taxa.dataset <- as.data.frame(matrix(sum(comm.out[, 4], na.rm = TRUE),
						nrow = length(4:ncol(comm.out)), ncol = 1))
					}
			rownames(taxa.dataset) <- names(comm.out)[4:ncol(comm.out)]
			colnames(taxa.dataset) <- full_data[[i]]$data.name
			}else{
				taxa.dataset <- cbind(taxa.dataset, matrix(NA, nrow(taxa.dataset), ncol = 1))
					names(taxa.dataset)[ncol(taxa.dataset)] <- full_data[[i]]$data.name
				matched.taxa <- match(names(comm.out)[4:ncol(comm.out)], rownames(taxa.dataset))
					matched.taxa2 <- matched.taxa[!is.na(matched.taxa)]
				if(length(matched.taxa2) > 0){	
					if(ncol(comm.out) > 4){
						taxa.dataset[matched.taxa2, ncol(taxa.dataset)] <- 
							colSums(comm.out[, 4:ncol(comm.out)], na.rm = TRUE)[which(!is.na(matched.taxa))]
						}else{
							taxa.dataset[matched.taxa2, ncol(taxa.dataset)] <- sum(comm.out[, 4], na.rm = TRUE)
							}
					}
				new.taxa <- names(comm.out)[4:ncol(comm.out)][is.na(matched.taxa)]
					new.taxa <- new.taxa[!is.na(new.taxa)]
				if(length(new.taxa) > 0){	
					new.dataset <- as.data.frame(matrix(NA, nrow = length(new.taxa), ncol = ncol(taxa.dataset)))
						names(new.dataset) <- names(taxa.dataset)
						rownames(new.dataset) <- new.taxa
				if(ncol(comm.out) > 4){	
					new.dataset[ , ncol(new.dataset)] <- colSums(comm.out[, 4:ncol(comm.out)], na.rm = TRUE)[match(new.taxa, names(comm.out))-3]
					}else{
						new.dataset[ , ncol(new.dataset)] <- sum(comm.out[, 4], na.rm = TRUE)[match(new.taxa, names(comm.out))-3]
						}
					taxa.dataset <- rbind(taxa.dataset, new.dataset)
					}
				}
		}
	return(taxa.dataset)
	}

###################################################################################
###################################################################################









	