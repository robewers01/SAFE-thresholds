
require(lme4)
require(pastecs)
require(strucchange)


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
		taxonind <- match(taxon, names(comm))
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

		
#Function to extract and align data from multiple surveys per functional group
align.data.func <- function(func, func_survey, full_data, lidar, predictor, min.obs = 1){
	#func = functional group name to extract
	#func_survey = output from call to function.dataset
	#full_data = full dataset containing all data to model
	#lidar = dataframe containing lidar estimates at point locations
	#predictor = dependent variable(s) for analyses
	#min.obs = minimum number of occurrences needed to model the taxon

	surveys <- names(func_survey)	#Which surveys contain that functional group
	data.out <- NULL
	for(k in 1:length(surveys)){		#For each of those surveys
		target <- full_data[[match(surveys[k], names(full_data))]]
		comm <- align.data(target$comm.out, lidar, predictor, min.obs = min.obs)	#Add lidar data
		
		if(nrow(comm) > 0){
			#Aggregate data across all relevant taxa
			taxonind <- match(rownames(func_survey), names(comm))
			taxonind <- taxonind[!is.na(taxonind)]
			if(length(taxonind) == 1){
				commB <- data.frame(comm[, taxonind])
				}else{
					commB <- data.frame(rowSums(comm[, taxonind]))
					}
			commB[commB > 0] <- 1
			names(commB) <- func
			commC <- cbind(comm[ , 1:(length(predictor)+ 3)], commB)
			rownames(commC) <- rownames(comm)
			commC$survey <- surveys[k]
			data.out <- rbind(data.out, commC)
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

#Function to fit univariate generalised mixed effect models
#Returns best model and name of the best univariate predictor
fit.glmer <- function(comm, predictor){
	#comm = output from align.data.multi
	#predictor = list of predictor variable(s) to be analysed

	if (!require(lme4)) install.packages("lme4") && require(lme4)   ## Check if required packages are installed

	responseInd <- length(predictor) + 4		#Index to taxon data
	
	#Fit GLMMs and extract AIC values
	aic <- numeric()
	for(scale in 1:length(predictor)){						#For each predictor...
		assign('pred', comm[, which(names(comm) == predictor[scale])])
		model <- try(glmer(comm[,responseInd] ~ pred + (1 + pred|factor(comm$survey)), family = binomial), FALSE)
		if(class(model) != "try-error"){
			assign(paste('mod', scale, sep = ""), model)
			aic <- c(aic, AIC(model))
			rm(model)
			}else{
			  aic <- c(aic, NA)
			}
		}
		mod0 <- try(glmer(comm[,responseInd] ~ 1 + (1 + 1|factor(comm$survey)), family = binomial), FALSE)	#Fit null model
	#Select best model
	if(sum(is.na(aic)) != length(aic)){	#Check at least one model was fitted
		if(sum(diff(aic), na.rm = TRUE) == 0){		#Where all models give same AIC or there's only a single model
			best <- get(paste('mod', (which(!is.na(aic)== TRUE)[1]), sep = ""))	#Select model with smallest spatial scale
			bestpred <- predictor[which(!is.na(aic)== TRUE)[1]]
			}else{
				best <- get(paste('mod', (which(aic == min(aic, na.rm = TRUE))), sep = ""))
				bestpred <- predictor[which(aic == min(aic, na.rm = TRUE))]
				}
				
		pval <- anova(best, mod0, test = 'Chi')$Pr[2]		#Likelihood ratio test
		}else{
			best <- bestpred <- pval <- NA
			}
	
	return(list(bestmod = best, bestpred = bestpred, pval = pval))
	}

###################################################################################
###################################################################################


#Function to extract summary statistics from glm
glm.extract <- function(glms_output, comm = comm, taxon, coefs = NA){
	#glms_output = output from call to fit.glm
	#comm = comm data used in analysis; output from align.data
	#taxon = taxon name
	#coefs = previously extracted stats from glm.extract

	if(!is.data.frame(coefs)){
		coefs <- data.frame(matrix(ncol=11, nrow=0))		#Empty dataframe to store summary statistics
		colnames(coefs) <- c('taxon', 'modtype', 'r2', 'num.obs', 
			'num.occs', 'best.pred', 'min.x', 'max.x', 'intercept', 'slope', 'pval')
		}
	
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
	pval <- glms_output$pval
	
	stats <- c(taxon, modtype, r2, num.obs, num.occs, best.pred,
		min.x, max.x, intercept, slope, pval)
	coefs[nrow(coefs)+1, ] <- stats
	
	return(coefs)
	}

###################################################################################
###################################################################################


#Function to extract summary statistics from glmer
glmer.extract <- function(glmer_output, comm = comm, taxon, coefs = NA){
	#glmer_output = output from call to fit.glmer
	#comm = comm data used in analysis; output from align.data
	#taxon = taxon name
	#coefs = previously extracted stats from glm.extract

	if (!require(MuMIn)) install.packages("MuMIn") && require(MuMIn)   ## Check if required packages are installed

	if(!is.data.frame(coefs)){
		coefs <- data.frame(matrix(ncol=11, nrow=0))		#Empty dataframe to store summary statistics
		colnames(coefs) <- c('taxon', 'modtype', 'r2', 'num.obs', 
			'num.occs', 'best.pred', 'min.x', 'max.x', 'intercept', 'slope', 'pval')
		}
	
	taxon <- taxon
	modtype <- 'glmer'
	modR2 <- try(r.squaredGLMM(glmer_output$bestmod)[1,1], silent = TRUE)
	if(class(modR2) != "try-error"){
		r2 <- modR2
		}else{
			r2 <- NA
			}
	num.obs <- nrow(comm)
	num.occs <- sum(comm[, match(taxon, names(comm))], na.rm = TRUE)
	best.pred <- glmer_output$bestpred
	min.x <- min(comm[, match(best.pred, names(comm))], na.rm = TRUE)
	max.x <- max(comm[, match(best.pred, names(comm))], na.rm = TRUE)
	intercept <- fixef(glmer_output$bestmod)[1]
	slope <- fixef(glmer_output$bestmod)[2]
	pval <- glmer_output$pval
	
	stats <- c(taxon, modtype, r2, num.obs, num.occs, best.pred,
		min.x, max.x, intercept, slope, pval)
	coefs[nrow(coefs)+1, ] <- stats
	
	return(coefs)
	}

###################################################################################
###################################################################################

#Function to fit binomial GLM or GLMER as appropriate
fit.models <- function(full_data, func_data = NA, lidar, predictor, min.observs = 5){
	#full_data = raw data to be analysed
	#func_data = functional traits data; only use when analysing functional groups rather than individual taxa
	#lidar = dataframe containing lidar estimates at point locations
	#predictor = vector of dependent variables
	#min.observs = minimum number of presences required to analyse taxa

	taxa_survey <- taxon.dataset(full_data)
	to.iterate <- 1:nrow(taxa_survey)
	#Extract list of taxa belonging to each functional group
	if(is.data.frame(func_data)){
		taxa.lists <- taxaXfunc(func_groups = func.groups)$taxon.lists
		to.iterate <- names(taxa.lists)
		}

#	models <- list()
	coefs <- NA

	for(i in 1:length(to.iterate)){
		#Get aligned dataset
		if(!is.data.frame(func_data)){
			print(paste('fitting models to taxon', i, 'of', length(to.iterate), ':', rownames(taxa_survey)[i]))
			comm_data <- align.data.multi(taxon = rownames(taxa_survey)[i], taxa_survey = taxa_survey, full_data = full_data,
				lidar = lidar, predictor = predictor,
				min.obs = 1)
			taxon <- rownames(taxa_survey)[i]
			}else{
				print(paste('fitting models to taxon', i, 'of', length(to.iterate), ':', to.iterate[i]))
				longlist <- match(taxa.lists[[i]], rownames(taxa_survey))
				if(length(longlist) == 1 & is.na(longlist[1])){
					comm_data <- NULL
					}else{
						shortlist <- longlist[!is.na(longlist)]
						taxdat <- taxa_survey[shortlist, ]		#All taxa in here belong to that functional group
						shorttaxa <- rownames(taxdat)
						#Extract list of relevant datasets
						colinds <- which(colSums(taxdat, na.rm = TRUE) > 0)
						if(length(colinds) > 1){
							funcdat <- taxdat[ , colinds]
							}else{
								funcdat <- as.data.frame(taxdat[ , colinds])
								names(funcdat) <- names(colinds)
								}
						comm_data = align.data.func(func = names(taxa.lists)[i], func_survey = funcdat, full_data = thresh.data,
							lidar = lidar.data, predictor = predictor, min.obs = 1)
						taxon <- to.iterate[i]
						}
				}
		
			
		if(!is.null(comm_data)){
			
			if(sum(comm_data[, (length(predictor)+4)], na.rm = TRUE) >= min.observs){		#Check there are enough presences to model
				#Decide if mixed effect model is needed
				mixed <- FALSE
				if(length(unique(comm_data$survey)) > 1)	mixed <- TRUE		#If there's >1 survey, try a GLMER
				
				#Fit models
				if(mixed){
					comm = comm_data
					model <- fit.glmer(comm = comm, predictor = predictor)
					if(is.na(model$bestmod)) mixed <- FALSE		#If the model failed, then try a GLM
					}
				if(!mixed){
					model <- fit.glm(comm = comm_data, predictor, taxon_ind = length(predictor) + 4)
					coefs <- glm.extract(glms_output = model, comm = comm_data,
						taxon = taxon, coefs = coefs)
					}
				if(mixed){
					coefs <- glmer.extract(glmer_output = model, comm = comm,
						taxon = taxon, coefs = coefs)					
					}
#				models[[i]] <- model
#				names(models)[[i]] <- rownames(taxa_survey)[i]	
				rm(mixed, model, comm_data, comm)
				}else{
					coefs[i ,] <- c(taxon, NA, NA, nrow(comm_data), sum(comm_data[, (length(predictor)+4)]), rep(NA, 6))
					}
			}else{
				coefs[i ,] <- c(taxon, rep(NA, 10))
				}
			
		}
	
#	return(list(coefs = coefs, models = models))
	return(coefs)
	}


###################################################################################
###################################################################################

#Function to calculate predicted values and derivatives
fitted.vals <- function(predx = seq(-1000, 1000, 1), mod_coef){
	#predx = x values for which to return fitted values
	#mod_coef = summary model information (passed for follow-up functions)

	cf <- list(a = as.numeric(mod_coef$intercept), b = as.numeric(mod_coef$slope))
	mod_expr <- expression((exp(a + b*predx)/(1 + exp(a + b*predx))))

	#Calculate expressions for derivatives of the model
	x_p <- D(mod_expr, 'predx')	#First derivative
	x_pp <- D(x_p, 'predx')		#Second derivative

	#Predicted values
	obs <- with(cf, eval(mod_expr))
	d1.vals <- with(cf, eval(x_p))
	d2.vals <- with(cf, eval(x_pp))
	
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
	yD1 <- round(fitted_vals$fits$d1,8)
	keep <- which(y != Inf & !is.na(y))
	x <- x[keep]
	y <- y[keep]
	yD1 <- abs(yD1[keep])
	
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
	try(rate <- turnpoints(yD1)$tp, silent = TRUE)
	if(!is.na(rate[1])){
		maxrate <- x[rate]
		#Constrain to fit within valid range of agb values
		if(maxrate > 100) maxrate <- 100
		if(maxrate < 0) maxrate <- 0
		}else{
			maxrate <- NA
			if(yD1[1] > yD1[length(yD1)]){
				maxrate <- 0
				}
			if(yD1[1] < yD1[length(yD1)]){
				maxrate <- 100
				}
			}

		
	return(list(tps = tps, peak = peak, maxrate = maxrate, coefs = fitted_vals$coefs))		
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

	#Strip out taxa that weren't modelled
	fitted_mod <- fitted_mod[!is.na(fitted_mod$modtype) , ]
	
	turnpoints <- matrix(NA, nrow = nrow(fitted_mod), ncol = 3)
	for(i in 1:nrow(fitted_mod)){			#For each fitted model....
		print(paste('getting turning points for model', i, 'of', nrow(fitted_mod)))
		#Check if fitted model was significant or not
		if(as.numeric(fitted_mod$pval[i]) < 0.05){	#If it was a significant model
			#Find turning points
			fits <- fitted.vals(mod_coef = fitted_mod[i , ])		#Estimate fitted values and derivatives
			turnings <- root.finder(fitted_vals = fits)
			turnsdat <- unlist(extract.turns(turns_estimate = turnings))
			turnpoints[i ,] <- c(turnsdat, turnings$maxrate)
			}
		}
	out <- data.frame(cbind(fitted_mod, turnpoints))
	names(out)[(ncol(out)-2):ncol(out)] <- c('turn.point', 'range.weight', 'maxrate')
	out[ , c(3:5, 7:13)] <- sapply(out[, c(3:5, 7:13)], as.numeric)
	
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


#Function to assign taxa to broad taxonomic categories and level of taxonomic precision
assign.taxon <- function(dataset, taxon_table = NA){
		#dataset = dataframe containing 'taxon' field
	
	if (!require(safedata)) install.packages("safedata") && require(safedata)   ## Check if required packages are installed

	
	if(is.data.frame(taxon_table)){
		taxa <- taxon_table		
		}else{
			#Import full taxa list from safedata
			taxa <- jsonlite::fromJSON('https://www.safeproject.net/api/taxa')
			taxa$worksheet_name <- taxa$taxon_name
			taxa <- safedata:::taxon_index_to_taxon_table(taxa)
			}

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
	dataset$Order <- as.factor(matched.taxa$order[match(dataset$taxon, matched.taxa$taxon_name)])
	
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
	#models = output from call to fit.models
	#predx = set of x-values to make predictions for

	rates <- obs <- data.frame(matrix(NA, nrow = length(predx), ncol = nrow(models)))
	for(i in 1:nrow(models)){
		target <- models[i, ]
		fits <- fitted.vals(predx = predx, mod_coef = target)$fits
		rates[ , i] <- fits$d1
		obs[ , i] <- fits$obs
		}
	names(rates) <- names(obs) <- models$taxon
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


#Function to extract lists of taxa belonging to all functional groups
taxaXfunc <- function(func_groups){
	#func_groups = functional traits data

	func_groups <- func_groups[!is.na(func_groups$taxon_name), ]
	taxon.names <- unique(func_groups$taxon_name)
	func_groups <- as.data.frame(func_groups[match(taxon.names, func_groups$taxon_name), ])
	
	#Generate vector of functional groups
	groups <- extract.func(func_groups)
	groups.unique <- groups$group
	
	#Extract taxon list
	out.list <- list()
	for(i in 1:length(groups.unique)){		#For each grouping factor
		target.group <- groups.unique[i]
		target.func <- as.data.frame(func_groups[ , match(target.group, names(func_groups))])
		rownames(target.func) <- rownames(func_groups)
	#		#Remove duplicates and rename rows by taxon_name
	#		taxon.names <- unique(func_groups$taxon_name)
	#		target.func <- as.data.frame(target.func[match(taxon.names, func_groups$taxon_name), ])
	#		rownames(target.func) <- taxon.names
		target.vals <- unique(c(unique(unlist(apply(target.func, 2, unique)))))
		target.vals <- target.vals[!is.na(target.vals)]
		out.names <- paste(target.group, '_', target.vals, sep = '')	#Unique functional groups to form
		for(j in 1:length(out.names)){	#For each of those functional groups...
			#Generate list of taxa
			taxanames <- NULL
			for(k in 1:ncol(target.func)){	#For each column...
				taxanames <- c(taxanames, rownames(target.func)[which(target.func[, k] == target.vals[j])])
				}
			out.list[[length(out.list)+1]] <- unique(taxanames)
			names(out.list)[length(out.list)] <- out.names[j]
			}
		}
	
	return(list(taxon.lists = out.list, functional.groups = groups))
	}

###################################################################################
###################################################################################


#Function to pull out relevant functional groups for analysis
extract.func <- function(func_groups){
	#func_groups = functional traits data

	groups <- names(func_groups)
	#Exclude numeric traits (these have all been categorised and are captured in other fiels)
	num.col <- NULL
	for(k in 1:ncol(func_groups)){
		num.col <- c(num.col, is.numeric(func_groups[, k]))
		}
	groups <- groups[!num.col]
	#Exclude unnecessary taxonomic information
	groups <- groups[-c(1:12,15)]
	
	to.retain <- NULL	#Vector to store indices to desired groups
	func.categ <- NULL	#Vector to record broad categories of groups
	keeps <- grep('_', groups)
	exclude.alltaxa <- (1:length(groups))[-keeps]	#functional groups relating to all taxa
	#Body mass by taxon
	to.retain <- c(to.retain, grep('BodyMass', groups))
		#Exclusions
		exclusions <- exclude.alltaxa
		exclusions <- c(exclusions, which(groups == 'BodyMass_Categorised'))
	to.retain <- to.retain[is.na(match(to.retain, exclusions))]
	to.retain <- c(to.retain, to.retain)
	func.categ <- rep('Body mass', length(to.retain))
	#Movement by taxon
	movement <- grep('Movement', groups)
		#Exclusions
		exclusions <- exclude.alltaxa
	movement <- movement[is.na(match(movement, exclusions))]
	to.retain <- c(to.retain, movement)
	func.categ <- c(func.categ, rep('Movement', length(movement)))
	#Strata types by taxon
	strata <- grep('Stra', groups)
		#Exclusions
		exclusions <- exclude.alltaxa
		#Exclude all taxa types combined
		#Exclude life history variation
		exclusions <- c(exclusions, grep('StraLifeHist', groups))
		#Exclude generalism index
		exclusions <- c(exclusions, grep('StrataGeneralism', groups))
	strata <- strata[is.na(match(strata, exclusions))]
	to.retain <- c(to.retain, strata)
	func.categ <- c(func.categ, rep('Habitat strata', length(strata)))
	#Strata life history variation
	slh <- grep('StraLifeHist', groups)
	to.retain <- c(to.retain, slh[-grep('_', groups[slh])])
	func.categ <- c(func.categ, rep('Habitat strata', length(slh[-grep('_', groups[slh])])))
	#Strata generalism
	sg <- grep('StrataGeneralism', groups)
	to.retain <- c(to.retain, sg[1])
	func.categ <- c(func.categ, rep('Habitat strata', length(sg[1])))
	#IUCNthreat and RedList
	iucn <- grep('IUCN', groups)
	to.retain <- c(to.retain, iucn[-grep('_', groups[iucn])])
	func.categ <- c(func.categ, rep('Red List status', length(iucn[-grep('_', groups[iucn])])))
	#Development
	dev <- grep('Development', groups)
	to.retain <- c(to.retain, dev[-grep('_', groups[dev])])
	func.categ <- c(func.categ, rep('Development', length(dev[-grep('_', groups[dev])])))
	#Physiology
	phys <- grep('Physiology', groups)
	to.retain <- c(to.retain, phys[-grep('_', groups[phys])])
	func.categ <- c(func.categ, rep('Physiology', length(phys[-grep('_', groups[phys])])))
	#Sociality
	social <- grep('Sociality', groups)
	to.retain <- c(to.retain, social[-grep('_', groups[social])])
	func.categ <- c(func.categ, rep('Sociality', length(social[-grep('_', groups[social])])))
	#Trophic
	trophic <- grep('Troph', groups)
	to.retain <- c(to.retain, trophic[-grep('_', groups[trophic])])
	func.categ <- c(func.categ, rep('Trophic', length(trophic[-grep('_', groups[trophic])])))
	#Trophic generalism
	tg <- grep('TrophicGen', groups)
	to.retain <- c(to.retain, tg[1])
	func.categ <- c(func.categ, rep('Trophic', length(tg[1])))
	#Diet
	diet <- grep('Diet', groups)
	to.retain <- c(to.retain, diet[-grep('_', groups[diet])])
	func.categ <- c(func.categ, rep('Diet', length(diet[-grep('_', groups[diet])])))
	#Diet generalism
	dg <- grep('DietGen', groups)
	to.retain <- c(to.retain, dg[1])
	func.categ <- c(func.categ, rep('Trophic', length(dg[1])))
	#Plants
	plant <- grep('Plant', groups)
	to.retain <- c(to.retain, plant[-grep('_plant', groups[plant])])
	func.categ <- c(func.categ, rep('Plant', length(plant[-grep('_plant', groups[plant])])))
	
	out <- data.frame(group = groups[to.retain], category = func.categ)
	#Remove generic vertebrate group
	out <- out[-grep('vertebrate', out$group), ]
	#Add taxon column
	out$taxon <- 'all taxa'
	out$taxon[grep('bird', out$group)] <- 'bird'
	out$taxon[grep('invert', out$group)] <- 'invertebrate'
	out$taxon[grep('amphibian', out$group)] <- 'amphibian'
	out$taxon[grep('reptile', out$group)] <- 'reptile'
	out$taxon[grep('mammal', out$group)] <- 'mammal'
	out$taxon[grep('plant', out$group)] <- 'plant'
	out$taxon[grep('fish', out$group)] <- 'fish'
	
	return(out)
	}

###################################################################################
###################################################################################


#	func_points <- readRDS('results/func_points.rds')
#	func_grousp = func.groups

#Function to categorise and name functional groups
rename.funcs <- function(func_groups, func_points){
		#func_groups = functional groups data
		#func_points = output from turns applied to functional groups
		
	newfuncs <- func_points
	funcs <- extract.func(func_groups)
	
	
	aligned <- sapply(funcs$group, grep, func_points$taxon, value = TRUE)
	aligned.data <- data.frame(matrix(NA, nrow = 0, ncol = 2))
	for(i in 1:length(aligned)){
		full <- aligned[[i]]
		short <- rep(names(aligned)[i], length(full))
		new <- cbind(full, short)
		aligned.data <- rbind(aligned.data, new)
		}
	aligned.data$TaxType <- funcs$taxon[match(aligned.data$short, funcs$group)]
	aligned.data$TaxType[grep('Plant', aligned.data$full)] <- 'plant'
	aligned.data$category <- funcs$category[match(aligned.data$short, funcs$group)]

	newfuncs$shortname <- aligned.data$short[match(newfuncs$taxon, aligned.data$full)]
	newfuncs$TaxType <- aligned.data$TaxType[match(newfuncs$taxon, aligned.data$full)]
	newfuncs$category <- aligned.data$category[match(newfuncs$taxon, aligned.data$full)]

	qualifier <- gsub('.*_Categorised_', '' , newfuncs$taxon)
	qualifier <- gsub('Diet.*_', '' , qualifier)
	qualifier <- gsub('Troph.*_', '' , qualifier)
	qualifier <- gsub('Sociality.*_', '' , qualifier)
	qualifier <- gsub('IUCN.*_', '' , qualifier)
	qualifier <- gsub('Stra.*_', '' , qualifier)
	qualifier <- gsub('Develo.*_', '' , qualifier)
	qualifier <- gsub('Physio.*_', '' , qualifier)
	qualifier <- gsub('Movem.*_', '' , qualifier)
	
	qualifier[grep('WoodDens', newfuncs$shortname)] <- paste('wood density', qualifier[grep('WoodDens', newfuncs$shortname)], sep = '-')
	qualifier[grep('Ruderal', newfuncs$shortname)] <- paste('ruderal', qualifier[grep('Ruder', newfuncs$shortname)], sep = '-')
	qualifier[grep('Stress', newfuncs$shortname)] <- paste('stress', qualifier[grep('Stress', newfuncs$shortname)], sep = '-')
	qualifier[grep('Comp', newfuncs$shortname)] <- paste('competitor', qualifier[grep('Comp', newfuncs$shortname)], sep = '-')
	qualifier[grep('PhotoCap', newfuncs$shortname)] <- paste('photosynthesis', qualifier[grep('Photo', newfuncs$shortname)], sep = '-')

	qualifier[grep('TrophicGen', newfuncs$shortname)] <- paste('generalism', qualifier[grep('TrophicGen', newfuncs$shortname)], sep = '-')
	qualifier[grep('DietGen', newfuncs$shortname)] <- paste('generalism', qualifier[grep('DietGen', newfuncs$shortname)], sep = '-')
	qualifier[grep('StrataGen', newfuncs$shortname)] <- paste('generalism', qualifier[grep('StrataGen', newfuncs$shortname)], sep = '-')
	qualifier[grep('Threatened', qualifier)] <- 'threatened'
	qualifier[grep('Not threatened', qualifier)] <- 'not threatened'
	
	newfuncs$qualifier = qualifier
	newfuncs$category[grep('diet gen', newfuncs$qualifier)] <- 'Diet'

	newfuncs$plotname <- gsub('_.*', '' , newfuncs$shortname)
	newfuncs$plotname <- gsub('Plant', '' , newfuncs$plotname)
	
	newfuncs$fullname <- paste(newfuncs$category, newfuncs$qualifier, newfuncs$TaxType, sep = '-')
	newfuncs$description <- paste(newfuncs$category, newfuncs$qualifier, sep = '-')
	
	return(newfuncs)
	}

###################################################################################
###################################################################################


#Function to find example taxa linked to functional groups
find.egs <- function(Function, function_qual, taxtype = NA, 
	turn_points, func_points, func_groups){
	#Function = functional group; column header in func_groups
	#function_qual = subfunction... varaible nested within Function
	#taxtype = target Taxa Type
	#turn_points = output from turns called on taxon data
	#func_points = output from turns called on functional group data
	#func_groups = functional traits data
	
	#Find slope of functional group
	funcs <- rename.funcs(func_groups = func.groups, func_points = func_points)
	if(is.na(taxtype)){
		fullname <- paste(Function, function_qual, sep = '_')
		}else{
			fullname <- paste(Function, taxtype, function_qual, sep = '_')
			}
	funcslope <- funcs$slope[match(fullname, funcs$taxon)]
	if(funcslope > 0) func.dir = 'positive' 
	if(funcslope < 0) func.dir = 'negative'


	column <- which(names(func_groups) == Function)
	rows <- which(func_groups[, column] == function_qual)
	if(!is.na(taxtype)){
		taxcols <- which(func_groups$TaxonType == taxtype)
		rows <- intersect(taxcols, rows)
		}
	possibles <- func_groups$taxon_name[rows]
	possmatch <- match(possibles, turn_points$taxon)
	possmatch <- possmatch[!is.na(possmatch)]

	poss_coefs <- turn_points[possmatch, ]
	if(func.dir == 'positive'){
		finalset <- poss_coefs$taxon[poss_coefs$slope > 0 & poss_coefs$pval < 0.05]
		}
	if(func.dir == 'negative'){
		finalset <- poss_coefs$taxon[poss_coefs$slope < 0& poss_coefs$pval < 0.05]
		}
	
	return(sort(finalset))
	}

###################################################################################
###################################################################################

#Function to find breakpoints in curvilinear bivariate plots
break.points <- function(x, y, ss = 10){
	#x = vector of x-axis values
	#y = vector of y-axis values
	#ss = what fraction of observations to use when scanning for breakpoints? 
		#This speeds up the fitting, but is only called if length(x) > 1000

	if (!require(strucchange)) install.packages("strucchange") && require(strucchange)   ## Check if required packages are installed

	dat = data.frame(x = x, y = y)
	
	#Subset every ss'th observations 
	if(length(x) > 1000){
		data = dat[c(seq(1, nrow(dat), ss)), ]
		}else{
			data = dat}
	
	#Find breakpoints
	bp <- breakpoints(data$y ~ data$x)
	d <- coef(bp)[, 2]
	#Record breakpoints
	bps<-data$x[bp$breakpoints]
	bps.y<-data$y[bp$breakpoints]
	
	#Classify into acceleration vs deceleration
	rate <- rep(NA, length(bps))
	rate[which(diff(d) > 0)] <- 'accelerate'
	rate[which(diff(d) < 0)] <- 'decelerate'

	out <- data.frame(bps = bps, bps.y = bps.y, type = rate)
	return(out)
	
	}
	
###################################################################################
###################################################################################


#Function to calculate and add breakpoints to figure
plot.breaks <- function(x, y, add_plot = TRUE, decel.col, ...){
	#x = x values
	#y = y values
	#add_plot = binary; add breakpoints to open plot or not?
	#decel.col = colour for deceleration points
	#... = parameters passed to plot function
	
	bp <- break.points(x = x, y = y)
	pt.col <- rep('white', nrow(bp))
		pt.col[bp$type == 'decelerate'] <- alpha(decel.col, 1)
#	for(i in 1:nrow(bp)){
#		arrows(bp$bps[i], bp$bps.y[i], bp$bps[i], -1, col = alpha(pal[3], 0.4), code = 0, lwd = 2, lty = 2)
#		}
	if(add_plot) points(bp$bps, bp$bps.y, bg = alpha(pt.col, 0.8), ...)
	return(bp)
	}

###################################################################################
###################################################################################


#Function to filter and arrange data for plotting functional groups
arrange.funcplot <- function(func_groups, func_points){
	#func_groups = functional traits data
	#func_points = output from turns called on functional group data

	funcs <- rename.funcs(func_groups = func_groups, func_points = func_points)
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
	
	return(funcs)
	}

###################################################################################
###################################################################################


#Function to summarise number of different functional groups in the dataset
func.summary <- function(func_groups){
	#func_groups = functional traits data
	
	diets <- names(func.groups[grep('Diet', names(func.groups))])
	numdiets <- length(diets[-c(grep('_', diets), grep('Gener', diets))])		#Number of diets
	trophic <- names(func.groups[grep('Troph', names(func.groups))])
	numtroph <- length(trophic[-c(grep('_', trophic), grep('Gener', trophic))])		#Number of trophic levels
	move <- names(func.groups)[grep('Movem', names(func.groups))]
	movetypes <- unique(func.groups[, move[-c(grep('_', move))]])
	movetypes <- movetypes[-is.na(movetypes)]
	physiol <- names(func.groups)[grep('Phys', names(func.groups))]
	physioltypes <- unique(func.groups[, physiol[-c(grep('_', physiol))]])
	physioltypes <- physioltypes[-is.na(physioltypes)]
	strata <- names(func.groups)[grep('Stra', names(func.groups))]
	stratatypes <- strata[-c(grep('_', strata))]
	stratatypes <- stratatypes[-c(is.na(stratatypes), grep('Gener', stratatypes))]
	social <- names(func.groups)[grep('Soci', names(func.groups))]
	socialtypes <- unique(func.groups[, social[-c(grep('_', social))]])
	socialtypes <- socialtypes[-is.na(socialtypes)]
	cons <- names(func.groups)[grep('IUCN', names(func.groups))]
	constypes <- unique(func.groups[, cons[-c(grep('_', cons))]])
	constypes <- constypes[-is.na(constypes)]
	
	print(paste('Number of diets:', numdiets))
	print(paste('Number of trophic levels:', numtroph))
	print(paste('Number of movement modes:', length(movetypes)))
	print(paste('Number of physiologies:', length(physioltypes)))
	print(paste('Number of habitat levels:', length(stratatypes)))
	print(paste('Number of sociality syndromes:', length(socialtypes)))
	print(paste('Number of conservation levels:', length(constypes)))

	}

###################################################################################
###################################################################################


##Function to expand functional groups by taxonomic groups
expand.funcs <- function(taxa){
	#taxa = functional groups raw data
	
	taxon_cat <- sort(unique(taxa$TaxonCat))
	taxon_groups <- c('amphibian', 'bird', 'fish', 'mammal', 'reptile')		#Only groups not already represented in taxon_cat
	for(i in 16:ncol(taxa)){
		orig <- taxa[, i]						#Select the trait
		for(k in 1:length(taxon_cat)){			#Iterate through each taxon category
			tax_ind <- which(taxa$TaxonCat == taxon_cat[k])
			new_entry <- rep(NA, nrow(taxa))
			new_entry[tax_ind] <- taxa[tax_ind, i]
			taxa <- cbind(taxa, new_entry)
			names(taxa)[ncol(taxa)] <- paste(names(taxa)[i], taxon_cat[k], sep = '_')
			rm(new_entry)
			}
		for(j in 1:length(taxon_groups)){	#Iterate through each taxon group
			tax_ind <- which(taxa$TaxonType == taxon_groups[j])
			new_entry <- rep(NA, nrow(taxa))
			new_entry[tax_ind] <- taxa[tax_ind, i]
			taxa <- cbind(taxa, new_entry)
			names(taxa)[ncol(taxa)] <- paste(names(taxa)[i], taxon_groups[j], sep = '_')
			rm(new_entry)
			}
		}
	#Delete empty trait columns
	emptycols <- sapply(taxa, function (k) all(is.na(k)))
	taxa <- taxa[, !emptycols]
	
	return(taxa)
	}
	
###################################################################################
###################################################################################


##Function to generate categorical versions of numeric functional traits
func.cats <- function(taxa){
	#taxa = functional traits data matrix

	num_traits <- names(taxa[sapply(taxa, is.numeric)])
	for(i in 1:length(num_traits)){
		col_ind <- match(num_traits[i], names(taxa))
		if(length(grep('Body', num_traits[i])) == 1){
			new_entry <- as.character(cat.data(taxa[, col_ind], num.cats = 3, reverse = FALSE))
			}else{
				new_entry <- as.character(cat.data(taxa[, col_ind], num.cats = 2, reverse = FALSE))
				}
		if(length(unique(new_entry)) >= 2){			#Only keep if multiple categories can be identified
			taxa <- cbind(taxa, new_entry)
			names(taxa)[ncol(taxa)] <- paste(names(taxa)[col_ind], 'Categorised', sep = '_')
			rm(new_entry)
			}
		}
	
	return(taxa)
	}
		
###################################################################################
###################################################################################















#####################################################################
##NOT USED###########################
#######################################
#######################################
#######################################
#######################################



#Function to extract thresholds per functional group
func.thresholds <- function(func_groups, turns_taxa){
	#func_groups = functional traits data
	#turns_taxa = output from call to turns that used individual taxa

	taxa.func <- taxaXfunc(func_groups = func_groups)	#list of taxa per functional group
	turns <- rates <- data.frame(matrix(NA, nrow = 0, ncol = 3))
	for(i in 1:length(taxa.func)){
		print(paste('Working on functional group', i, 'out of', length(taxa.func)))
		#Extract data
		taxa <- taxa.func[[i]]
		taxinds <- match(taxa, turns_taxa$taxon)
		taxinds <- taxinds[!is.na(taxinds)]
		models <- turns_taxa[taxinds, ]
		
		#Construct CDFs
		turnsCDF <- cdf(data = models$turn.point, x_range = c(0:100))
		ratesCDF <- cdf(data = models$maxrate, x_range = c(0:100))
		
		#Breakpoints
		turn_breaks <- try(plot.breaks(x = turnsCDF$x, y = turnsCDF$prop, add_plot = FALSE), silent = TRUE)
		rate_breaks <- try(plot.breaks(x = ratesCDF$x, y = ratesCDF$prop, add_plot = FALSE), silent = TRUE)
		
		if(class(turn_breaks) != "try-error"){			
			turns <- rbind(turns, turn_breaks)
			}
		if(class(rate_breaks) != "try-error"){			
			rates <- rbind(rates, rate_breaks)
			}
		rm(turn_breaks, rate_breaks)
		}
	
	return(list(turn.breaks = turns, rate.breaks = rates))
	}

###################################################################################
###################################################################################




