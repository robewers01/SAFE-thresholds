
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


#Function to generate single, united dataset
full.data <- function(full_data, lidar, predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000')){
	#full_data = full dataset containing all data to model
	#lidar = dataframe containing lidar estimates at point locations
	#predictor = dependent variable(s) for analyses

	combined_data <- NULL
	for(i in 1:length(names(full_data))){
		print(paste('Adding dataset', i, 'out of', length(full_data), sep = ' '))
		targ_comm <- full_data[[i]]$comm.out
		targ_data <- align.data(comm.out = targ_comm, lidar, predictor = predictor, min.obs = 5)
		if(nrow(targ_data) != 0 & ncol(targ_data) > 8){
			restr_data <- restruct.aligned(aligned = targ_data, data.name = full_data[[i]]$data.name, sample.year = full_data[[i]]$sample.year)
			combined_data <- rbind(combined_data, restr_data)
			}
		}
		
	return(combined_data)
	}

###################################################################################
###################################################################################


#Function to convert aligned site x species matrix to vertical table structure
restruct.aligned <- function(aligned, data.name, sample.year){
	#aligned = output from call to align.data
	#data.name = name of the study
	#sample.year = year the data were collected
	
	new_data <- NULL
	for(i in 9:ncol(aligned)){
		dat_extract <- aligned[, c(1:8, i), ]
		dat_extract$taxon <- names(aligned)[i]
		names(dat_extract)[9] <- 'presence'
		dat_extract$year <- sample.year
		dat_extract$studyID <- data.name
		new_data <- rbind(new_data, dat_extract)
		}
	
	return(new_data)
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
				if((length(longlist) == 1 & is.na(longlist[1])) | sum(is.na(longlist)) == length(longlist)){
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
	
	#Remove any taxa with no taxonomic information
	todelete <- grep('Unknown', coefs$taxon)
	if(length(todelete) > 0){
		coefs <- coefs[-todelete, ]
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
	TaxonType[matched.taxa$kingdom == 'Plantae'] <- 'Plant'
	TaxonType[matched.taxa$class == 'Mammalia'] <- 'Mammal'
	TaxonType[matched.taxa$class == 'Aves'] <- 'Bird'
	TaxonType[matched.taxa$class == 'Reptilia'] <- 'Reptile'
	TaxonType[matched.taxa$class == 'Amphibia'] <- 'Amphibian'
	TaxonType[matched.taxa$class == 'Actinopterygii'] <- 'Fish'
	TaxonType[matched.taxa$phylum == 'Platyhelminthes' | matched.taxa$phylum == 'Nematoda' | matched.taxa$phylum == 'Annelida' | matched.taxa$phylum == 'Arthropoda' | matched.taxa$phylum == 'Mollusca' | matched.taxa$phylum == 'Nematomorpha'] <- 'Invertebrate'
		TaxonType[matched.taxa$class == 'Arachnida'] <- 'Arachnid'
		TaxonType[matched.taxa$class == 'Insecta'] <- 'Insect'
		TaxonType[matched.taxa$family == 'Formicidae'] <- 'Ant'
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
	func_groupsEXP <- expand.funcs(taxa = func_groups)
	
	#Extract taxon list
	out.list <- list()
	for(i in 1:length(groups.unique)){		#For each grouping factor
		target.group <- groups.unique[i]
		target.func <- as.data.frame(func_groupsEXP[ , match(target.group, names(func_groupsEXP))])
		rownames(target.func) <- rownames(func_groupsEXP)
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

	func_groups <- expand.funcs(taxa = func_groups)

	groups <- names(func_groups)
	#Exclude numeric traits (these have all been categorised and are captured in other fields)
	num.col <- NULL
	for(k in 1:ncol(func_groups)){
		num.col <- c(num.col, is.numeric(func_groups[, k]))
		}
	groups <- groups[!num.col]
	#Exclude unnecessary taxonomic information
	groups <- groups[-c(1:12)]
	
	to.retain <- NULL	#Vector to store indices to desired groups
	func.categ <- NULL	#Vector to record broad categories of groups
	keeps <- grep('_', groups)
	exclude.alltaxa <- (1:length(groups))[-keeps]	#functional groups relating to individual taxa

	#Body mass by taxon
	to.retain <- c(to.retain, grep('BodyMass', groups))
		#Exclusions
		exclusions <- exclude.alltaxa
		exclusions <- c(exclusions, which(groups == 'BodyMass_Categorised'), grep('Reptile', groups))
	to.retain <- to.retain[is.na(match(to.retain, exclusions))]
	func.categ <- rep('Body mass', length(to.retain))
	#Movement by taxon
	movement <- grep('Movement', groups)
		#Exclusions
		exclusions <- c(grep('Amphibian', groups), grep('Bird', groups), grep('Ant', groups), grep('Arachnid', groups), grep('Fish', groups), grep('Plant', groups))
	movement <- movement[is.na(match(movement, exclusions))]
	to.retain <- c(to.retain, movement)
	func.categ <- c(func.categ, rep('Movement', length(movement)))
	#Strata types by taxon
	strata <- grep('Stra', groups)
		#Exclusions
		#Exclude life history variation
		exclusions <- grep('StraLifeHist', groups)
		#Exclude generalism index
		exclusions <- c(exclusions, grep('StrataGeneralism', groups))
		#Exclude taxa with no meaningful comparisons
		exclusions <- c(exclusions, grep('Ant', groups), grep('Arachnid', groups), grep('Fish', groups), grep('Invertebrate', groups), grep('Plant', groups))
		
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
		#Exclusions
		excl <- c(grep('Amphibian', groups), grep('Bird', groups), grep('Mammal', groups), grep('Arachnid', groups), grep('Ant', groups), grep('Reptile', groups), grep('Fish', groups))
	dev <- dev[is.na(match(dev, excl))]
	to.retain <- c(to.retain, dev)
	func.categ <- c(func.categ, rep('Development', length(dev)))
	#Physiology
	phys <- grep('Physiology', groups)
	to.retain <- c(to.retain, phys[-grep('_', groups[phys])])
	func.categ <- c(func.categ, rep('Physiology', length(phys[-grep('_', groups[phys])])))
	#Sociality
	social <- grep('Sociality', groups)
		#Exclusions
		excl <- c(grep('Ant', groups), grep('Arachnid', groups), grep('Fish', groups), grep('Invertebrate', groups), grep('Plant', groups), grep('Reptile', groups))
	social <- social[is.na(match(social, excl))]
	to.retain <- c(to.retain, social)
	func.categ <- c(func.categ, rep('Sociality', length(social)))
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
	func.categ <- c(func.categ, rep('Diet', length(dg[1])))
	#Plants
	plant <- grep('Plant', groups)
	to.retain <- c(to.retain, plant[-grep('_Plant', groups[plant])])
	func.categ <- c(func.categ, rep('Plant', length(plant[-grep('_Plant', groups[plant])])))
	
	#Create dataframe for output
	out <- data.frame(group = groups[to.retain], category = func.categ)

	#Add taxon column
	out$taxon <- rep('all taxa', nrow(out))
	out$taxon[grep('Bird', out$group)] <- 'Bird'
	out$taxon[grep('Invertebrate', out$group)] <- 'Invertebrate'
	out$taxon[grep('Amphibian', out$group)] <- 'Amphibian'
	out$taxon[grep('Reptile', out$group)] <- 'Reptile'
	out$taxon[grep('Mammal', out$group)] <- 'Mammal'
	out$taxon[grep('Plant', out$group)] <- 'Plant'
	out$taxon[grep('Fish', out$group)] <- 'Fish'
	out$taxon[grep('Arachnid', out$group)] <- 'Arachnid'
	out$taxon[grep('Insect', out$group)] <- 'Insect'
	out$taxon[grep('Ant', out$group)] <- 'Ant'
	
	return(out)
	}

###################################################################################
###################################################################################


#Function to ensure all references to TaxonTypes have capitalised first letter
rename.taxtypes <- function(taxtypes){
	#taxtypes = vector containing strings that embed TaxonTypes
	
	taxtypes <- gsub('amphib', 'Amphib', taxtypes)
	taxtypes <- gsub('plant', 'Plant', taxtypes)
	taxtypes <- gsub('reptile', 'Reptile', taxtypes)
	taxtypes <- gsub('mammal', 'Mammal', taxtypes)
	taxtypes <- gsub('bird', 'Bird', taxtypes)
	taxtypes <- gsub('fish', 'Fish', taxtypes)
	taxtypes <- gsub('invert', 'Invert', taxtypes)
	
	return(taxtypes)
	}
	
###################################################################################
###################################################################################


#Function to categorise and name functional groups
rename.funcs <- function(func_groups, func_points = NA){
		#func_groups = functional groups data
		#func_points = output from turns applied to functional groups
			#NB: Not used anymore...
			
	taxa.lists <- taxaXfunc(func_groups = func_groups)
	fullfuncs <- names(taxa.lists$taxon.lists)
	funcs <- taxa.lists$functional.groups
	newfuncs <- data.frame(taxon = fullfuncs)

	aligned <- sapply(funcs$group, grep, fullfuncs, value = TRUE)
	#Manual exclusions to the alignments
	aligned$Movement <- aligned$Movement[1:3]
	aligned$Development <- aligned$Development[1:2]
	aligned$StraUnd <- aligned$StraUnd[1]
	aligned$StraAqu <- aligned$StraAqu[1]
	aligned$StraAer <- aligned$StraAer[1]
	aligned$StraArb <- aligned$StraArb[1]
	aligned$StraTerr <- aligned$StraTerr[1]
	aligned$StraSub <- aligned$StraSub[1]
	aligned$Sociality <- aligned$Sociality[1:4]
	
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
	funcs <- merge(funcs, func_points, by = 'taxon')
	if(is.na(taxtype)){
		fullname <- paste(Function, function_qual, sep = '_')
		}else{
			fullname <- paste(Function, taxtype, function_qual, sep = '_')
			}
	funcslope <- funcs$slope[match(fullname, funcs$taxon)]
	if(funcslope > 0) func.dir = 'positive' 
	if(funcslope < 0) func.dir = 'negative'

	#Extract set of taxa to check
	d <- taxaXfunc(func_groups)
	possibles <- d$taxon.lists[which(names(d$taxon.lists) == fullname)]

#	column <- which(names(func_groups) == Function)
#	rows <- which(func_groups[, column] == function_qual)
#	if(!is.na(taxtype)){
#		taxcols <- which(func_groups$TaxonType == taxtype)
#		rows <- intersect(taxcols, rows)
#		}
#	possibles <- func_groups$taxon_name[rows]
	possmatch <- match(possibles[[1]], turn_points$taxon)
	possmatch <- possmatch[!is.na(possmatch)]

	poss_coefs <- turn_points[possmatch, ]
	if(func.dir == 'positive'){
		finalset <- poss_coefs$taxon[poss_coefs$slope > 0]# & poss_coefs$pval < 0.1]
		}
	if(func.dir == 'negative'){
		finalset <- poss_coefs$taxon[poss_coefs$slope < 0]# & poss_coefs$pval < 0.1]
		}
	if(func.dir == 'positive'){
		strongfinalset <- poss_coefs$taxon[poss_coefs$slope > 0 & poss_coefs$pval < 0.05]
		}
	if(func.dir == 'negative'){
		strongfinalset <- poss_coefs$taxon[poss_coefs$slope < 0 & poss_coefs$pval < 0.05]
		}
	
	return(list(strong = sort(strongfinalset), weak = sort(finalset)))
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
	out <- out[out$bps > 0, ]
	out <- out[out$type == 'accelerate', ]
	out <- out[!is.na(out$type == '<NA>'), ]
	return(out)
	
	}
	
###################################################################################
###################################################################################


#Function to calculate and add breakpoints to figure
plot.breaks <- function(x, y, add_plot = TRUE, decel.col = 1, ...){
	#x = x values
	#y = y values
	#add_plot = binary; add breakpoints to open plot or not?
	#decel.col = colour for deceleration points
	#... = parameters passed to plot function
	
	bp <- break.points(x = x, y = y)
	bp <- bp[bp$type == 'accelerate', ]		#Remove the deceleration points
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
arrange.funcplot <- function(func_groups, func_points, pal = NA){
	#func_groups = functional traits data
	#func_points = output from turns called on functional group data
	#pal = colour palette for plotting

	funcs <- rename.funcs(func_groups = func_groups, func_points = func_points)
	
	funcs <- merge(funcs, func_points, by = 'taxon')
	
	funcs <- funcs[!is.na(funcs$turn.point) ,]
	funcs <- funcs[!is.na(funcs$category), ]
	funcs <- funcs[!duplicated(funcs), ]
	#Exclude groups with impacts that only start/end at the ends of the gradient
	funcs <- funcs[-which(funcs$turn.point == funcs$maxrate), ]
	funcs <- funcs[-which(funcs$turn.point > funcs$maxrate), ]
	
	#Sort data
	funcs$category <- factor(funcs$category, levels = c('Red List status', 'Habitat strata', 
		'Plant', 'Physiology', 'Development', 'Sociality', 'Movement', 'Diet', 'Trophic',  'Body mass'))
	funcs$qualifier <- factor(funcs$qualifier, levels = c('high', 'medium', 'low', 
		'aerial', 'arboreal', 'understory', 'terrestrial', 'subterranean', 'aquatic', 
		'endotherm', 'ectotherm', 'direct', 'indirect', 
		'eusocial', 'social', 'pair', 'solitary', 
		'winged', 'legged', 'legless', 
		'threatened', 'not threatened', 
		'parasitoid', 'parasite', 'hematophage', 
		'carnivore', 'piscivore', 'vertivore', 'invertivore','bacteriophage', 
		'herbivore', 'frugivore', 'granivore', 'florivore', 'nectarivore', 'palynivore', 'folivore', 
		'phloeophage', 'xylophage', 'rhizophage',
		'algivore','mycophage', 
		'saprophage', 'coprophage', 'necrophage', 'detritivore', 
		'saproxylic', 'soilphage', 
		'producer', 
		'constant', 'variable',	'generalism','generalism-high','generalism-low',
		'wood density', 'photosynthesis', 'competitor', 'stress', 'ruderal', 		
		rev(funcs$qualifier[funcs$TaxType == 'plant'])))
	
	funcs <- funcs[order(funcs$category, funcs$qualifier, funcs$TaxType, funcs$turn.point), ]
	
	#Add colours
	funcs$col <- pal[as.numeric(factor(funcs$category))]
	#Add line types
	funcs$lty <- 1
	funcs$lty[funcs$slope < 0] <- 2
	#symbol types
	types <- c('all taxa', 'invertebrate', 'insect', 'amphibian', 'bird', 'fish', 'mammal', 'plant')
	pchtypes <- c(21:25,8,12, 13)
	pchmatch <- data.frame(cbind(types, pchtypes))
	funcs$pch <- as.numeric(pchmatch$pchtypes[match(funcs$TaxType, pchmatch$types)])
	
	return(list(summary = funcs, pchmatch = pchmatch))
	}

###################################################################################
###################################################################################


#Function to summarise number of different functional groups in the dataset
func.summary <- function(func_groups){
	#func_groups = functional traits data
	
	diets <- names(func.groups[grep('Diet', names(func.groups))])
	numdiets <- length(diets[-c(grep('Gener', diets))])		#Number of diets
	trophic <- names(func.groups[grep('Troph', names(func.groups))])
	numtroph <- length(trophic[-c(grep('Gener', trophic))])		#Number of trophic levels
	move <- names(func.groups)[grep('Movem', names(func.groups))]
	movetypes <- unique(func.groups[, move])
	movetypes <- movetypes[-is.na(movetypes)]
	physiol <- names(func.groups)[grep('Phys', names(func.groups))]
	physioltypes <- unique(func.groups[, physiol])
	physioltypes <- physioltypes[-is.na(physioltypes)]
	strata <- names(func.groups)[grep('Stra', names(func.groups))]
	stratatypes <- strata
	stratatypes <- stratatypes[-c(is.na(stratatypes), grep('Gener', stratatypes))]
	social <- names(func.groups)[grep('Soci', names(func.groups))]
	socialtypes <- unique(func.groups[, social])
	socialtypes <- socialtypes[-is.na(socialtypes)]
	cons <- names(func.groups)[grep('IUCN', names(func.groups))]
	constypes <- unique(func.groups[, cons])
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


#Function to convert first letter of a string to upper case
firstup <- function(x) {
	substr(x, 1, 1) <- toupper(substr(x, 1, 1))
	x
	}

###################################################################################
###################################################################################


##Function to expand functional groups by taxonomic groups

expand.funcs <- function(taxa){
	#taxa = functional groups raw data
	
	taxa$TaxonType <- firstup(taxa$TaxonType)
	
#	taxon_cat <- sort(unique(taxa$TaxonCat))
	taxon_groups <- c('Plant', 'Invertebrate', 'Insect', 'Ant', 'Arachnid', 'Amphibian', 'Bird', 'Fish', 'Mammal', 'Reptile')		#Only groups not already represented in taxon_cat
	for(i in 16:ncol(taxa)){
		orig <- taxa[, i]						#Select the trait
#		for(k in 1:length(taxon_cat)){			#Iterate through each taxon category
#			tax_ind <- which(taxa$TaxonCat == taxon_cat[k])
#			new_entry <- rep(NA, nrow(taxa))
#			new_entry[tax_ind] <- taxa[tax_ind, i]
#			taxa <- cbind(taxa, new_entry)
#			names(taxa)[ncol(taxa)] <- paste(names(taxa)[i], taxon_cat[k], sep = '_')
#			rm(new_entry)
#			}
		for(j in 1:length(taxon_groups)){	#Iterate through each taxon group
			tax_ind <- which(taxa$TaxonType == taxon_groups[j])
			if(length(tax_ind) > 0){
				new_entry <- rep(NA, nrow(taxa))
				new_entry[tax_ind] <- taxa[tax_ind, i]
				taxa <- cbind(taxa, new_entry)
				names(taxa)[ncol(taxa)] <- paste(names(taxa)[i], taxon_groups[j], sep = '_')
				rm(new_entry)
				}
			}
		}
	#Delete empty trait columns
	emptycols <- sapply(taxa, function (k) all(is.na(k)))
	full_ind <- match(names(which(emptycols==FALSE)), names(taxa))
	taxa <- taxa[, full_ind]
	names(taxa) <- sub('.1', '', names(taxa))
	
	#Create categorical versions of numeric traits
	taxa_cat <- func.cats(taxa = taxa)
	
	return(taxa_cat)
	}
	
###################################################################################
###################################################################################


##Function to generate categorical versions of numeric functional traits
func.cats <- function(taxa){
	#taxa = functional traits data matrix

	num_traits <- names(taxa[sapply(taxa, is.numeric)])
	for(i in 1:length(num_traits)){
		col_ind <- match(num_traits[i], names(taxa))
		if(!is.na(col_ind)){	
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
		}
	
	return(taxa)
	}
		
###################################################################################
###################################################################################


#Function to categorise vector of numeric values into quantiles
cat.data <- function(x, num.cats = 3, reverse = FALSE){
	#x = vector of numeric values to be categorised
	#num.cats = number of categories to return
	#reverse = should high categories be given to low values?
	
	if(reverse == TRUE){
		x <- 0 - x
		}
	
	#rank values
	stand <- rank(x, na.last = 'keep')

	#category boundaries
	probs <- (2:num.cats - 1)/num.cats		#quantile boundaries
	bounds <- quantile(x, probs = probs, na.rm = TRUE)
	
	#category names
	if(num.cats == 3 & length(unique(stand)) >= 3){
		name <- c('high', 'medium', 'low') 
		}
	if(num.cats == 2 | length(unique(stand)) == 2){
		name <- c('high', 'low') 
		}
	if(num.cats > 3){
		name <- c('1st', '2nd', '3rd', paste(c(4:num.cats), 'th', sep = ""))
		}
	
	#categorise data
	cats <- rep(NA, length(x))
	cats[x >= bounds[length(bounds)]] <- name[1]
	cats[x < bounds[1]] <- name[length(name)]
	if(num.cats > 2){
		for (i in length(bounds):2)
			cats[x < bounds[i] & x >= bounds[i-1]] <- name[i]
		}
	
	return(as.factor(cats))
	}

###################################################################################
###################################################################################


#Function to return vector of taxa represented in aggregated datasets
taxa.list <- function(dup_data){
	#dup_data = dataset containing raw data
	
	taxa <- character()
	for(i in 1:length(dup_data)){		#For each dataset...
		taxa <- c(taxa, names(dup_data[[i]]$comm.out)[4:ncol(dup_data[[i]]$comm.out)])	#Extract vector of taxa names
		}
		
	return(unique(taxa))
	}

###################################################################################
###################################################################################


#Function to arrange phylogenetic data for plotting
arrange.phylo <- function(timetree, raw_data, taxa_safe, tt_map, coefs, palette_col){
	#timetree = phylogeny generated from timetree.org
	#raw_data = raw data used in analysis
	#taxa_safe = full list of all taxa in SAFE database
	#tt_map = map to connect taxa in data to taxa in timetree
	#coefs = results from fit.models
	#palette_col = colour palette
	
	tr <- timetree
	taxa <- taxa_safe
	map <- tt_map
	taxa_list <- taxa.list(dup_data = raw_data)

	#Extract number of taxa per order
	taxa_mod <- taxa[match(taxa_list, taxa$taxon_name), ]
	orders <- sort(unique(taxa_mod$order))
	tax.order <- taxa_mod[match(orders, taxa_mod$order), ]
		tax.order <- tax.order[, -c(5:13)]			#Delete extraneous columns
	num_taxa <- summary(factor(taxa_mod$order), maxsum = length(orders)+5)
	tax.order$num_taxa <- num_taxa[match(orders, names(num_taxa))]
	
	#Number modelled taxa per order
	coefs2 <- coefs[!is.na(coefs[,2]), ]
	modelled_tax <-sort(unique(coefs2$taxon))
	coef_sub <- taxa_mod[match(modelled_tax, taxa_mod$taxon_name), ]
	numbers <- summary(factor(coef_sub$order), maxsum = length(orders)+5)
	tax.order$num.modelled <- numbers[match(orders, names(numbers))]
		tax.order$num.modelled[is.na(tax.order$num.modelled)] <- 0

	#Update taxonomy to align with timetree.org
	tax.order <- tax.order[-(which(tax.order$order == 'Hypocreales')), ]		#No fungi were included in the analyses
	tax.order <- tax.order[-(which(tax.order$order == 'Collembola')), ]			#No taxon with enough observations to be modelled
	tax.order[tax.order$order == 'Passeriformes', 5:6]	<- tax.order[tax.order$order == 'Passeriformes', 5:6] + tax.order[tax.order$order == 'Hydrornis', 5:6]	#Merge Hydrornis with Passeriformes
		tax.order <- tax.order[-(which(tax.order$order == 'Hydrornis')),]
	tax.order[tax.order$order == 'Blattodea', 5:6]	<- tax.order[tax.order$order == 'Blattodea', 5:6] + tax.order[tax.order$order == 'Isoptera', 5:6]	#Merge Isoptera with Blattodea
		tax.order <- tax.order[-(which(tax.order$order == 'Isoptera')),]
	tax.order$order[tax.order$order == 'Dinosauria'] <- 'Galliformes'		#They were chickens, not dinosaurs...

	#Identify one example species per order that exists on TimeTree.org
	tax.order$tt_fam <- map$V3[match(tax.order$order, map$V1)]
	tax.order$tt_fam <- sub(' ', '_', tax.order$tt_fam)
	
	#Trim to keep just the orders in the dataset
	tr <- drop.tip(tr, which(tr$tip.label =='Polydesmus_complanatus'))			#No millipedes in final dataset
	ords_to_keep <- tr$tip.label[match(tax.order$tt_fam, tr$tip.label)]
	tr$tip.label <- tax.order$order[match(tr$tip.label, tax.order$tt_fam)]

	#Add plotting details
	num1 <- tax.order$num.modelled[match(tr$tip.label, tax.order$order)]
		number <- matrix(num1)
		rownames(number) <- tr$tip.label
		colnames(number) <- 'num_mod'
	num_tax <- tax.order$num_taxa[match(tr$tip.label, tax.order$order)]
		number <- cbind(number, num_tax)
	prop_mod <- num1 / num_tax
		number <- cbind(number, prop_mod)
	prop_not_mod <- 1 - prop_mod
		number <- cbind(number, prop_not_mod)
	class1 <- tax.order$class[match(tr$tip.label, tax.order$order)]
		number <- cbind(number, class1)

	for(i in 1:4){
		number[is.na(number[, i]), i] <- 0
		}

	numbers2 <- data.frame(numTax = log10(as.numeric(number[, 2]) + 1))
	numbers3<- data.frame(propTax = log10(as.numeric(number[, 1]) + 1))	#Number of taxa that were modelled
		numbers3[numbers3 == 0] <- NA
	rownames(numbers2) <- rownames(numbers3) <- tr$tip.label

	classes <- unique(class1)
	col_index <- match(class1, classes)

	return(list(tr = tr, numbers2 = numbers2, numbers3 = numbers3, 
		pal = palette_col, col_index = col_index))
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


#Function to generate taxon resilience summary data
resil.dat <- function(turns, bayes, grouping = 'TaxonType', full_taxa = NA){
	#turns = output from call to 'turns'
	#bayes = output from Bayesian analysis (replicability MS)
	#grouping
	#full_taxa = full taxonomic details of all taxa in analyses

	#Define internal functions
	proportion <- function(x) sum(!is.na(x)) / length(x)
	count.samp <- function(x) sum(!is.na(x))
	prop.generalist <- function(x) sum(x > 0.05) / length(x)

	#Join Bayes probability of correct slope to main results
	turns$bayes.prob <- bayes$prob.slope[match(turns$taxon, bayes$taxon)]
	
	#Add grouping details
	if(grouping != 'TaxonType'){
		coltarg <- match(grouping, names(taxa))
		turns_group <- full_taxa[match(turns$taxon, taxa$taxon_name), coltarg]
		turns <- cbind(turns, turns_group)
		names(turns)[ncol(turns)] <- grouping
		}
	
	#Assign generalists a 'turning point' at 100 % AGB removal
	gen_turns <- turns$turn.point
	gen_turns[turns$pval > 0.05] <- 100

	#Calculate summaries
	group_targ <- match(grouping, names(turns))
	prop_imp <- by(turns$turn.point, factor(turns[, group_targ]), FUN = proportion)	#Proportion taxa with turning points
	mean_turn <- by(gen_turns, factor(turns[, group_targ]), FUN = mean, na.rm = TRUE) 	
	bayes_prob <- by(turns$bayes.prob, factor(turns[, group_targ]), FUN = mean, na.rm = TRUE) 	
	samp_size <- by(turns$bayes.prob, factor(turns[, group_targ]), FUN = count.samp) 	

	#Scale data
	scale_imp <- array(prop_imp)
	scale_turn <- 1-array(mean_turn/100)
	scale_resil <- 1 - (scale_imp * scale_turn)
	new_samp <- array(samp_size)
	
	out <- data.frame(prop_imp = array(prop_imp), mean_turn = array(mean_turn), 
		susceptibility = scale_imp, sensitivity = scale_turn, 
		resilience = scale_resil, 
		bayes_prob = array(bayes_prob), bayes_samp = new_samp)
	rownames(out) <- names(prop_imp)
	
	return(out)
	}

###################################################################################
###################################################################################


#Function to calculate resilience
#resil.calc <- function(turns){
#	#Turns = output from call to 'turns', subsetted to contain just the taxa you want to group and calculate resilience for
#	
#	proportion <- function(x) sum(!is.na(x)) / length(x)
#
#	#Assign generalists a 'turning point' at 100 % AGB removal
#	gen_turns <- turns$turn.point
#	gen_turns[turns$pval > 0.05] <- 100
#	
#	#Calculate resilience
#	scale_imp <- proportion(x = turns$turn.point)
#	mean_turns <- mean(gen_turns)
#	scale_turn <- 1- mean_turns/100
#	scale_resil <- (scale_imp + scale_turn) / 2
#	new_samp <- nrow(turns)
#	
#	outNew <- data.frame(susceptibility = scale_imp, sensitivity = scale_turn, 
#		resilience = scale_resil, num.taxa = new_samp)
#	
#	return(outNew)
#		}

###################################################################################
###################################################################################


#Function to extract susceptibility, sensitivity and resilience per functional group
resil.func <- function(func_groups, func_points, turns){
	#func_groups = functional traits data
	#func_points = output from turns called on functional group data
	#turns = output from call to 'turns'

	#Arrange and extract key data
	funcs <- arrange.funcplot(func_groups = func_groups, func_points = func_points)$summary
	func_groups <- expand.funcs(taxa = func_groups)
	
	func_set <- unique(funcs$shortname)		#Set of functions to iterate through
	
	#Output matrix
	sens <- susc <- data.frame(matrix(, nrow = 0, ncol = 2))
	for(i in 1:length(func_set)){		#For each function type
		#Extract target functional category
		func_targ <- func_set[i]
		targ_ind <- match(func_targ, names(func_groups))
		targ_col <- func_groups[ , targ_ind]
		func_type <- unique(targ_col)
		func_type <- func_type[!is.na(func_type)]
		for(j in 1:length(func_type)){		#For each functional type
			targ_taxa <- rownames(func_groups)[which(targ_col == func_type[j])]		#List of taxa belonging to that functional type
			indices <- match(targ_taxa, turns$taxon)
			if(sum(is.na(indices)) < length(indices)){
				indices <- indices[!is.na(indices)]			#Row numbers for taxa that were modelled
				sub_turns <- turns[indices, ]
				sub_turns$TaxonType <- factor(paste(func_targ, func_type[j], sep = '_'))
				#Bootstrap sensitivity and susceptibility values
				susc_boot <- boot.susc(turns_out = sub_turns)	
				#Calculate susceptibility
				func <- function(x) sum(!is.na(x))/length(x)
				prop_imp <- by(sub_turns$turn.point, factor(sub_turns$TaxonType), FUN = func)	#Proportion taxa with turning points
				susc_boot$susc$prop_imp <- prop_imp
				susc_boot$susc$num.taxa <- nrow(sub_turns)
				#Calculate sensitivity
				gen_turns <- sub_turns$turn.point
				gen_turns[is.na(gen_turns)] <- 100
				mean.turn <- 1 - by(gen_turns, sub_turns$TaxonType, FUN = mean, na.rm = TRUE) / 100	#Mean turning point
				susc_boot$sens$mean.turn <- mean.turn
				susc_boot$sens$num.taxa <- nrow(sub_turns)

				sens <- rbind(sens, susc_boot$sens)
				susc <- rbind(susc, susc_boot$susc)
				
				rm(susc_boot)
				}
			}
		}
	#Only retain functional types with >= 5 taxa
	keeps <- which(sens$num.taxa >= 5)
	sens <- sens[keeps, ]
	susc <- susc[keeps, ]
	
	#Obtain display names and key info
	align <- match(rownames(sens), funcs$taxon)
	keeps2 <- which(!is.na(align))
	sens <- sens[keeps2, ]
	susc <- susc[keeps2, ]
	align <- align[!is.na(align)]
	
	funcs_keep <- funcs[align, ]
	funcs_keep$resilience <-  (susc$prop_imp * sens$mean.turn) 
	
	out <- list(sens = sens, susc = susc, funcs = funcs_keep)
	return(out)
	}

###################################################################################
###################################################################################


#Function to bootstrap susceptibility of taxa
boot.susc <- function(turns_out){
	#turns = output from call to 'turns'; NB: must include TaxonType column
	
	types <- levels(turns_out$TaxonType)
	sens <- susc <- matrix(NA, nrow = length(types), ncol = 3)
	for(i in 1:length(types)){
		rawvals <- turnvals <- turns_out$turn.point[turns_out$TaxonType == types[i]]
		binvals <- is.na(rawvals)
		turnvals[is.na(turnvals)] <- 100
		sums <- means <- NULL
		for(k in 1:1000){
			sums[k] <- sum(sample(binvals, size = length(binvals), replace = TRUE)) / length(binvals)
			means[k] <- 1 - mean(sample(turnvals, size = length(binvals), replace = TRUE)) / 100
			}
		sens[i, ] <- 1-quantile(sums, probs = c(0.025, 0.5, 0.975))
		susc[i, ] <- quantile(means, probs = c(0.025, 0.5, 0.975))
	#	sens[i,2] <- 1-mean(binvals)
	#	susc[i,2] <- 1 - mean(turnvals)/100
		}
	sens <- data.frame(sens)
	susc <- data.frame(susc)
	rownames(sens) <- rownames(susc) <- types
	names(sens) <- names(susc) <- c('z025', 'z5', 'z975')
	
	return(list(sens = sens, susc = susc))
	}

###################################################################################
###################################################################################



#Generic bootstrapping function that returns 2.5 and 97.5 % quantiles
bootstrap <- function(x){
	#x = vector of data to bootstrap
	
	means <- numeric()
	for(i in 1:1000){
		means[i] <- mean(sample(x, length(x), replace = TRUE))
		}
	out <- quantile(means, probs = c(0.027, 0.975))
	return (out)
	}

###################################################################################
###################################################################################


#Function to estimate location of ecological thresholds
estimate.thresholds <- function(turn_points, func_points){
	#turn_points = output from call to 'fitted_mod' applied to taxa
	#func_points = output from call to 'fitted_mod' applied to functional groups
	
	#Thresholds in taxa turning points
		taxa_turn <- turn_points$dataset$turn.point[!is.na(turn_points$dataset$turn.point)]
		tt_bp <- break.points(density(taxa_turn)$x, density(taxa_turn)$y)
	#Thresholds in functional group turning points
		func_turn <- func_points$turn.point[!is.na(func_points$turn.point)]
		ft_bp <- break.points(density(func_turn)$x, density(func_turn)$y)
	#Thresholds in taxa peak rate
		taxa_rate <- turn_points$dataset$maxrate[!is.na(turn_points$dataset$maxrate)]
		tr_bp <- break.points(density(taxa_rate)$x, density(taxa_rate)$y)
	#Thresholds in functional group peak rate
		func_rate <- func_points$maxrate[!is.na(func_points$maxrate)]
		fr_bp <- break.points(density(func_rate)$x, density(func_rate)$y)

	combined <- list(tt = tt_bp, ft = ft_bp, tr = tr_bp, fr = fr_bp)
	#Extract acceleration points
	ttA <- which(combined$tt$type == 'accelerate')
	ftA <- which(combined$ft$type == 'accelerate')
	trA <- which(combined$tr$type == 'accelerate')
	frA <- which(combined$fr$type == 'accelerate')
	
	
	tt_pro <- combined$tt$bps[ttA[1]]
	tr_pro <- combined$tr$bps[trA[1]]
	ft_pro <- combined$ft$bps[ftA[1]]
	fr_pro <- combined$fr$bps[frA[1]]
	
	tt_rea <- combined$tt$bps[ttA[length(ttA)]]
	tr_rea <- combined$tr$bps[trA[length(trA)]]
	ft_rea <- combined$ft$bps[ftA[length(ftA)]]
	fr_rea <- combined$fr$bps[frA[length(frA)]]
	

	#Proactive threshold
	pros <- c(tt_pro, tr_pro, ft_pro, fr_pro)
	pros <- c(pros, mean(pros))
	
	#Reactive threshold
	react <- c(tt_rea, tr_rea, ft_rea, fr_rea)
	react <- c(react, mean(react))

	out <- data.frame(proactive = pros, reactive = react)
	out$type <- c('taxa_turn', 'taxa_rate', 'func_turn', 'func_rate', 'mean')
	
	return(out)
	}


##########################################################################################################
##########################################################################################################

#Function to work out how many datasets had repeat site visits within a survey
repeat.visits <- function(full.data){
	visits <- numeric()
	for(i in 1:length(full.data)){
		comm <- full.data[[i]]$comm.out
		reps <- summary(factor(comm$site), maxsum = nrow(comm))
		visits[i] <- mean(reps)
		}
	
	out <- data.frame(survey = names(full.data), mean.visits = visits)
	
	return(out)
	}


###################################################################################
###################################################################################


#Function to sort functional groups
sort.funcs <- function(func_groups){
	#func_groups = functional traits dataset

	funcs <- rename.funcs(func_groups = func.groups)


	funcs_sub <- funcs[ , c('category', 'qualifier', 'TaxType')]
	funcs_sub$level <- NA
		funcs_sub$level[grep('high', funcs_sub$qualifier)] <- 'high'
		funcs_sub$level[grep('medium', funcs_sub$qualifier)] <- 'medium'
		funcs_sub$level[grep('low', funcs_sub$qualifier)] <- 'low'
		funcs_sub$level <- factor(funcs_sub$level, levels = c('low', 'medium', 'high'))
	funcs_sub$qualifier <- gsub('-', '', funcs_sub$qualifier)
	funcs_sub$qualifier <- gsub('high', '', funcs_sub$qualifier)
	funcs_sub$qualifier <- gsub('medium', '', funcs_sub$qualifier)
	funcs_sub$qualifier <- gsub('low', '', funcs_sub$qualifier)
	funcs_sub$TaxType <- tolower(funcs_sub$TaxType)
	#Set factor levels for sorting
	funcs_sub$TaxType <- factor(funcs_sub$TaxType, levels = c('plant', 'invertebrate', 'arachnid', 'insect', 'ant', 'mammal', 'bird', 'reptile', 'amphibian', 'fish', 'all taxa'))
	funcs_sub$category <- factor(funcs_sub$category, levels = rev(c('Red List status', 
		'Plant', 'Habitat strata', 'Physiology', 'Development', 'Sociality', 'Movement', 'Diet', 'Trophic',  'Body mass')))
	funcs_sub$qualifier <- factor(funcs_sub$qualifier, levels = c(#'high', 'medium', 'low', 
		'aerial', 'arboreal', 'understory', 'terrestrial', 'subterranean', 'aquatic', 
		'endotherm', 'ectotherm', 'direct', 'indirect', 
		'eusocial', 'social', 'pair', 'solitary', 
		'winged', 'legged', 'legless', 
		'threatened', 'not threatened', 
		'parasitoid', 'parasite', 'hematophage', 
		'carnivore', 'piscivore', 'vertivore', 'invertivore','bacteriophage', 
		'herbivore', 'frugivore', 'granivore', 'florivore', 'nectarivore', 'palynivore', 'folivore', 
		'phloeophage', 'xylophage', 'rhizophage',
		'algivore','mycophage', 
		'saprophage', 'coprophage', 'necrophage', 'detritivore', 
		'saproxylic', 'soilphage', 
		'producer', 
		'constant', 'variable',	'generalism',
		'wood density', 'photosynthesis', 'competitor', 'stress', 'ruderal', 		
		rev(funcs$qualifier[funcs_sub$TaxType == 'plant'])))
	
	funcs_sub <- funcs_sub[, c(1,3,2,4)]
	funcs_sub <- funcs_sub[order(funcs_sub$category, funcs_sub$TaxType, funcs_sub$qualifier, funcs_sub$level), ]
	
	return(funcs_sub)
	}

###################################################################################
###################################################################################


#Bootstrap confidence intervals around ecological threshold estimates
boot.thresholds <- function(turn_points, func_points, reps = 100){
	#turn_points = output from call to 'fitted_mod' applied to taxa
	#func_points = output from call to 'fitted_mod' applied to functional groups
	#reps = number of bootstrap samples
	
	pro <- rea <- NULL
	for(i in 1:reps){
		print(paste('Running bootstrap', i, 'of', reps))
		dummy_tp <- turn_points
		tp_samp <- sample(1:nrow(turn_points), nrow(turn_points), replace = TRUE)
		fp_samp <- sample(1:nrow(func_points), nrow(func_points), replace = TRUE)
		dummy_tp$dataset <- turn_points[tp_samp, ]
		dummy_fp <- func_points[fp_samp, ]
		boot_est <- estimate.thresholds(turn_points = dummy_tp, func_points = dummy_fp)
		pro[i] <- boot_est$proactive[5]
		rea[i] <- boot_est$reactive[5]
		}
	proquant <- quantile(pro, probs = c(0.025, 0.975))
	reaquant <- quantile(rea, probs = c(0.025, 0.975))
	out <- data.frame(rbind(proquant, reaquant))
	rownames(out) <- c('proactive', 'reactive')
	return(out)
	}

###################################################################################
###################################################################################


#Function to determine which taxa are present in primary forest
primary.species <- function(thresh_data = thresh.data, lidar_data = lidar.data, 
	primary_cutoff = 95){
	#thresh_data = Site x species matrix for each dataset, containing only taxa that are duplicated across >1 dataset
	#lidar_data = Lidar data for all sites in full dataset
	#primary_cutoff = percentage of AGB required to consider a site as primary forest
	
	maxAGB <- apply(lidar_data[, 3:7], 1, max, na.omit = TRUE)
	primary <- data.frame(site = lidar_data$site, primary = as.numeric(maxAGB) >= primary_cutoff)
	prim_taxa <- NULL
	for(i in 1:length(thresh_data)){	#For each dataset
		targs <- match(thresh_data[[i]]$comm.out$site, primary$site)
		prim_sites <- thresh_data[[i]]$comm.out[which(primary$primary[targs] == TRUE), ]	#Rows that represent sites with primary forest levels of AGB
		if(nrow(prim_sites) >= 1){
			if(ncol(prim_sites) == 4){
				prim_sites$dummy <- 0
				}
			prim_taxa <- c(prim_taxa, names(which(colSums(prim_sites[, 4:ncol(prim_sites)]) >= 1)))
			}
		}
	prim_taxa <- unique(prim_taxa)
	return(prim_taxa)
	}

###################################################################################
###################################################################################
