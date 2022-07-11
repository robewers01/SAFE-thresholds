
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
fitted.vals <- function(model, mod_coef){
	#model = fitted model for which values should be predicted
	#mod_coef = summary model information (passed for follow-up functions)

	cf <- as.list(coef(model$bestmod))
	names(cf) <- c("a","b")
	mod_expr <- expression((exp(a + b*predx)/(1 + exp(a + b*predx))))

	predx <- seq((-150), 150, 0.1)
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
	
	past <- NA
	try(past <- turnpoints(y)$tp, silent = TRUE)
	try(peak <- turnpoints(y)$firstispeak, silent = TRUE)
	if(!is.na(past[1])){
		first <- x[past]
		}else{
			first <- numeric()
			peak <- NA
			}
		
	return(list(tps = min(first), peak = peak, coefs = fitted_vals$coefs))		
	}

###################################################################################
###################################################################################


#Function to extract turning point data
extract.turns <- function(turns_estimate){
	#turns_estimate = output from root.finder

	obs.x.range <- as.numeric(c(turns_estimate$coefs[, 10:11]))
	turnings <- turns$tps
	
	#If no turning point was detected:
	if(length(turnings) == 0 & as.numeric(turns_estimate$coefs$pval) < 0.05){			
		turnings <- 0	#Set turning points to be the min (because fitted model is significant)
		}

	#Check for turning points that fall outside range of valid AGB values
	if(as.numeric(turns_estimate$coefs$slope) > 0){
		if(turns_estimate$peak == FALSE) turnings <- 0		#Turning point returned has missed initial acceleration point
		}
	if(as.numeric(turns_estimate$coefs$slope) < 0){ 
		if(turns_estimate$peak == TRUE) turnings <- 0		#Turning point returned has missed initial deceleration point
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



