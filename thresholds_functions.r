


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
				rm(glms)
				}
			}
		}
	
	#Remove any models that weren't fully fitted or tested
	if(!is.na(coefs)[1]){
		coefs <- coefs[!is.na(coefs$pval), ]
		}
		
	return(coefs)
	}

###################################################################################
###################################################################################
