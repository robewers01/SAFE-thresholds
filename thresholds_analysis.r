
#Import functions
	setwd("C:/Users/robew/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	setwd("C:/Users/rewers/OneDrive - Imperial College London/work/papers/in prep/ewers et al - SAFE meta-analysis/SAFE-thresholds/")
	source("thresholds_functions.r")

	require(safedata)
	set_safe_dir("C:/Users/rewers/OneDrive - Imperial College London/work/SAFE project - private/safedata",update = FALSE)

#Import data
	thresh.data <- readRDS("data/threshold_taxa_data.rds")	#Site x species matrix for each dataset, containing only taxa that are duplicated across >1 dataset
	lidar.data <- readRDS("data/lidar_percent.rds")		#Lidar data for all sites in full dataset

#Fit models and calculate summaries
	fitted_thresh <- fit.models(taxa_data = thresh.data, lidar = lidar.data, min.obs = 5,
		predictor = c('agb250', 'agb500', 'agb1000', 'agb2000', 'agb4000'))
		saveRDS(fitted_mod, 'results/fitted_thresh.rds')
		


#Function to estimate turning points for fitted models

fitted_mod = fitted_thresh

function(fitted_mod){

	turnpoints <- matrix(NA, nrow = length(fitted_mod$models), ncol = 2)
	for(i in 1:length(fitted_mod$models)){			#For each fitted model....
		target = fitted_mod$models[[i]]
		#Check if fitted model was significant or not
		if(as.numeric(fitted_mod$coef$pval[i]) > 0.05){
			turnpoints[i, ] <- rep(NA, 2)		#No informative turning points
			}else{
				#Find turning points
				
				
				}
		
		
		
		}
	
	
	
	
	}


####
###

#Function to calculate predicted values and derivatives

fitted.vals <- function(model){
	#model = fitted model for which values should be predicted

	cf <- as.list(coef(model$bestmod))
	names(cf) <- c("a","b")
	mod_expr <- expression((exp(a + b*predx)/(1 + exp(a + b*predx))))

	predx <- 0:100
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
	
	return(out)
	}
	
	