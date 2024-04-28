# Run system identification analysis and print results/figures

# E = 6 means complex, E = 1 means simple
# Theta = 0 means linear, Theta > 0 may or may not be linear depending on significance distribution
# MAE = mean absolute error, bigger differences between linear/non-linear point to more significance for either linear/non-linear

# Copied over essentials from Paul's 4/8 script for the class
# and turned into a function to run per station/parameter

run_systemID_analysis <- function(swmp_data_station_param, nRandVectors = 10) {
  
  # Source the helper file that has key functions
  source('6_system_identification_helpers.R', local = TRUE)
  
  # Set up plotting region to have three rows
  par(mfrow = c(3,1))
  
  # nRandVectors = Number of random vectors for one of the statistical tests. 
  # Typically set to 1000, but that takes a long time to run. 
  # Set to 10 if just testing out the script
  
  # Declare the current station and parameter to analyze
  target_station <- unique(swmp_data_station_param$station)
  target_var <- unique(swmp_data_station_param$param)
  
  edmData <- swmp_data_station_param %>% 
    # Normalize the data
    mutate(targetVar_norm = zscore(value)) %>% 
    # Detrend the data
    mutate(targetVar_norm_detrend = astsa::detrend(targetVar_norm))
  
  message('_________________________________________\n', 
          sprintf('--> Start System ID for %s\n', unique(swmp_data_station_param$station_param)))
  
  # Plot the different versions of the data
  myMin = min(c(edmData$value, edmData$targetVar_norm, edmData$targetVar_norm_detrend))
  myMax = max(c(edmData$value, edmData$targetVar_norm, edmData$targetVar_norm_detrend))
  plot(edmData$year_frac,
       edmData$value,
       ylim=c(myMin,myMax), type='l',
       xlab='Time',
       ylab=unique(swmp_data_station_param$paramf),
       main = unique(swmp_data_station_param$station_param))
  abline(h=0,lty=2)
  lines(edmData$year_frac,edmData$targetVar_norm,col='green')
  lines(edmData$year_frac,edmData$targetVar_norm_detrend,col='blue')
  legend('topright',legend=c('Original','Normalized','Detrended'),
         col=c('black','green','blue'),lty=c(1,1,1))
  
  ##### 1: Calculate stats #####
  
  # Calculate E, embedding dimension selected as highest corr coeff, rho, 
  # between pred and obs
  edmData$E = simplex_extra_fun(edmData)
  unique(edmData$E)
  
  # Calculate Theta (non-linear tuning parameter), given E, again using rho to select
  # A test for non-linear dynamics in the data
  edmData$theta = theta_fun(edmData)
  
  # Calculate Theta (non-linear tuning parameter), given E, again using rho to select
  # A test for non-linear dynamics in the data
  edmData$theta = theta_fun(edmData)
  
  # Calculate MAE comparing theta==0 to theta for best MAE
  # as MAE (thata==0) - MAE (besttheta/smallest MAE) = how much better S-map did than just simplex
  # if MAE(theta==0) > MAE(best), then non-linear model helped
  # Get MAE of linear model and lowest MAE of non-linear model
  edm2MAEs = MAEsOnly_fun(edmData)
  edmDataMAE = nonlin_fun(edmData) # Note that this is a scalar
  
  ##### 2: Calc stats for null model #####
  
  # Calculate null distribution of ΔMAE to compare our original ΔMAE against
  # generate phase-randomized surrogates for assessing significance of nonlinear vs. linear model
  
  # Get e.g., 1000, surrogate samples from the original timeseries
  temp.series <- surrogates(edmData$targetVar_norm_detrend, nsim = nRandVectors, verbose = F)
  
  # Calculate the MAE distribution for the null model
  message("    Calculating MAEs for null model. This can take awhile if null model has lots of samples...\n")
  null_maes <- matrix(NA, ncol=nRandVectors, nrow=1)
  for(j in 1:nRandVectors){
    # Simplex to get E
    t_simp <- simplex(time_series = temp.series[,j],
                      E=1:10)
    E <- t_simp$E[which.max(t_simp$rho)]
    # s_map to get ΔMAE
    t_s_map <- s_map(time_series = temp.series[,j],
                     E = E)
    # Paul's addition
    t_s_map_theta = t_s_map$theta[which.max(t_s_map$rho)]
    # Ty's version follows
    null_maes[1, j] <- as.numeric(t_s_map$mae)[which(as.numeric(t_s_map$theta)==0)] - 
      min(as.numeric(t_s_map$mae))
  }
  
  # Get the bounds of significance from the analysis of the null model
  quant <- apply(null_maes, 1, quantile, 0.95)
  
  # Merge with estimated MAEs from the observed data
  MAEsSig <- data.frame(edmDataMAE, quant)
  
  # Categorize lakes as linear or nonlinear based on significance of ΔMAE compared to null dist'n of ΔMAEs
  linear_nonlinear <- MAEsSig %>%
    mutate(nonlin=ifelse(
      edmDataMAE > quant,
      "nonlinear",
      "linear"))
  
  # Plot of the density of null maes plus quant plus edmDataMAE
  myDensity = density(null_maes)
  plot(density(null_maes),xlim=c(min(myDensity$x),max(max(myDensity$x),edmDataMAE)),main="",ylab="Density of null MAEs")
  abline(v=quant,lty=2)
  abline(v=edmDataMAE,lty=1,col='blue')
  legendText = paste('MAE(l-nl) ', signif(edmDataMAE,3),', ',linear_nonlinear$nonlin,sep="")
  legend('topright',c(legendText,"Quantile"),lty=c(1,2),col=c('blue','black'))
  
  ##### 3. Calculate predictability #####
  
  edmSummaryInfo <- edmData %>%
    # Determine rho, which is the out of sample forecast skill
    mutate(rho=forecast_fun(.data),
           p=pred_fun(.data)) %>%
    select(station, param, station_param_f, E, theta, rho, p) %>% 
    # Get one record per station then merge with significance of nonlinearity
    distinct() %>% 
    bind_cols(linear_nonlinear) %>% 
    # Determine whether classification is significant
    mutate(nonlin = factor(nonlin),
           significant = factor(ifelse(p < 0.05, "sig", "not sig")))
  
  # Run S-map in model fitting mode (stats_only = FALSE) to get predictions
  modelPrediction <- s_map(time_series = edmData$targetVar_norm_detrend,
                           E = min(edmData$E),
                           theta = min(edmData$theta),
                           tp = 1,
                           stats_only = FALSE)
  
  # Get length of prediction vector, which can differ from the observation vector
  nModel <- length(modelPrediction$model_output[[1]]$Predictions)
  nObs <- length(edmData$year_frac)
  if (nModel>nObs){
    nStart <- 1
    nFinish <- nObs
  }else{
    nStart <- nObs-nModel+1
    nFinish <- nModel
  }
  
  # Plot the observations and predictions from the model 
  plot(edmData$year_frac[nStart:length(edmData$year_frac)], 
       modelPrediction$model_output[[1]]$Observations[1:nFinish],type='l',
       xlab='Time', 
       ylab=unique(swmp_data_station_param$paramf))
  lines(edmData$year_frac[nStart:length(edmData$year_frac)],
        modelPrediction$model_output[[1]]$Predictions[1:nFinish],col='blue')
  legend('topright',legend=c('Observed', 'Predicted',
                             sprintf('rho: %s, p: %s', 
                                     signif(edmSummaryInfo$rho,3),
                                     signif(edmSummaryInfo$p,3))),
         lty=c(1,1,NA),col=c('black','blue'))
  
  final_model_results <- paste(
    '----\n',
    sprintf('Results for %s, %s \n', target_station, target_var),
    sprintf('Simplex & s_map, E: %s, theta: %s \n', 
            edmSummaryInfo$E, 
            edmSummaryInfo$theta),
    sprintf('    MAEs(lin,nonlin,diff): %s, %s, %s \n',
            signif(edm2MAEs[1],3), 
            signif(edm2MAEs[2],3),
            signif(edm2MAEs[1]-edm2MAEs[2],3)),
    sprintf('    Model fit, rho: %s, p: %s \n', 
            signif(edmSummaryInfo$rho,3), 
            signif(edmSummaryInfo$p,3)),
    sprintf('Linear/nonlin, sig: %s, %s \n', 
            edmSummaryInfo$nonlin, 
            edmSummaryInfo$significant),
    '----\n', 
    collapse='')
  
  #message(final_model_results)
  message('    COMPLETE\n')
  return(final_model_results)
}

# Map over station and parameter to model system identification
systemID_results <- swmp_data_NAfilled %>% 
  split(.$station_param) %>% 
  map(run_systemID_analysis,
      nRandVectors = 1000)

# Clean up the global environment
rm(run_systemID_analysis)
