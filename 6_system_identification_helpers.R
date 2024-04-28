# Helper functions for SWMP system identification

# Stdz function
zscore <- function(x){
  (x - mean(x))/sd(x)
}

# FUNCTION - Calculate E 
simplex_extra_fun<-function(x){
  output <- simplex(time_series = x$targetVar_norm_detrend,
                    lib=c(1,length(x$targetVar_norm_detrend)),
                    pred=c(1,length(x$targetVar_norm_detrend)),
                    E=1:10,
                    tau=-1)
  output$E[which.max(output$rho)]
  #output$E[which.min(output$mae)]
}

# FUNCTION = CALCULATE THETA - NONLINEARITY
theta_fun <- function(x){
  output <- s_map(time_series = x$targetVar_norm_detrend,
                  norm = 2,
                  lib=c(1,length(x$targetVar_norm_detrend)),
                  pred=c(1,length(x$targetVar_norm_detrend)),
                  E = min(x$E))
  output$theta[which.max(output$rho)]
  #output$theta[which.min(output$mae)]
}

# FUNCTION = CALCULATE NONLINEARITY FROM THETA's MAE
nonlin_fun <- function(x){
  output <- s_map(time_series = x$targetVar_norm_detrend,
                  norm = 2,
                  lib=c(1,length(x$targetVar_norm_detrend)),
                  pred=c(1,length(x$targetVar_norm_detrend)),
                  E = min(x$E))
  as.numeric(output$mae)[which(as.numeric(output$theta)==0)] - min(as.numeric(output$mae)) # This from tye's code
  # I *think* the following makes more sense?
  # as.numeric(output$mae)[which(as.numeric(output$theta)==0)] - as.numeric(output$mae)[which(as.numeric(output$theta)==x$theta[1])]
}

# FUNCTION = CALCULATE return linear and non-linear MAEs
MAEsOnly_fun <- function(x){
  output <- s_map(time_series = x$targetVar_norm_detrend,
                  norm = 2,
                  lib=c(1,length(x$targetVar_norm_detrend)),
                  pred=c(1,length(x$targetVar_norm_detrend)),
                  E = min(x$E))
  c(as.numeric(output$mae)[which(as.numeric(output$theta)==0)],min(as.numeric(output$mae))) # This from tye's code
  # I *think* the following makes more sense?
  # as.numeric(output$mae)[which(as.numeric(output$theta)==0)] - as.numeric(output$mae)[which(as.numeric(output$theta)==x$theta[1])]
}

# FUNCTION: CALCULATE Tp - FORECAST SKILL
forecast_fun <- function(x) {
  output <- s_map(time_series = x$targetVar_norm_detrend,
                  E = min(x$E),
                  theta = min(x$theta),
                  tp = 1)
  
  as.numeric(output$rho)[which(output$tp==1)]
}

#FUNCTION: Calculate p values for forecast skill
pred_fun <- function(x) {
  output <- s_map(time_series = x$targetVar_norm_detrend,
                  E = min(x$E),
                  theta = min(x$theta),
                  tp = 1)
  
  output$p_val[which(output$tp==1)]
}

# FUNCTION: Generate a known signal to be tested by algorithms
GenerateSignal = function(SignalType = "noise",nPoints){
  
  # SignalType = noise, linear, cycle, twocurves, breakpoint, discretelogistic, randomwalk
  if (SignalType == "noise"){
    b = 8
    myNoise = rnorm(nPoints,0,2)
    mySignal = b+myNoise
  }
  
  if (SignalType == "linear"){
    m = 0.1 # units/year
    b = 8 # intercept units
    x = 1:nPoints
    myNoise = rnorm(nPoints,0,2)
    mySignal = m*x + b
    mySignal = mySignal+myNoise
  }
  if (SignalType == "flatline"){
    m = 0 # units/year
    b = 8 # intercept units
    x = 1:nPoints
    myNoise = rnorm(nPoints,0,0.001)
    mySignal = b+myNoise
    #print(mySignal)
  }
  
  if (SignalType == "cycle"){
    m = 0.1 # units/year
    b = 5 # intercept units
    x = 1:nPoints
    # myNoise = rnorm(nPoints,0,sd(siteData$targetVar))
    myNoise = rnorm(nPoints,0,0.5)
    mySignal = sin(seq(7,6*pi,length.out=nPoints)) + b
    # plot(x,sin(seq(0,8*pi,length.out=176)),type='l')
    mySignal = mySignal+myNoise
  }
  
  # Create random vector for testing
  testCurve = FALSE
  if (SignalType == "twocurves"){
    x = 1:nPoints
    x1 = 1:floor((0.5)*nPoints)
    x2 = (floor((0.5)*nPoints)+1):nPoints
    myNoise1 = rnorm(nPoints,0,2)
    myNoise2 = runif(nPoints,0,2)
    a1 = 20
    r1 = 0.05
    r2 = 0.02
    # a(1+r)^x
    mySignal1 = a1*(1-r1)^x1
    a2 = 0.2
    mySignal = rep(NA,nPoints)
    mySignal2 = rep(NA, nPoints)
    mySignal2 = a2*(1+r2)^x2
    mySignal[x1] = mySignal1+myNoise1[x1]
    mySignal[x2] = mySignal2+myNoise2[x2]*0.1
  }
  
  if (SignalType == 'breakpoint'){
    m1 = 0.01 # units/year
    m2 = 0.5
    b1 = 8 # intercept units
    x1 = 1:floor((0.5)*nPoints)
    x2 = (floor((0.5)*nPoints)+1):nPoints
    mySignal = rep(NA,nPoints)
    mySignal2 = rep(NA, nPoints)
    myNoise1 = rnorm(length(x1),0,2)
    myNoise2 = rnorm(length(x2),0,2)
    mySignal1 = m1*x1 + b1
    b2 = -35
    mySignal2 = m2*x2 + b2
    mySignal[x1] = mySignal1 +myNoise1
    mySignal[x2] = mySignal2 + myNoise2
  }
  if (SignalType == 'discretelogistic'){
    K = 20
    r = 2.3 #1.2, 1.88, 1.9 bifurcation, 2.1 quad, 2.3 chaos  
    N = rep(NA,nPoints)
    N[1] = 3
    for (i in 2:nPoints){
      N[i] = N[i-1]^(r*(1-(N[i-1]/K)))
    }
    # plot(1:nPoints,N)
    mySignal = N
  }
  
  if (SignalType == 'randomwalk'){
    mySignal = rep(NA,nPoints)
    #myNoise = rnorm(nPoints,0,2)
    myNoise = runif(nPoints,-5,5)
    
    mySignal[1] = 10
    for (i in 2:nPoints){
      mySignal[i] = mySignal[i-1] + myNoise[i]
    }
  }
  return(mySignal)
}