
library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)

#--------------------------------------------------------------------------------------------
# For  Ricker AR1, mulitiple stock
SRDat <- read.csv("../data/Harrison_simples_Apr18.csv")
#SRDat_ar <- filter(SRDat, CU_Name == "Lower_Thompson"|CU_Name == "South_Thompson")
#SRDat_std <- filter(SRDat, CU_Name != "Lower_Thompson" & CU_Name != "South_Thompson")

TMB_Inputs <- list(Scale = 1000, logA_Start = 1, rho_Start = 0.1, Sgen_sig = 1)


#Set up data and parameter lists for input into TMB model
Scale <- TMB_Inputs$Scale
data <- list()
# Data 
data$S <- SRDat$S_adj/Scale 
data$logR <- log(SRDat$R/Scale)
data$stk <- rep(1,length(SRDat$R))#as.numeric(SRDat$CU_ID)
N_Stocks <- 1#length(unique(SRDat$CU_Name))
data$yr <- seq_along(data$stk)#SRDat$yr_num
#data$N_Stks <- N_Stocks
#data$model <- rep(0,N_Stocks)
#data$model[3] <- 1 #3rd stock has AR(1)
#data$Sgen_sig <- TMB_Inputs$Sgen_sig

param <- list()
# Parameters for stocks with AR1
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat  %>% summarise(x=quantile(S_adj, 0.8)))$x/Scale) )
param$rho <- rep(TMB_Inputs$rho_Start, N_Stocks)
param$logSigma <- rep(-2, N_Stocks)

# Parameters for stocks without AR1
#param$logA_std <- rep(TMB_Inputs$logA_Start, N_Stocks)
#param$logB_std <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
#param$logSigma_std <- rep(-2, N_Stocks)
#param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 

# Compile model if changed:
dyn.unload(dynlib("TMB_Files/Ricker_ar1_single"))
compile("Ricker_ar1_single.cpp")

dyn.load(dynlib("Ricker_ar1_single"))

# For Phase 1, fix Sgen
#map <- list(logSgen=factor(rep(NA, N_Stocks))) # Determine which parameters to fix
#obj <- MakeADFun(data, param, DLL="Ricker_ar1", silent=TRUE, map=map)
obj <- MakeADFun(data, param, DLL="Ricker_ar1_single", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
summary(sdreport(obj))
exp(summary(sdreport(obj))["logSigma",])
(summary(sdreport(obj))["rho",])

names(sdreport(obj))

# For Phase 2, pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
upper <- c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf )

pl$logSgen <- log(0.3*SMSYs)


obj <- MakeADFun(data, pl, DLL="Ricker_ar1", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )
pl <- obj$env$parList(opt$par) 
summary(sdreport(obj))


#--------------------------------------------------------------------------------------------
# For  Ricker AR1, single stock
#SRDat <- read.csv("DataIn/SRDat.csv")
SRDat <- read.csv("DataIn/SRDat_Harrison.csv")
SRDat <- SRDat %>% filter (CU_ID==0)

TMB_Inputs <- list(Scale = 1000, logA_Start = 1, rho_Start = 0.1)


#Set up data and parameter lists for input into TMB model
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$yr <- SRDat$yr_num


param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$rho <- rep(TMB_Inputs$rho_Start, N_Stocks)
param$logSigma <- rep(-2, N_Stocks)


# Compile model if changed:
dyn.unload(dynlib("TMB_Files/Ricker_ar1_single"))
compile("TMB_Files/Ricker_ar1_single.cpp")

dyn.load(dynlib("TMB_Files/Ricker_ar1_single"))

obj <- MakeADFun(data, param, DLL="Ricker_ar1_single", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) 
summary(sdreport(obj))


#--------------------------------------------------------------------------------------------
#R functions for Ricker AR1, single stock, with SMSY, SGEN, SREP
RickerAR1.model <- function(theta,R,S){
  a <- exp(as.numeric((theta[1])))
  b <- exp(as.numeric(theta[2]))
  rho <- as.numeric(theta[3])
  sig <- exp(as.numeric(theta[4]))
  n <- length(S)
  logPR <- NA
  dev <- NA
  for (i in 1:n){ 
    if (i==1) {
      logPR[i] <- log(a) + log(S[i]) - b * S[i]
      dev[i] <- log(R[i]) - logPR[i]
      }
    if (i>=2) {
      logPR[i] <- log(a) + log(S[i]) -b * S[i] + rho * dev[i-1]
      dev[i] <- log(R[i]) - logPR[i]
      }
  }
  epsilon.wna <- log(R) - logPR	#residuals
  epsilon <- as.numeric(na.omit(epsilon.wna))
  nloglike <- -sum(dnorm(epsilon,0,sig, log=T))
  #return(list(PR=PR, epsilon=epsilon, nloglike=nloglike)) 
  return(nloglike=nloglike) 
}

RickerAR1.solver <- function(R,S){
  #init.vals<-c(1,-2.4,0.1,-2)
  init.vals<-c(1,-7,0.1,-2)#good starting values for Harrison
  
  SRfit=optim(par=init.vals,fn=RickerAR1.model,R=R,S=S, method="BFGS", hessian=T)	#CH: hessian=2nd derivative, optim good for up to 40 parameters
  V=solve(SRfit$hessian) #covariance matrix
  std=sqrt(abs(diag(V)))
  X=V/(std%o%std)	#outer product of standard dev matrix (CH comment)
  return(list(SRfit=SRfit, etheta=SRfit$par, V=V, X=X))
}

SRDat <- read.csv("DataIn/SRDat.csv")
SRDat <- SRDat %>% filter (CU_ID==2)
SRDat <- read.csv("DataIn/SRDat_Harrison.csv")
SRDat <- SRDat %>% filter (CU_ID==0)
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
ch <- RickerAR1.solver(exp(data$logR), data$S)$SRfit$par
ch
sm <- (1 - lambert_W0( exp (1  - ch[1]))) / exp(ch[2]) #SMSY
sm.adj <- (1 - lambert_W0( exp (1  - (ch[1] + exp(ch[4])^2/2 )))) / exp(ch[2]) #SMSY adjusted for log-normal transformation bias
sg <- Sgen.solver(exp(ch[1]), exp(ch[2]),1) #SGEN
sg.adj <- Sgen.solver(exp(ch[1]+ exp(ch[4])^2/2), exp(ch[2]),1) #SGEN
sr <- ch[1]/exp(ch[2]) #SREP
sr.adj <- (ch[1] + exp(ch[4])^2/2 )/exp(ch[2]) #SREP

#--------------------------------------------------------------------------------------------
# For  Ricker, multiple stocks
# gives same results at IFRoutTo2013rec.rds for Lower Thompson (CU_ID=2)

SRDat <- read.csv("DataIn/SRDat.csv")
#SRDat <- SRDat %>% filter(Spawners!=9522) #Try removing 4th year of CU_ID 2 (i.e., if it were NA)

TMB_Inputs <- list(Scale = 1000, logA_Start = 1)


#Set up data and parameter lists for input into TMB model
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Spawners/Scale 

data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
#data$N_Stks <- N_Stocks
#data$yr <- SRDat$yr_num


param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma <- rep(-2, N_Stocks)


# Compile model if changed:
dyn.unload(dynlib("TMB_Files/Ricker"))
compile("TMB_Files/Ricker.cpp")

dyn.load(dynlib("TMB_Files/Ricker"))

obj <- MakeADFun(data, param, DLL="Ricker", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))


#--------------------------------------------------------------------------------------------
# R code for Ricker, single stock
Ricker.model <- function(theta,R,S){
  a <- exp(as.numeric((theta[1])))
  b <- exp(as.numeric(theta[2]))
  sig <- exp(as.numeric(theta[3]))
  PR=a*S*exp(-b*S)
  epsilon.wna=log(R)-log(PR)	#residuals
  epsilon=as.numeric(na.omit(epsilon.wna))
  nloglike = -sum(dnorm(epsilon,0,sig, log=T))
  #return(list(PR=PR, epsilon=epsilon, nloglike=nloglike)) #actually returns postive loglikelihood
  return(nloglike=nloglike)
}


Ricker.solver <- function(R,S){
  init.vals <- c(1,-2.4,-2)
  SRfit = optim(par=init.vals,fn=Ricker.model,R=R,S=S, method="BFGS", hessian=T)	#CH: hessian=2nd derivative, optim good for up to 40 parameters
  V = solve(SRfit$hessian) #covariance matrix
  std = sqrt( abs( diag(V) ) )
  X = V / (std %o% std)	#outer product of standard dev matrix (CH comment)
  return(list(SRfit=SRfit, etheta=SRfit$par, V=V, X=X))
}

SRDat <- read.csv("DataIn/SRDat.csv")
SRDat <- SRDat %>% filter (CU_ID==2)
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Spawners/Scale 
#data$S[4] <- NA #Try putting NA in 4th year of Lower Thompson (CU_ID=2) 
data$logR <- log(SRDat$Recruits/Scale)

ch <- Ricker.solver(exp(data$logR), data$S)$SRfit$par
ch

SRDat <- read.csv("DataIn/SRinputfile.csv")
SRDat <- SRDat %>% filter (Stocknumber==0)#2)
Scale <- TMB_Inputs$Scale
data <- list()
data$S <- SRDat$Sp/Scale 
data$logR <- log(SRDat$Rec/Scale)


ch <- Ricker.solver(exp(data$logR), data$S)$SRfit$par
ch

#------------------------------------------------------------------------------------------------------------------
require(gsl)

Sgen.model<-function(S,a,b,sig){
  PR<-a*S*exp(-b*S)
  SMSY <- ( 1 - lambert_W0( exp(1 - log(a)))) / b 
  epsilon.wna <- log(SMSY)-log(PR)	#residuals
  epsilon <- as.numeric(na.omit(epsilon.wna))
  nloglike <- -sum(dnorm(epsilon,0,sig, log=T))
  if(is.na(sum(dnorm(epsilon,0,sig, log=T)))==TRUE) print(c(a,b,sig))
  return(nloglike)
  
}

Sgen.solver <- function(a,b,sig) {
  SMSY <- (1 - lambert_W0( exp (1  - log(a)))) / b #(log(a)/b)*(0.5-0.07*log(a))
  SRfit=optimize(f=Sgen.model,interval=c(0, SMSY), a=a, b=b, sig=sig)	 # nb: not optim() !!
  return(list(SRfit=SRfit$minimum))  # returns the minimum S
}
