#=================================================================
#Simple fit exploration of SR data
#Author: Catarina Wor
#Date: April 10th 2018
#=================================================================

#load in required packages and libraries
library(ggplot2)
library(TMB)

#load in directories.R
source("C:/Users/worc/Documents/HarrisonSR/R/directories.R")

#read in simple data set
setwd(data_dir)
SR<-read.csv("Harrison_simples_Apr18.csv")

summary(SR)
mydata<-list(obs_logR=log(SR$R),obs_S=SR$S_adj)

setwd(model_dir)
compile("Rickerkf.cpp")

dyn.load("Rickerkf.so")


setwd("/Users/catarinawor/Documents/work/Chinook/KalmanFilterStudy/TMB")


data <- mydata


parameters <- list(
  alphao=2.0,
  Smax = 29000,
  logSigalpha=0,
  logSigObs= -0.9,
  alpha=rep(2,length(data$obs_logR))
  )

#newtonOption(smartsearch=FALSE)
obj<-MakeADFun(data,parameters,random="alpha",DLL="Rickerkf")
obj$fn()
obj$gr()
opt<-nlminb(obj$par,obj$fn,obj$gr)
rep<-obj$report()
save(rep,file="kf.RData")













