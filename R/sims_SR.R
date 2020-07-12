#=================================================================
#simulations to evaluate performance of TMB routines
#Author: Catarina Wor
#Date: July 12th 2018
#=================================================================

#load in required packages and libraries
library(ggplot2)
library(TMB)

#load in directories.R
source("C:/Users/worc/Documents/HarrisonSR/R/directories.R")
#source("/Users/catarinawor/Documents/work/Chinook/srkf/R/directories.R")

setwd(model_dir)
source("calc_quantile.R")
source("TMB_functions.R")


#read in simple data set
setwd(data_dir)
SR<-read.csv("Harrison_simples_Apr18.csv")

# LM version
#simple model
srm<-lm(log(SR$R/SR$S_adj)~ SR$S_adj)
a_sim<-srm$coefficients[1]
b_sim<--srm$coefficients[2]

predR<- SR$S_adj*exp(a_sim-b_sim*SR$S_adj)
sd_sim<-0.8

log(sd_sim)

ricker_simple <- function(x,a,b,sde){
  y<-x*exp(a-b*x+(rnorm(length(x),0,sd=sde)))
  return(y)
}


ricker_simple_biascorr <- function(x,a,b,sde){
  y<-x*exp(a-b*x+(rnorm(length(x),0,sd=sde))-.5*(sde^2))
  return(y)
}

k=6
nsims<-200
hr=.0

simsR<-matrix(NA, ncol=length(SR$R)+3*k,nrow=nsims)
simsS<-matrix(NA, ncol=length(SR$R)+3*k,nrow=nsims)

simsS[,1:k]<-1/b_sim
mat<-c(0.0,0.01,0.2,0.6,0.9,1)

aest<-NULL
best<-NULL
sdest<-NULL

for(i in 1:nsims){


  #initialization
   
  for(y in (1+k):ncol(simsS)){


    simsR[i,y]<-ricker_simple_biascorr(simsS[i,y-k],a_sim,b_sim,sd_sim)
    simsS[i,y]<-simsR[i,y]*(1-hr)

  }

  simdata<-list(obs_logR=log(simsR[i,(1+k*3):ncol(simsS)]),obs_S=simsS[i,(1+k*2):(ncol(simsS)-(k))])
  simparam <- list(
  alpha=(a_sim),
  logbeta = log(b_sim),
  logSigObs= log(sd_sim)  )

  simpl<-list(
  dat=simdata,
  params=simparam,
  rndm=NULL,
  dll="Ricker_simple",
  DIR=model_dir
  )

  simpleobj<-runTMB(simpl)

  sdest[i]<-simpleobj$obj$report()$SigObs
  aest[i]<-simpleobj$obj$report()$alpha
  best[i]<-simpleobj$obj$report()$beta

}


DF<-data.frame(estimated=c(sdest,aest,best), 
  simulated=rep(c(sd_sim,a_sim,b_sim),each=nsims),
  nomes=rep(c("sd","a","b"),each=nsims), 
  dummy=" ")

p<-ggplot(DF,aes(x=dummy,y=estimated))
p<-p+geom_boxplot()
p<-p+facet_wrap(~ nomes, scales="free")
p<-p+geom_hline(aes(yintercept=simulated))
p<-p+theme_bw(16)
p<-p+xlab(" parameter")
p

