#=================================================================
#Simple fit exploration of SR data
#Author: Catarina Wor
#Date: April 10th 2018
#=================================================================

#load in required packages and libraries
library(ggplot2)
library(TMB)
#library(TMBhelper)
library(bayesplot)
library(tmbstan)
library(reshape)

#load in directories.R
source("C:/Users/worc/Documents/HarrisonSR/R/directories.R")
#source("/Users/catarinawor/Documents/work/Chinook/srkf/R/directories.R")

#read in simple data set
setwd(data_dir)
SR<-read.csv("Harrison_simples_Apr18.csv")

#simple model
srm<-lm(log(SR$R/SR$S_adj)~ SR$S_adj)
a_srm<-srm$coefficients[1]
b_srm<--srm$coefficients[2]
alpha<-exp(a_srm)

predR1<- SR$S_adj*exp(a_srm-b_srm*SR$S_adj)
df<-data.frame(S=SR$S_adj,predR=predR1,R=SR$R,model="simple",BroodYear=SR$BroodYear)
df<-df[order(sort(SR$S_adj)),]

#========================================
#kalman filter model
mydata<-list(obs_logR=log(SR$R),obs_S=SR$S_adj)
setwd(model_dir)
compile("Rickerkf_ratiovar.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="")

dyn.load(dynlib("Rickerkf_ratiovar"))



data <- mydata


parameters <- list(
  alphao=a_srm,
  logbeta = log(b_srm),
  rho=.5,
  logvarphi= 0.1,
  alpha=rep(0.9,length(data$obs_logR))
  )

objkf<-MakeADFun(data,parameters,random="alpha",DLL="Rickerkf_ratiovar")
newtonOption(objkf, smartsearch=FALSE)
objkf$fn()
objkf$gr()

optkf<-nlminb(objkf$par,objkf$fn,objkf$gr)
repkf<-objkf$report()
repkf


#TMBAIC(optkf, p = 2, n = Inf)
predR2<-matrix(NA,nrow=length(SR$S_adj),ncol=length(repkf$alpha))

for(i in 1:length(repkf$alpha)){
	predR2[,i]<- SR$S_adj*exp(repkf$alpha[i]-repkf$beta*SR$S_adj)
}


df<-data.frame(S=rep(SR$S_adj,length(repkf$alpha)),
		predR=c(predR2),
		R=rep(SR$R,length(repkf$alpha)),
		a=rep(repkf$alpha,each=length(repkf$alpha)),
		ayr=as.factor(rep(SR$BroodYear,each=length(repkf$alpha))),
		model="Kalman filter",
		BroodYear=rep(SR$BroodYear,length(repkf$alpha)))




###################
#MCMC
###################


fitmcmc2 <- tmbstan(objkf, chains=3,
              iter=100000, init="random",
              lower=c(0.1,-13.9,0.0,-3.0), upper=c(4.0,-9.214608,1.0,5.0),
               control = list(adapt_delta = 0.98))

mc <- extract(fitmcmc2, pars=names(objkf$par),
              inc_warmup=TRUE, permuted=FALSE)
npar<-dim(mc)[3]

rhats <- stan_rhat(fitmcmc2,par=c(names(objkf$par),"alpha"))
stan_ess(fitmcmc2,par=c(names(objkf$par),"alpha"))

stan_diag(fitmcmc2, 
            information = 'sample', 
            chain = 0)




fit_summary <- summary(fitmcmc2)
fit_summary$summary

posterior <- as.array(fitmcmc2)
dim(posterior)

head(posterior[,1,2])

posterior[1:3,,]


mainrun<-melt(posterior)
head(mainrun)
summary(mainrun)
mainrun$type="posdata"
mainrun$type<-as.factor(mainrun$type)


pm<-ggplot(mainrun)
pm<-pm+geom_density(aes(x=value, color=chains))
pm<-pm+facet_wrap(~parameters, scales="free")
pm
?density

#mcmc_dens_overlay(posterior, pars =  c(names(objkf$par)))
#mcmc_pairs(posterior, pars =  c(names(objkf$par)),
#          off_diag_args = list(size = 1.5))


#==================================================================
#Post model pre data
#==================================================================


setwd(model_dir)
compile("Rickerkf_ratiovar_predata.cpp")

dyn.load(dynlib("Rickerkf_ratiovar_predata"))

data <- mydata

parameters <- list(
  alphao=a_srm,
  logbeta = log(b_srm),
  rho=.5,
  logvarphi= 0.1,
  alpha=rep(0.9,length(data$obs_logR))
  )

objkfpd<-MakeADFun(data,parameters,random="alpha",DLL="Rickerkf_ratiovar_predata")
newtonOption(objkfpd, smartsearch=FALSE)
objkfpd$fn()
objkfpd$gr()

optkfpd<-nlminb(objkfpd$par,objkfpd$fn,objkfpd$gr)
repkfpd<-objkfpd$report()
repkfpd

fitmcmc2pd <- tmbstan(objkfpd, chains=3,
              iter=100000, init="random",
              lower=c(0.1,-13.9,0.0,-3.0), upper=c(4.0,-9.214608,1.0,5.0),
               control = list(adapt_delta = 0.98))



mc <- extract(fitmcmc2pd, pars=names(objkfpd$par),
              inc_warmup=TRUE, permuted=FALSE)
npar<-dim(mc)[3]

dim(mc)


posteriorpd <- as.array(fitmcmc2pd)
dim(posteriorpd)
head(posteriorpd)

predata<-melt(posteriorpd)
predata$type="predata"
predata$type<-as.factor(predata$type)


allruns<-rbind(predata,mainrun)


myparam<-unique(allruns$parameters)



for(i in 1:length(myparam)){
  

  if(myparam[i]=="logbeta"){
    resalpha<-allruns[allruns$parameters==myparam[i],]

    pao<-ggplot(resalpha)
    pao<-pao+geom_density(aes(x=exp(value), color=type,fill=type),alpha=0.5)
    pao<-pao+xlab("beta")
    print(pao)

    setwd(figs_dir)

    ggsave(paste("pmpd","beta",".pdf", sep=""), plot = pao, width = 7, height = 7, units = "cm")
    
    pao<-ggplot(resalpha)
    pao<-pao+geom_density(aes(x=1/exp(value), color=type,fill=type),alpha=0.5)
    pao<-pao+xlab("Smax")
    print(pao)

    setwd(figs_dir)

    ggsave(paste("pmpd_gamma","Smax",".pdf", sep=""), plot = pao, width = 10, height = 10, units = "cm")
    

  }else{


  resalpha<-allruns[allruns$parameters==myparam[i],]

  pao<-ggplot(resalpha)
  pao<-pao+geom_density(aes(x=value, color=type,fill=type),alpha=0.5)
  pao<-pao+xlab(myparam[i])
  print(pao)

  setwd(figs_dir)

  ggsave(paste("pmpd_gamma",myparam[i],".pdf", sep=""), plot = pao, width = 10, height = 10, units = "cm")
  }
}












mcmc_dens_overlay(posteriorpd, pars =  c(names(objkfpd$par)))
mcmc_pairs(posteriorpd, pars =  c(names(objkfpd$par)),
           off_diag_args = list(size = 1.5))


