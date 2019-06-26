#=================================================================
#Simple fit exploration of SR data
#Author: Catarina Wor
#Date: April 10th 2018
#=================================================================

#load in required packages and libraries
library(ggplot2)
library(TMB)
library(tmbstan)

#load in directories.R

#read in simple data set
SR<-read.csv("../data/Harrison_simples_Apr18.csv")


#come up with guesses based on simple ricker model
srm <- lm(log(SR$R/SR$S_adj)~ SR$S_adj)
a_srm <- srm$coefficients[1]
b_srm <- -srm$coefficients[2]
alpha <- exp(a_srm)



#set  up TMB inputs
mydata<-list(obs_logR=log(SR$R),obs_S=SR$S_adj,prbeta1=3.0,prbeta2=3.0)


parameters_recursive <- list(
  alphao=a_srm,
  logbeta = log(b_srm),
  rho=.2,
  logvarphi= 0.1,
  alpha=rep(0.9,length(SR$R))
  )


#compile TMB code
compile("Rickerkf_ratiovar.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="")
dyn.load(dynlib("Rickerkf_ratiovar"))


obj<-MakeADFun(mydata,parameters_recursive,DLL="Rickerkf_ratiovar")
opt<-nlminb(obj$par,obj$fn,obj$gr)


repkf<-obj$report()
#this is not working - hessian is not positive definite - not a big deal as we will use bayes estimates
sdreport(obj)



#===================================================================================================
#Bayes run


fitmcmc1 <- tmbstan(obj, chains=4,
          iter=100000, init="random",
          lower=c(0.1,-13.5,0.0,-3.0,rep(0.001,length(SR$R))),
          upper=c(4.0,-9.,1.0,5.0,rep(4.0,length(SR$R))),
          control = list(adapt_delta = 0.98))

mc <- extract(fitmcmc1, pars=names(obj$par),
              inc_warmup=TRUE, permuted=FALSE)
    

fit_summary <- summary(fitmcmc1)

posterior <- as.array(fitmcmc1)

mainrun <- reshape::melt(posterior)



#plot bayes time series
dfac<-data.frame(broodyear=SR$BroodYear,
  a=c(fit_summary$summary[5:34,6]), #,umsyposteriorsummary3$x[,2]
  lower=fit_summary$summary[5:34,4],
  upper=fit_summary$summary[5:34,8])


pa<-ggplot(dfac)
pa<-pa+geom_ribbon(aes(x=broodyear,ymin=lower,ymax=upper),alpha=.4)
pa<-pa+geom_line(aes(x=broodyear,y=a), size=1.2)
pa<-pa+theme_bw(16)
pa<-pa+labs(title = "Recursive Bayes model - a time series", y = expression(a[t]), x = "Brood year") 
pa

