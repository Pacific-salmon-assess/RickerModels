#=================================================================
#Simple fit exploration of SR data
#Author: Catarina Wor
#Date: April 10th 2018
#=================================================================

#load in required packages and libraries

library(ggplot2)
library(TMB)
library(tmbstan)
#library(TMBhelper)
library(xtable)
#read in simple data set
SR <- read.csv("../data/Harrison_simples_Apr18.csv")

#=================================================================
#Simple Ricker model

# LM version
srm<-lm(log(SR$R/SR$S_adj)~ SR$S_adj)
a_srm<-srm$coefficients[1]
b_srm<--srm$coefficients[2]
alpha<-exp(a_srm)

#=============================================================================================================
#Autocorrelated recruitment model - and prior sensitivity

mydata<-list(obs_logR=log(SR$R),obs_S=SR$S_adj)

parameters_autocorr <- list(
  alpha=(a_srm),
  logbeta = log(b_srm),
  logSigObs= log(.4), 
  rho=0.3,
  delta=rep(0.0,length(SR$R))
  )


compile("Ricker_autocorr.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="",tracesweep = TRUE)
dyn.load(dynlib("Ricker_autocorr"))

obj<-MakeADFun(mydata,parameters_autocorr,random="delta",DLL="Ricker_autocorr")
  newtonOption(obj, smartsearch=FALSE)
#obj<-MakeADFun(mydata,parameters_autocorr,random=NULL,DLL="Ricker_autocorr")
#  newtonOption(obj, smartsearch=FALSE)


opt<-nlminb(obj$par,obj$fn,obj$gr)


obj$report()

sqrt(
(obj$report()$SigObs)^2+
(obj$report()$SigObs*sqrt(1-obj$report()$rho^2))^2
)


#Carrie's version
parameters_autocorr_ch <- list(
  alpha=(a_srm),
  logbeta = log(b_srm),
  logSigObs= log(.4), 
  rho=0.3
  #delta=rep(0.0,length(SR$R))
  )

dyn.unload(dynlib("Ricker_autocorr_ch"))
compile("Ricker_autocorr_ch.cpp", libtmb=FALSE, "-O1 -g", DLLFLAGS="",tracesweep = TRUE)
dyn.load(dynlib("Ricker_autocorr_ch"))

  
#removed annual random year effect:
objch<-MakeADFun(mydata,parameters_autocorr_ch,DLL="Ricker_autocorr_ch")
  newtonOption(obj, smartsearch=FALSE)
  
optch<-nlminb(objch$par,objch$fn,objch$gr)

autoRickerAIC <- 2*4-2*-opt$objective
autoRickerAICch <- 2*4-2*-optch$objective


restab <- data.frame(Parameter=c("$b$","$S_{max}$","\\alpha","$\\rho$",paste("$\\epsilon_{",SR$BroodYear,"}$")),
          catarina=c(obj$report()$beta,obj$report()$Smax,obj$report()$alpha,obj$report()$rho,obj$report()$epsilon),
          	carrie=c(objch$report()$beta,objch$report()$Smax,objch$report()$alpha,objch$report()$rhoo,objch$report()$epsilon),
          	diff=c(obj$report()$beta,obj$report()$Smax,obj$report()$alpha,obj$report()$rhoo,obj$report()$epsilon)-
          	c(objch$report()$beta,objch$report()$Smax,objch$report()$alpha,objch$report()$rhoo,objch$report()$epsilon))

xtable(restab, digits=4,caption = "Ricker recruitment autocorrelaton MLE estimates")


plot(SR$BroodYear,objch$report()$residuals, type="b",lwd=2)
lines(SR$BroodYear,obj$report()$residuals+obj$report()$delta, type="b",lwd=2, col="grey80")


#====================================================================
#projections tests

names(obj$report())
nsim<-200

et <-matrix(NA,ncol=20,nrow=nsim)
etch <- matrix(NA,ncol=20,nrow=nsim)
et[,1]<-obj$report()$delta[length(obj$report()$delta)]+obj$report()$residuals[length(obj$report()$delta)]
etch[,1]<-objch$report()$residual[length(obj$report()$residual)]
dp<-matrix(NA,ncol=20,nrow=nsim)
dpch<-matrix(NA,ncol=20,nrow=nsim)
for(i in 1:nsim){
  set.seed(i)
  dp[i,2:20] <- rnorm(19,-obj$report()$SigObs^2/2,obj$report()$SigObs)+rnorm(19,-obj$report()$SigAR^2/2,obj$report()$SigAR)
    set.seed(i)
  dpch[i,2:20] <- rnorm(19,-objch$report()$SigObs^2/2,objch$report()$SigObs)
  for(j in 2:20){
  	et[i,j] <- et[i,j-1] * obj$report()$rhoo +  dp[i,j] *sqrt(1-obj$report()$rhoo^2)
    etch[i,j] <- etch[i,j-1] * objch$report()$rhoo +  dpch[i,j] *sqrt(1-objch$report()$rhoo^2)
  }
   

}

df<-rbind(data.frame(reshape::melt(etch),model="carrie"),
	data.frame(reshape::melt(et),model="catarina"))


ggplot(df) +
geom_point(aes(x=X2,y=value),alpha=.6) +
geom_line(aes(x=X2,y=value,group =X1),alpha=.6) +
facet_wrap(~model) +
theme_bw(16)



SRdiagauto<-SR
SRdiagauto$residuals <- obj$report()$residuals
SRdiagauto$epsilon <- obj$report()$epsilon



rp <- ggplot(SRdiagauto)
rp <- rp + geom_line(aes(x=BroodYear,y=residuals),size=1.2)
rp <- rp + geom_point(aes(x=BroodYear,y=residuals, col=residuals),stroke=3)
rp <- rp + geom_hline(yintercept = 0)
rp <- rp + theme_bw(16)
rp <- rp + scale_color_viridis_c(end = 0.8)
rp



ep <- ggplot(SRdiagauto)
ep <- ep + geom_line(aes(x=BroodYear,y=epsilon),size=1.2)
ep <- ep + geom_point(aes(x=BroodYear,y=epsilon, col=epsilon),stroke=3)
ep <- ep + geom_hline(yintercept = 0)
ep <- ep + theme_bw(16)
ep <- ep + scale_color_viridis_c(end = 0.8)
ep



Splot<-seq(0,obj$report()$Smax*2)
Rpred<-Splot*exp(obj$report()$alpha)*exp(-obj$report()$beta*Splot)
df<-data.frame(Spawners=Splot,Rpred=Rpred)


p <- ggplot(df)
p <- p + geom_line(aes(x=Spawners,y=Rpred), size=2, alpha=0.6)
p <- p + geom_text(data=SR, aes(x=S_adj,y=R,label=BroodYear ),hjust=0, vjust=0)
p <- p + theme_bw(16)
p <- p + labs( x = "Spawners", y = "Recruits") 
p


#=====================================================
#Carrie




<-SR
SRdiagautoch$residuals <- objch$report()$residuals
SRdiagautoch$epsilon <- objch$report()$epsilon

SRdiagautoch$model <- "Carrie"
SRdiagauto$model <- "Catarina"

dff<-rbind(SRdiagautoch,SRdiagauto)

rp <- ggplot(SRdiagauto)
rp <- rp + geom_line(aes(x=BroodYear,y=obj$report()$delta),size=1.2)
rp <- rp + geom_point(aes(x=BroodYear,y=obj$report()$delta, col=obj$report()$delta),stroke=3)
rp <- rp + geom_hline(yintercept = 0)
rp <- rp + theme_bw(16) + facet_wrap(~model)
rp <- rp + scale_color_viridis_c(end = 0.8)
rp


data.frame(catarina=c(obj$report()$SigObs,obj$report()$rhoo),
	carrie=c(objch$report()$SigObs,objch$report()$rhoo),row.names=c("SigObs","rho"))