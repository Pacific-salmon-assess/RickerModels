#=================================================================
#Simple fit exploration of SR data
#Author: Catarina Wor
#Date: April 10th 2018
#=================================================================

#load in required packages and libraries
library(ggplot2)
library(TMB)
library(TMBhelper)

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

#============================================================
#Rmissing<-SR$R
#Rmissing[which(is.na(SR$Age2Surv))]<-NA
#Smissing<-SR$S_adj
#Smissing[which(is.na(SR$Age2Surv))]<-NA
#summary(SR)
#

#========================================
#kalman filter model
mydata<-list(obs_logR=log(SR$R),obs_S=SR$S_adj)
setwd(model_dir)
compile("Rickerkf.cpp", "-O1 -g", DLLFLAGS="")

dyn.load("Rickerkf.so")


data <- mydata


parameters <- list(
  alphao=a_srm,
  logbeta = log(b_srm),
  logSigalpha=log(.5),
  logSigObs= log(.4),
  alpha=rep(0.9,length(data$obs_logR))
  )

obj<-MakeADFun(data,parameters,random="alpha",DLL="Rickerkf")
newtonOption(obj, smartsearch=FALSE)
obj$fn()
obj$gr()

optkf<-nlminb(obj$par,obj$fn,obj$gr)
rep<-obj$report()
rep


TMBAIC(optkf, p = 2, n = Inf)
predR2<-matrix(NA,nrow=length(SR$S_adj),ncol=length(rep$alpha))

for(i in 1:length(rep$alpha)){
	predR2[,i]<- SR$S_adj*exp(rep$alpha[i]-rep$beta*SR$S_adj)
}


df<-data.frame(S=rep(SR$S_adj,length(rep$alpha)),
		predR=c(predR2),
		R=rep(SR$R,length(rep$alpha)),
		a=rep(rep$alpha,each=length(rep$alpha)),
		ayr=as.factor(rep(SR$BroodYear,each=length(rep$alpha))),
		model="Kalman filter",
		BroodYear=rep(SR$BroodYear,length(rep$alpha)))


df<-df[sort(order(df$S)),]

head(df)


p <- ggplot(df)
#p <- p + geom_point(aes(x=S_adj,y=R))
p <- p + geom_line(aes(x=S,y=predR, color=ayr), size=2, alpha=0.6)
p <- p + geom_text(aes(x=S,y=R,label=BroodYear ),hjust=0, vjust=0)
p <- p + theme_bw(16)
p <- p + labs(title = "Kalman Filter model", x = "Spawners", y = "Recruits", color = "Year\n") 
#ylab("Recruits") + xlab("Spawners")
p



dfa<-data.frame(broodyear=SR$BroodYear,alpha=exp(rep$alpha))


pa<-ggplot(dfa)
pa<-pa+geom_line(aes(x=broodyear,y=alpha), size=2)
pa<-pa+geom_point(aes(x=broodyear,y=alpha),size=4)
pa<-pa+theme_bw(16)
pa<-pa+labs(title = "Kalman Filter model - alpha time series", y = expression(alpha), x = "Brood year") 
pa


###################
#MCMC
###################


fitmcmc2 <- tmbstan(obj, chains=4,
              iter=1000000, init=list(opt$par),
              lower=c(0.1,-13.9,-3.0), upper=c(4.0,-9.214608,0.25),
               control = list(adapt_delta = 0.89))


#===================================================================
#kf with process error

mydata<-list(obs_logR=log(SR$R),obs_S=SR$S_adj)
setwd(model_dir)
compile("Rickerkf_procerr.cpp")

dyn.load("Rickerkf_procerr.so")


data <- mydata


parameters <- list(
  alphao=1.27,
  beta =  4.025214e-06,
  #148434,
  logSigalpha=log(.9),
  logSigObs= log(1.4),
  logSigProc= log(1.5),
  alpha=rep(0.9,length(data$obs_logR)),
  v=rep(0.,length(data$obs_logR))
  )

#newtonOption(smartsearch=FALSE)
objpe<-MakeADFun(data,parameters,random=c("alpha","v"),DLL="Rickerkf_procerr")
objpe$fn()
objpe$gr()
optkfpe<-nlminb(objpe$par,objpe$fn,objpe$gr)
reppe<-objpe$report()

reppe

plot(SR$BroodYear,exp(rep$alpha), lwd=2, type="o")
lines(SR$BroodYear,exp(reppe$alpha), lwd=2, type="o",col="blue")
#save(rep,file="kf.RData")

TMBAIC(optkf, p = 2, n = Inf)









#===================================================================
#model with survival up t age 2 as a covariate



mydata_surv<-list(obs_logR=log(SR$R[!is.na(SR$Age2Surv)]),obs_S=SR$S_adj[!is.na(SR$Age2Surv)],
	surv=SR$Age2Surv[!is.na(SR$Age2Surv)])

setwd(model_dir)
compile("Ricker_surv.cpp")

dyn.load("Ricker_surv.so")


 mydata_surv

parameters_surv <- list(
  alphao=0.91298,
  Smax = 384818.7,
  logSigObs= 0
  )

#newtonOption(smartsearch=FALSE)
obj_surv<-MakeADFun(mydata_surv,parameters_surv,DLL="Ricker_surv")
obj_surv$fn()
obj_surv$gr()
opt<-nlminb(obj_surv$par,obj_surv$fn,obj_surv$gr)
rep_surv<-obj_surv$report()

rep_surv
exp(rep_surv$alphao)


TMBAIC(opt, p = 2, n = Inf)








