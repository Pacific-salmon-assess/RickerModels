#=================================================================
#Simple fit exploration of SR data for Yukon
#Author: Catarina Wor
#Date: April 10th 2018
#=================================================================

#load in required packages and libraries
library(ggplot2)
library(TMB)
#library(TMBhelper) -- not available in most recent R versions
#library(bayesplot)
#library(tmbstan)


#load in directories.R
#source("C:/Users/worc/Documents/HarrisonSR/R/directories.R")
#source("/Users/catarinawor/Documents/work/Chinook/srkf/R/directories.R")


#read in simple data set

SR<-read.csv("../data/yukon.csv")
SR <- SR[!is.na(SR$recruits)&!is.na(SR$escapement.total),]
# LM version
#simple model

srm<-lm(log(SR$recruits/SR$escapement.total)~ SR$escapement.total)
a_srm<-srm$coefficients[1]
b_srm<--srm$coefficients[2]
alpha<-exp(a_srm)

u_msy=.5*a_srm-0.07*a_srm^2





mydata<-list(obs_logR=log(SR$recruits),obs_S=SR$escapement.total)
parameters_simple <- list(
  alpha=(a_srm*.8),
  logbeta = log(b_srm*1.2),
  logSigObs= log(.4)  )
#=======================================================================================================
#=======================================================================================================
#simple model on TMB
#=======================================================================================================
#=======================================================================================================


compile("Ricker_simple.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="")
dyn.load("Ricker_simple")
  
  
obj<-MakeADFun(mydata,parameters_simple,DLL="Ricker_simple")
  
opt<-nlminb(obj$par,obj$fn,obj$gr)

rep<-obj$report()
S <-seq(0,82000, by=2000)
predR1<- S*exp(rep$alpha-rep$beta*S)

plot(S,predR1, type="l", lwd=2, ylim=c(0,max(SR$recruits)))
points(SR$escapement.total,SR$recruits)
#=============================================================================================================
#Recursive Bayes model

mydata_recursive<-list(obs_logR=log(SR$recruits),obs_S=SR$escapement.total,prbeta1=1,prbeta2=1)
parameters_recursive <- list(
  alphao=a_srm,
  logbeta = log(b_srm),
  rho=.5,
  logvarphi= 0.1,
  alpha=rep(1.9,length(mydata_recursive$obs_logR))
  )


compile("Rickerkf_ratiovar.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="")
dyn.load("Rickerkf_ratiovar")
  
  
obj<-MakeADFun(mydata_recursive,parameters_recursive,DLL="Rickerkf_ratiovar", random="alpha")

  
opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(eval.max = 1e4, iter.max = 1e4))

rep<-obj$report()
repsd<-sdreport(obj, getJointPrecision=TRUE)




predR2<-matrix(NA,nrow=length(SR$escapement.total),ncol=length(rep$alpha))

s_pred<-seq(0,max(SR$escapement.total)*1.05,length=length(SR$escapement.total))

for(i in 1:length(rep$alpha)){
 
  predR2[,i]<- s_pred*exp(rep$alpha[i]-rep$beta*s_pred)
}


df<-data.frame(S=rep(SR$escapement.total,length(rep$alpha)),
  predS=rep(s_pred,length(rep$alpha)),
    predR=c(predR2),
    R=rep(SR$recruits,length(rep$alpha)),
    a=rep(rep$alpha,each=length(rep$alpha)),
  #umsy=rep(rep$umsy,each=length(rep$umsy)),
    ayr=as.factor(rep(SR$broodyear,each=length(rep$alpha))),
    model="Kalman filter",
    BroodYear=rep(SR$broodyear,length(rep$alpha)))


df<-df[sort(order(df$S)),]



p <- ggplot(SR)
#p <- p + geom_point(aes(x=S_adj,y=R))
p <- p + geom_point(aes(x=broodyear,y=recruits), size=2, alpha=0.6)
p


p <- ggplot(df)
p <- p + geom_point(aes(x=S,y=R))
p <- p + geom_line(aes(x=predS,y=predR, color=ayr), size=2, alpha=0.6)
p <- p + geom_text(aes(x=S,y=R,label=BroodYear ),hjust=0, vjust=0)
p <- p + theme_bw(16) +scale_colour_viridis_d()
p <- p + coord_cartesian(xlim=c(0,max(df$S)))

p <- p + labs(title = "Recursive Bayes model", x = "Spawners ", y = "Recruits ", color = "Year\n") 
p
ggsave("recursive_bayes_fit.pdf", plot=p, width=10,height=7)