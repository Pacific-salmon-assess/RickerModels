#=================================================================
#Simple fit exploration of SR data
#Author: Catarina Wor
#Date: April 10th 2018
#=================================================================

#load in required packages and libraries
library(ggplot2)
library(TMB)
#library(TMBhelper) -- not available in most recent R versions
library(bayesplot)
library(tmbstan)
library(reshape)
library(xtable)

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
a_srm<-srm$coefficients[1]
b_srm<--srm$coefficients[2]
alpha<-exp(a_srm)

predR1<- SR$S_adj*exp(a_srm-b_srm*SR$S_adj)

mydata<-list(obs_logR=log(SR$R),obs_S=SR$S_adj)
parameters_simple <- list(
  alpha=(a_srm),
  logbeta = log(b_srm),
  logSigObs= log(.4)  )
#=======================================================================================================
#=======================================================================================================
#simple model on TMB
#=======================================================================================================
#=======================================================================================================
simpl<-list(
  dat=mydata,
  params=parameters_simple,
  rndm=NULL,
  dll="Ricker_simple",
  DIR=model_dir
  )

simpleobj<-runTMB(simpl)

simpleobj$report()

#MCMC
simpleB<-list(
  obj=simpleobj,
  nchain=3,
  iter=1000000,
  lowbd=c(0.1,-13.5,-6.0),
  hibd=c(4.0,-8.5,5.0)
)

posterior_simple<-posteriorsdf(simpleB)

#plots
plot_posteriors(posterior_simple$posteriors)

#posteriors of derived quantities -- the interesting ones
simpdf<-posterior_simple$posteriors

a<-(simpdf$value[simpdf$parameters=="alpha"])
alpha<-exp(simpdf$value[simpdf$parameters=="alpha"])
beta<-exp(simpdf$value[simpdf$parameters=="logbeta"])
Smax<-1/exp(simpdf$value[simpdf$parameters=="logbeta"])
sig<-exp(simpdf$value[simpdf$parameters=="logSigObs"])

deriv_posteriors<-data.frame(chains=rep(simpdf$chains[simpdf$parameters=="logbeta"],4),
                             parameters = rep(c("a","b","Smax","sig"),each=length(a)),
                             value = c(a,beta,Smax,sig)
                             )

plot_posteriors(deriv_posteriors,salvar=TRUE,DIR=figs_dir,nome="posterior_simple_model_fit.pdf")


Di<-list(
  DIR=tex_dir,
  param_names=c("a","b","$\\alpha$","$S_{max}$"),
  MLE=c(a_srm,b_srm,exp(a_srm),1/b_srm),
  MCMC=cbind(a,beta,alpha,Smax),
  caption = "Parameter estimates for traditional Ricker model.",
  digits=matrix(c(2,-2,2,2),ncol=6,nrow=4),
  other=NULL,
  filename="simple_tab.tex")

 results_table(Di) 

# mode -- should be the equivalent to the MLE
density(a)$x[which.max(density(a)$y)]


#Model predictions
pred_bayes<-matrix(NA,ncol=length(SR$S_adj),nrow=length(a))

SR$S_adj*exp(a[1]+beta[1]*SR$S_adj)

for(i in 1:length(a)){
  pred_bayes[i,]<-SR$S_adj*exp(a[i]-beta[i]*SR$S_adj)
}

M<-list(
  predBayes=pred_bayes,
  predFreq=predR1,
  orig_data=SR
  )

model_pred_plot(M, salvar=TRUE,DIR=figs_dir,filename="simple_model_fit.pdf")

#=============================================================================================================
#Recursive Bayes model


parameters_recursive <- list(
  alphao=a_srm,
  logbeta = log(b_srm),
  rho=.5,
  logvarphi= 0.1,
  alpha=rep(0.9,length(mydata$obs_logR))
  )

recursive<-list(
  dat=mydata,
  params=parameters_recursive,
  rndm="alpha",
  dll="Rickerkf_ratiovar",
  DIR=model_dir
  )

recursiveobj<-runTMB(recursive)
repkf<-recursiveobj$report()


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


df<-df[sort(order(df$S)),]

head(df)

p <- ggplot(df)
#p <- p + geom_point(aes(x=S_adj,y=R))
p <- p + geom_line(aes(x=S,y=predR, color=ayr), size=2, alpha=0.6)
p <- p + geom_text(aes(x=S,y=R,label=BroodYear ),hjust=0, vjust=0)
p <- p + theme_bw(16)
p <- p + labs(title = "Recursive Bayes model", x = "Spawners", y = "Recruits", color = "Year\n") 
#ylab("Recruits") + xlab("Spawners")
p



recursiveB<-list(
  obj=recursiveobj,
  nchain=3,
  iter=1000000,
  lowbd=c(0.1,-13.5,0.0,-3.0),
  hibd=c(4.0,-9.,1.0,5.0)
)


posterior_recursive<-posteriorsdf(recursiveB)

plot_posteriors(posterior_recursive$posteriors)

recrsdf<-posterior_recursive$posteriors


a<-(recrsdf$value[recrsdf$parameters=="alphao"])
beta<-exp(recrsdf$value[recrsdf$parameters=="logbeta"])
Smax<-1/exp(recrsdf$value[recrsdf$parameters=="logbeta"])
rho<-(recrsdf$value[recrsdf$parameters=="rho"])

deriv_posteriors<-data.frame(chains=rep(recrsdf$chains[recrsdf$parameters=="logbeta"],4),
                             parameters = rep(c("ao","b","Smax","rho"),each=length(a)),
                             value = c(a,beta,Smax,rho)
                             )

plot_posteriors(deriv_posteriors,salvar=TRUE,DIR=figs_dir,nome="posterior_recursive_model.pdf")

#prior on rho
pdf("prior_rho.pdf")
plot(seq(0,1,by=.05),dbeta(seq(0,1,by=.05),3,3), type="l",
  xlab=expression(paste(rho," value")), ylab= "density", las=1, lwd=3)
dev.off()



Dr<-list(
  DIR=tex_dir,
  param_names=c("b","$S_{max}$","$\\rho$",paste("a",SR$BroodYear)),
  MLE=c(repkf$beta,repkf$Smax,repkf$rho,repkf$alpha),
  MCMC=cbind(beta,Smax,rho),
  caption = "Parameter estimates for recursive Bayes Ricker model.",
  digits=matrix(c(2,-2,2,2),ncol=6,nrow=length(SR$BroodYear)+2),
  other=rbind(posterior_recursive$fit_summary$summary[5:34,4],posterior_recursive$fit_summary$summary[5:34,6],posterior_recursive$fit_summary$summary[5:34,8]),
  filename="recursive_tab.tex")

 results_table(Dr) 



dfa<-data.frame(broodyear=SR$BroodYear,a=c(repkf$alpha,posterior_recursive$fit_summary$summary[5:34,6]), 
  type=rep(c("MLE","Bayes"),each=length(SR$BroodYear)),lower=c(rep(NA,length(SR$BroodYear)),posterior_recursive$fit_summary$summary[5:34,4]),
  upper=c(rep(NA,length(SR$BroodYear)),posterior_recursive$fit_summary$summary[5:34,8]))


pa<-ggplot(dfa)
pa<-pa+geom_line(aes(x=broodyear,y=a,col=type), size=2)
pa<-pa+geom_ribbon(aes(x=broodyear,ymin=lower,ymax=upper, fill=type),alpha=.4)
pa<-pa+theme_bw(16)
pa<-pa+labs(title = "Recursive Bayes model - time series", y = expression(a), x = "Brood year") 
pa
setwd(figs_dir)
ggsave("recursive_a.pdf", plot=pa, width=10,height=7)


#=============================================================================================================
#=============================================================================================================
#Model with survival


mydatas<-list(obs_logR=log(SR$R),obs_S=SR$S_adj, obs_survival=SR$Age2Surv)


parameters_survival <- list(
  alpha=a_srm,
  logbeta = log(b_srm),
  logSigObs= log(0.9)
  )

survival<-list(
  dat=mydatas,
  params=parameters_survival,
  rndm=NULL,
  dll="Ricker_survival",
  DIR=model_dir
  )


survivalobj<-runTMB(survival)
repsurv<-survivalobj$report()


#TMBAIC(optkf, p = 2, n = Inf)
predR3<-matrix(NA,nrow=length(SR$S_adj),ncol=length(repsurv$pred_logR))

for(i in 1:length(repsurv$pred_logR)){
  predR3[,i]<- SR$S_adj*exp(repsurv$alpha+log(SR$Age2Surv[i])-repsurv$beta*SR$S_adj)
}


df<-data.frame(S=rep(SR$S_adj,length(repsurv$pred_logR)),
    predR=c(predR3),
    R=rep(SR$R,length(repsurv$pred_logR)),
    a=rep(repkf$alpha,each=length(repsurv$pred_logR)),
    ayr=as.factor(rep(SR$BroodYear,each=length(repsurv$pred_logR))),
    model="recursive Bayes",
    BroodYear=rep(SR$BroodYear,length(repsurv$pred_logR)))

df<-df[sort(order(df$S)),]

head(df)


p <- ggplot(df)
#p <- p + geom_point(aes(x=S_adj,y=R))
p <- p + geom_line(aes(x=S,y=predR, color=ayr), size=2, alpha=0.6)
p <- p + geom_text(aes(x=S,y=R,label=BroodYear ),hjust=0, vjust=0)
p <- p + theme_bw(16)
p <- p + labs(title = "Ricker model with explicit survival", x = "Spawners", y = "Recruits", color = "Year\n") 
#ylab("Recruits") + xlab("Spawners")
p
setwd(figs_dir)
ggsave("survival_model_fit.pdf", plot=p, width=10,height=7)



#MCMC
survivalB<-list(
  obj=survivalobj,
  nchain=3,
  iter=1000000,
  lowbd=c(0.1,-13.7,-6.0),
  hibd=c(8.0,-8.5,5.0)
)

posterior_survival<-posteriorsdf(survivalB)

#plots
plot_posteriors(posterior_survival$posteriors)

#posteriors of derived quantities -- the interesting ones
survdf<-posterior_survival$posteriors

a<-(survdf$value[survdf$parameters=="alpha"])
alpha<-exp(survdf$value[survdf$parameters=="alpha"])
beta<-exp(survdf$value[survdf$parameters=="logbeta"])
Smax<-1/exp(survdf$value[survdf$parameters=="logbeta"])
sig<-exp(survdf$value[survdf$parameters=="logSigObs"])

deriv_posteriors<-data.frame(chains=rep(survdf$chains[survdf$parameters=="logbeta"],4),
                             parameters = rep(c("a","b","Smax","sig"),each=length(a)),
                             value = c(a,beta,Smax,sig)
                             )

plot_posteriors(deriv_posteriors,salvar=TRUE,DIR=figs_dir,nome="posterior_survival_model_fit.pdf")


Ds<-list(
  DIR=tex_dir,
  param_names=c("a","b","$\\alpha$","$S_{max}$"),
  MLE=c(repsurv$alpha,repsurv$beta,exp(repsurv$alpha),1/repsurv$beta),
  MCMC=cbind(a,beta,alpha,Smax),
  caption = "Parameter estimates for Ricker model with survival as a covariate.",
  digits=matrix(c(2,-2,2,2),ncol=6,nrow=4),
  other=NULL,
  filename="surv_tab.tex")

 results_table(Ds) 





ps <- ggplot(SR)
ps <- ps + geom_point(aes(x=S_adj,y=Age2Surv, col=as.factor(BroodYear)),size=2)
ps <- ps + geom_line(aes(x=S_adj,y=Age2Surv))
ps

py <- ggplot(SR)
py<- py + geom_point(aes(x=BroodYear,y=Age2Surv))
py<- py + geom_line(aes(x=BroodYear,y=Age2Surv))
py


###################
#MCMC
###################


fitmcmc2 <- tmbstan(objkf, chains=3,
              iter=100000, init="random",
              lower=c(0.1,-13.5,0.0,-3.0), upper=c(4.0,-9.,1.0,5.0),
               control = list(adapt_delta = 0.99))

mc <- extract(fitmcmc2, pars=names(objkf$par),
              inc_warmup=TRUE, permuted=FALSE)
npar<-dim(mc)[3]

rhats <- stan_rhat(fitmcmc2,par=c(names(objkf$par),"alpha"))
stan_ess(fitmcmc2,par=c(names(objkf$par),"alpha"))
stan_diag(fitmcmc2, 
            information = 'sample', )



fit_summary <- summary(fitmcmc2)
fit_summary$summary[5:34,6]


posterior <- as.array(fitmcmc2)
mainrun<-melt(posterior)
head(mainrun)
summary(mainrun)
mainrun$type="posdata"
mainrun$type<-as.factor(mainrun$type)


pm<-ggplot(mainrun)
pm<-pm+geom_density(aes(x=value, color=chains))
pm<-pm+facet_wrap(~parameters, scales="free")
pm

unique(mainrun$parameters)


a<-(mainrun$value[mainrun$parameters=="alphao"])
beta<-exp(mainrun$value[mainrun$parameters=="logbeta"])
Smax<-1/exp(mainrun$value[mainrun$parameters=="logbeta"])
rho<-(mainrun$value[mainrun$parameters=="rho"])


deriv_posteriors<-data.frame(chains=rep(mainrun$chains[mainrun$parameters=="logbeta"],4),
                             parameters = rep(c("ao","b","Smax","rho"),each=length(a)),
                             value = c(a,beta,Smax,rho)
                             )

pm<-ggplot(deriv_posteriors)
pm<-pm+geom_density(aes(x=value, color=chains))
pm<-pm+facet_wrap(~parameters, scales="free")
pm
setwd(figs_dir)
ggsave("posterior_recursive_model.pdf", plot=pm, width=10,height=7)

#prior on rho
pdf("prior_rho.pdf")
plot(seq(0,1,by=.05),dbeta(seq(0,1,by=.05),3,3), type="l",
  xlab=expression(paste(rho," value")), ylab= "density", las=1, lwd=3)
dev.off()




recursive_tab<-data.frame(Parameter=c(paste("a",SR$BroodYear),"b","$S_{max}$","$\\rho$"),
                      MLE=c(repkf$alpha,repkf$beta,repkf$Smax,repkf$rho),
                      Median=c(fit_summary$summary[5:34,6],quantile(beta,.5),quantile(Smax,.5),quantile(rho,.5)),
                      Lower=c(fit_summary$summary[5:34,4],quantile(beta,.025),quantile(Smax,.025),quantile(rho,.025)),
                      Upper=c(fit_summary$summary[5:34,8],quantile(beta,.975),quantile(Smax,.975),quantile(rho,.975)))

setwd(tex_dir)
print(xtable(recursive_tab, digits=2,caption = "Parameter estimates for recursive Bayes Ricker model.")
  ,sanitize.text.function = function(x) {x},include.rownames = FALSE, 
  file="recursive_tab.tex")





dfa<-data.frame(broodyear=SR$BroodYear,a=c(repkf$alpha,fit_summary$summary[5:34,6]), 
  type=rep(c("MLE","Bayes"),each=length(SR$BroodYear)),lower=c(rep(NA,length(SR$BroodYear)),fit_summary$summary[5:34,4]),
  upper=c(rep(NA,length(SR$BroodYear)),fit_summary$summary[5:34,8]))


pa<-ggplot(dfa)
pa<-pa+geom_line(aes(x=broodyear,y=a,col=type), size=2)
pa<-pa+geom_ribbon(aes(x=broodyear,ymin=lower,ymax=upper, fill=type),alpha=.4)
pa<-pa+theme_bw(16)
pa<-pa+labs(title = "Recursive Bayes model - time series", y = expression(a), x = "Brood year") 
pa
setwd(figs_dir)
ggsave("recursive_a.pdf", plot=pa, width=10,height=7)



























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


