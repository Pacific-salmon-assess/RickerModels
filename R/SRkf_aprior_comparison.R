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


#setwd(model_dir)
source("calc_quantile.R")
source("TMB_functions.R")

#read in simple data set
SR<-read.csv("../data/Harrison_simples_Apr18.csv")

iteracs=100000
# LM version
#simple model

srm<-lm(log(SR$R/SR$S_adj)~ SR$S_adj)
a_srm<-srm$coefficients[1]
b_srm<--srm$coefficients[2]
alpha<-exp(a_srm)

u_msy=.5*a_srm-0.07*a_srm^2
predR1<- SR$S_adj*exp(a_srm-b_srm*SR$S_adj)

mydata<-list(obs_logR=log(SR$R),obs_S=SR$S_adj)
parameters_simple <- list(
  alpha=(a_srm),
  logbeta = log(b_srm),
  logSigObs= log(.4)  )

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
  iter=iteracs,
  lowbd=c(0.1,-13.5,-6.0),
  hibd=c(4.0,-8.5,5.0)
)

posterior_simple<-posteriorsdf(simpleB)

#plots
#plot_posteriors(posterior_simple$posteriors)
#posteriors of derived quantities -- the interesting ones
simpdf<-posterior_simple$posteriors

a<-(simpdf$value[simpdf$parameters=="alpha"])
alpha<-exp(simpdf$value[simpdf$parameters=="alpha"])
beta<-exp(simpdf$value[simpdf$parameters=="logbeta"])
Smax<-1/exp(simpdf$value[simpdf$parameters=="logbeta"])
sig<-exp(simpdf$value[simpdf$parameters=="logSigObs"])
umsy_simple<-.5*a-0.07*a^2





Rprdbsimple<-matrix(NA,ncol=length(SR$BroodYear),nrow=length(sig))

meanRsimple<-NULL

varRsimple<-NULL


for(i in 1:length(sig)){

  Rprdbsimple[i,]<-as.numeric(SR$S_adj*exp(a[i])*exp(-SR$S_adj* beta[i]))
  
  meanRsimple[i]<-mean(Rprdbsimple[i,]/SR$R) 

  varRsimple[i]<-var(Rprdbsimple[i,]/SR$R)

  
}


#=============================================================================================================
#Recursive Bayes model - testing the scenario with low observation error only


mydata3<-list(obs_logR=log(SR$R),obs_S=SR$S_adj,prbeta1=2.0,prbeta2=3.0)


parameters_recursive <- list(
  alphao=a_srm,
  logbeta = log(b_srm),
  rho=.2,
  logvarphi= 0.1,
  alpha=rep(0.9,length(SR$R))
  )


recursive3<-list(
  dat=mydata3,
  params=parameters_recursive,
  rndm="alpha",
  dll="Rickerkf_ratiovar",
  DIR=model_dir
  )

compile("Rickerkf_ratiovar.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="")
dyn.load(dynlib("Rickerkf_ratiovar"))
obj3<-MakeADFun(mydata3,parameters_recursive,DLL="Rickerkf_ratiovar")
newtonOption(obj3, smartsearch=FALSE)

opt<-nlminb(obj3$par,obj3$fn,obj3$gr)
repkf3<-obj3$report()
sdreport(obj3)

recursiveB3<-list(
  obj=obj3,
  nchain=3,
  iter=iteracs,
  lowbd=c(0.1,-13.5,0.0,-3.0,rep(0.001,length(SR$R))),
  hibd=c(4.0,-9.,1.0,5.0,rep(4.0,length(SR$R)))
)

posterior_recursive3<-posteriorsdf(recursiveB3)
recrsdf3<-posterior_recursive3$posteriors

a_rb3<-(recrsdf3$value[recrsdf3$parameters=="alphao"])
beta_rb3<-exp(recrsdf3$value[recrsdf3$parameters=="logbeta"])
Smax_rb3<-1/exp(recrsdf3$value[recrsdf3$parameters=="logbeta"])
rho_rb3<-(recrsdf3$value[recrsdf3$parameters=="rho"])


soa3<-recrsdf3[grep("alpha",recrsdf3$parameters),]
soa3$umsy<-.5*soa3$value-0.07*soa3$value^2
umsyposteriorsummary3<-aggregate(soa3$umsy,list(soa3$parameters),function(x){quantile(x,probs=c(0.025,.5,.975))})
umsyposteriorsummary3<-umsyposteriorsummary3[c(1,12,23,25:30,2:11,13:22,24),]

#====================================================================================================================
#model with posfun on alpha



compile("Rickerkf_ratiovar_aprior.cpp",libtmb=FALSE, "-O1 -g", DLLFLAGS="")
dyn.load(dynlib("Rickerkf_ratiovar_aprior"))
obj3_aprior<-MakeADFun(mydata3,parameters_recursive,random="alpha",DLL="Rickerkf_ratiovar_aprior")
newtonOption(obj3_aprior, smartsearch=FALSE)

opt<-nlminb(obj3_aprior$par,obj3_aprior$fn,obj3_aprior$gr)
repkf3_aprior<-obj3_aprior$report()



recursiveB3_aprior<-list(
  obj=obj3_aprior,
  nchain=3,
  iter=iteracs,
  lowbd=c(0.1,-13.5,0.0,-3.0),
  hibd=c(4.0,-9.,1.0,5.0)
)

posterior_recursive3_aprior<-posteriorsdf(recursiveB3_aprior)



#==============================================================================
#plots



dfac<-data.frame(broodyear=rep(SR$BroodYear,8),
  a=c(posterior_recursive3$fit_summary$summary[5:34,6],
    posterior_recursive3_aprior$fit_summary$summary[5:34,6]
    ), #,umsyposteriorsummary3$x[,2]
  scn=rep(c("no bounds","bounds"),each=length(SR$BroodYear)),
  lower=c(posterior_recursive3$fit_summary$summary[5:34,4],
    posterior_recursive3_aprior$fit_summary$summary[5:34,4]),
  upper=c(posterior_recursive3$fit_summary$summary[5:34,8],
    posterior_recursive3_aprior$fit_summary$summary[5:34,8]))


pa<-ggplot(dfac)
pa<-pa+geom_ribbon(aes(x=broodyear,ymin=lower,ymax=upper, fill=scn),alpha=.4)
pa<-pa+geom_line(aes(x=broodyear,y=a,col=scn), size=1.2)
pa<-pa+theme_bw(16)
pa<-pa+labs(title = "Recursive Bayes model - a time series", y = expression(a[t]), x = "Brood year") 
pa






















moltentmp<-melt(soa[,-c(1,2)])
moltentmp$iterations<-paste(soa$iterations,soa$chains,sep=".")
tmp<-cast(moltentmp[moltentmp$variable=="value",],formula =  iterations~ parameters , mean, value = 'value')
tmpu<-cast(moltentmp[moltentmp$variable=="umsy",],formula =  iterations~ parameters , mean, value = 'value')

moltentmp1<-melt(soa1[,-c(1,2)])
moltentmp1$iterations<-paste(soa1$iterations,soa1$chains,sep=".")
tmp1<-cast(moltentmp1[moltentmp1$variable=="value",],formula =  iterations~ parameters , mean, value = 'value')
tmpu1<-cast(moltentmp1[moltentmp1$variable=="umsy",],formula =  iterations~ parameters , mean, value = 'value')

moltentmp2<-melt(soa2[,-c(1,2)])
moltentmp2$iterations<-paste(soa2$iterations,soa2$chains,sep=".")
tmp2<-cast(moltentmp2[moltentmp2$variable=="value",],formula =  iterations~ parameters , mean, value = 'value')
tmpu2<-cast(moltentmp2[moltentmp2$variable=="umsy",],formula =  iterations~ parameters , mean, value = 'value')

moltentmp3<-melt(soa3[,-c(1,2)])
moltentmp3$iterations<-paste(soa3$iterations,soa3$chains,sep=".")
tmp3<-cast(moltentmp3[moltentmp3$variable=="value",],formula =  iterations~ parameters , mean, value = 'value')
tmpu3<-cast(moltentmp3[moltentmp3$variable=="umsy",],formula =  iterations~ parameters , mean, value = 'value')
sum(moltentmp$value[moltentmp$variable=="umsy"]>0)
summary(moltentmp)
dim(moltentmp[moltentmp$variable=="umsy",])

dim(tmpu)
head(tmp[,c(2,13,24,26:31,3:12,14:23,25)])
head(tmpu[,c(2,13,24,26:31,3:12,14:23,25)])

Rprdb0<-matrix(NA,ncol=length(SR$BroodYear),nrow=nrow(tmp))
Rprdb1<-matrix(NA,ncol=length(SR$BroodYear),nrow=nrow(tmp))
Rprdb2<-matrix(NA,ncol=length(SR$BroodYear),nrow=nrow(tmp))
Rprdb3<-matrix(NA,ncol=length(SR$BroodYear),nrow=nrow(tmp))

meanR0<-NULL
meanR1<-NULL
meanR2<-NULL
meanR3<-NULL

varR0<-NULL
varR1<-NULL
varR2<-NULL
varR3<-NULL

avga0<-NULL
avga1<-NULL
avga2<-NULL
avga3<-NULL

avgu0<-NULL
avgu1<-NULL
avgu2<-NULL
avgu3<-NULL

i=1

for(i in 1:nrow(tmp)){#

  Rprdb0[i,]<-as.numeric(SR$S_adj*exp(tmp[i,c(2,13,24,26:31,3:12,14:23,25)])*exp(-SR$S_adj* beta_rb[i]))
  Rprdb1[i,]<-as.numeric(SR$S_adj*exp(tmp1[i,c(2,13,24,26:31,3:12,14:23,25)])*exp(-SR$S_adj* beta_rb1[i]))
  Rprdb2[i,]<-as.numeric(SR$S_adj*exp(tmp2[i,c(2,13,24,26:31,3:12,14:23,25)])*exp(-SR$S_adj* beta_rb2[i]))
  Rprdb3[i,]<-as.numeric(SR$S_adj*exp(tmp3[i,c(2,13,24,26:31,3:12,14:23,25)])*exp(-SR$S_adj* beta_rb3[i]))

  meanR0[i]<-mean(Rprdb0[i,]/SR$R)
  meanR1[i]<-mean(Rprdb1[i,]/SR$R)
  meanR2[i]<-mean(Rprdb2[i,]/SR$R)
  meanR3[i]<-mean(Rprdb3[i,]/SR$R)

  varR0[i]<-var(Rprdb0[i,]/SR$R)
  varR1[i]<-var(Rprdb1[i,]/SR$R)
  varR2[i]<-var(Rprdb2[i,]/SR$R)
  varR3[i]<-var(Rprdb3[i,]/SR$R)

  avga0[i]<-mean(exp(unlist(tmp[i,c(21:23,25)])))
  avga1[i]<-mean(exp(unlist(tmp1[i,c(21:23,25)])))
  avga2[i]<-mean(exp(unlist(tmp2[i,c(21:23,25)])))
  avga3[i]<-mean(exp(unlist(tmp3[i,c(21:23,25)])))

  avgu0[i]<-mean(unlist(tmpu[i,c(21:23,25)]))
  avgu1[i]<-mean(unlist(tmpu1[i,c(21:23,25)]))
  avgu2[i]<-mean(unlist(tmpu2[i,c(21:23,25)]))
  avgu3[i]<-mean(unlist(tmpu3[i,c(21:23,25)]))
  
}

meanu<-c(quantile(umsy_simple,.5),quantile(avgu0,.5),quantile(avgu2,.5),quantile(avgu3,.5)) #
meana<-c(quantile(alpha,.5),quantile(avga0,.5),quantile(avga2,.5),quantile(avga3,.5)) #
meanR<-c(quantile(meanRsimple,.5),quantile(meanR0,.5),quantile(meanR2,.5),quantile(meanR3,.5))#quantile(meanR1,.5),
varR<-c(quantile(varRsimple,.5),quantile(varR0,.5),quantile(varR2,.5),quantile(varR3,.5)) #,quantile(varR1,.5)
bs<-c(quantile(Smax,.5),quantile(Smax_rb,.5),quantile(Smax_rb2,.5),quantile(Smax_rb3,.5))



vraptab<- data.frame(scenario=rep(c("simple Ricker", "RB base 3, 3","RB high obs 3, 2","RB low obs 2, 3"),1),aavg4=meana,
  b=bs,meanR=meanR,varR=varR,uavg4=meanu)

colnames(vraptab)<-c("scenario","$\\widetilde{a_{avg4}}$","$\\widetilde{b}$","meanR","varR","$U_{MSY avg}$")
 vraptabtmp<-xtable(vraptab, digits=matrix(c(2,rep(2,nrow(vraptab)-1)),ncol=ncol(vraptab)+1,nrow=nrow(vraptab),byrow=F),sanitize.text.function=function(x){x},
,caption = "Median parameter estimates across simple Ricker recruitment model and 
 three time-varying productivity scenarios (RB) with different $\\rho$ priors. Parameters were transformed to meet VRAP requirements.",
 label="vraptab" )


print(vraptabtmp,sanitize.text.function = function(x) {x},
      include.rownames = FALSE, 
file="../tex/vrap_params2.tex",caption.placement = "top", label="vraptab")



#========================================================================
#Comparison tables
#parameter estimates

tab<-data.frame(Parameter=c("b","$S_{max}$","$\\rho$",paste("$\\alpha$",SR$BroodYear)),
                      "base" =c(apply(cbind(beta_rb,Smax_rb,rho_rb),2,function(x) quantile(x, .5)),exp(posterior_recursive$fit_summary$summary[5:34,6])),
                      "uninformative"=c(apply(cbind(beta_rb1,Smax_rb1,rho_rb1),2,function(x) quantile(x, .5)),exp(posterior_recursive1$fit_summary$summary[5:34,6])),
                      "low"=c(apply(cbind(beta_rb2,Smax_rb2,rho_rb2),2,function(x) quantile(x, .5)),exp(posterior_recursive2$fit_summary$summary[5:34,6])),
                      "high"=c(apply(cbind(beta_rb3,Smax_rb3,rho_rb3),2,function(x) quantile(x, .5)),exp(posterior_recursive3$fit_summary$summary[5:34,6])))


 tabtmp<-xtable(tab, digits=matrix(c(-2,rep(2,nrow(tab)-1)),ncol=ncol(tab)+1,nrow=nrow(tab),byrow=F)
,caption = "Median parameter estimates for the recursive Bayes model across four $\\rho$ prior scenarios. " )

  print(tabtmp,sanitize.text.function = function(x) {x},
      include.rownames = FALSE, 
  file="../tex/recursive_tab_params_comparison.tex",caption.placement = "top",
  label="parreccomp")


#rho estimates, priors, posteriors
#table for vrap





#========================================================================
#Comparison figures
#MLE versus posteriors


rhodfp<-rbind(rhodf0p,rhodf1p,rhodf2p,rhodf3p)
names(rhodfp)

linerhodf<-data.frame(scn=rep(c("uninformative 1, 1","base 3, 3","low obs 2, 3","high obs 3, 2"),each=2),
  vals=c(quantile(rho_rb1,probs=c(.5)),repkf1$rho,
    quantile(rho_rb,probs=c(.5)),repkf$rho,
    quantile(rho_rb3,probs=c(.5)),repkf3$rho, 
    quantile(rho_rb2,probs=c(.5)),repkf2$rho),
  type=as.factor(rep(c("Median posterior","MLE"),4)))

rhodf3p$median<-quantile(rho_rb3,probs=c(.5))
rhodf3p$mle<-repkf3$rho


pr<-ggplot(rhodfp,aes(value))
pr<-pr+geom_density(size=1.2,aes(fill=distribution),alpha=.6)
pr<-pr+facet_wrap(~ scn)
pr<-pr+geom_vline(data=linerhodf,aes(xintercept=vals,color=type),size=1.3)
pr<-pr+ scale_color_manual(values=c("blue4","red"))
#+ scale_color_manual(values=c("red","black","#E69F00", "#56B4E9"))scale_fill_brewer(palette="Dark2")+
pr<-pr+theme_bw(16)
pr<-pr+labs(x = expression(paste(rho,"values"))) 
pr
ggsave("../figs/priors_posteriors_rho.pdf", plot=pr, width=10,height=8)



dfu<-data.frame(broodyear=rep(SR$BroodYear,8),
  umsy=c(umsyposteriorsummary$x[,2],umsyposteriorsummary1$x[,2],
    umsyposteriorsummary2$x[,2],umsyposteriorsummary3$x[,2],
    repkf$umsy,repkf1$umsy,repkf2$umsy,repkf3$umsy
    ), #,umsyposteriorsummary3$x[,2]
  scn=rep(rep(c("base 3, 3","uninformative 1, 1","low obs 3, 2","high obs 2, 3"),each=length(SR$BroodYear)),2),#,"2,3"
  lower=c(umsyposteriorsummary$x[,1],umsyposteriorsummary1$x[,1],
    umsyposteriorsummary2$x[,1],umsyposteriorsummary3$x[,1],
    rep(NA,each=length(SR$BroodYear)*4)),
  upper=c(umsyposteriorsummary$x[,3],umsyposteriorsummary1$x[,3],
    umsyposteriorsummary2$x[,3],umsyposteriorsummary3$x[,3],
    rep(NA,each=length(SR$BroodYear)*4)),
  type=rep(c("Posterior", "MLE"),each=length(SR$BroodYear)*4)
  )




pu<-ggplot(dfu)
pu<-pu+geom_ribbon(aes(x=broodyear,ymin=lower,ymax=upper, fill=scn),alpha=.4)
pu<-pu+geom_line(aes(x=broodyear,y=umsy,col=scn), size=1.2)
pu<-pu+ facet_wrap(~type)
pu<-pu+theme_bw(16)
pu<-pu+labs(title = "Recursive Bayes model -  Umsy time series", y = expression(u[MSY]), x = "Brood year") 
pu





dfa<-data.frame(broodyear=rep(SR$BroodYear,8),
  a=c(posterior_recursive$fit_summary$summary[5:34,6],
    posterior_recursive1$fit_summary$summary[5:34,6],
    posterior_recursive2$fit_summary$summary[5:34,6],
    posterior_recursive3$fit_summary$summary[5:34,6],
    repkf$alpha,repkf1$alpha,repkf2$alpha,repkf3$alpha
    ), #,umsyposteriorsummary3$x[,2]
  scn=rep(rep(c("base 3, 3","uninformative 1, 1","low obs 3, 2","high obs 2, 3"),each=length(SR$BroodYear)),2),#,"2,3"
  lower=c(posterior_recursive$fit_summary$summary[5:34,4],
    posterior_recursive1$fit_summary$summary[5:34,4],
    posterior_recursive2$fit_summary$summary[5:34,4],
    posterior_recursive3$fit_summary$summary[5:34,4],
    rep(NA,each=length(SR$BroodYear)*4)),
  upper=c(posterior_recursive$fit_summary$summary[5:34,8],
    posterior_recursive1$fit_summary$summary[5:34,8],
    posterior_recursive2$fit_summary$summary[5:34,8],
    posterior_recursive3$fit_summary$summary[5:34,8],
    rep(NA,each=length(SR$BroodYear)*4)),
  type=rep(c("Posterior", "MLE"),each=length(SR$BroodYear)*4)
  )


pa<-ggplot(dfu)
pa<-pa+geom_ribbon(aes(x=broodyear,ymin=lower,ymax=upper, fill=scn),alpha=.4)
pa<-pa+geom_line(aes(x=broodyear,y=umsy,col=scn), size=1.2)
pa<-pa+ facet_wrap(~type)
pa<-pa+theme_bw(16)
pa<-pa+labs(title = "Recursive Bayes model - a time series", y = expression(a[t]), x = "Brood year") 
pa
ggsave("../recursive_a_priors.pdf", plot=pa, width=10,height=7)








