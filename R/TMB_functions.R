#==============================================================
#Functions to run TMB stock recruitment functions
#Author: Catarina Wor
#Date: 14th of June 2018
#==============================================================


runTMB<-function(A,comps=TRUE){

	####
	#This function compiles and runs a standard TMB model
	###
  setwd(A$DIR)

  cppfile<-paste(A$dll,".cpp",sep="")

  if(comps==TRUE){
    compile(cppfile,libtmb=FALSE, "-O1 -g", DLLFLAGS="")
    dyn.load(dynlib(A$dll))
  }
  
  obj<-MakeADFun(A$dat,A$params,random=A$rndm,DLL=A$dll)
  newtonOption(obj, smartsearch=FALSE)
  
  opt<-nlminb(obj$par,obj$fn,obj$gr)
  

  return(list(obj=obj,
    opt=opt))
}




posteriorsdf<-function(B){

	####
	#This function runs TMB stan and produces clean posteriors in data-frame format
	###

    fitmcmc1 <- tmbstan(B$obj, chains=B$nchain,
              iter=B$iter, init="random",
              lower=B$lowbd, upper=B$hibd,
               control = list(adapt_delta = 0.98))

    mc <- extract(fitmcmc1, pars=names(B$obj$par),
              inc_warmup=TRUE, permuted=FALSE)
    

    fit_summary <- summary(fitmcmc1)

    posterior <- as.array(fitmcmc1)

    mainrun <- melt(posterior)

    poslist <- list(fit=fitmcmc1,
        fit_summary=fit_summary,
        posteriors=mainrun,
        mcmcobj=mc
      ) 

    return(poslist) 
}



plot_posteriors<-function(df,salvar=FALSE,DIR="",nome=""){

	####
	#This function plots posterior distributions by chain
	###

  pm<-ggplot(df)
  pm<-pm+geom_density(aes(x=value, color=chains))
  pm<-pm+facet_wrap(~parameters, scales="free")
  print(pm)

  if(salvar){
    setwd(DIR)
    ggsave(nome, plot=pm, width=10,height=7)

  }

}



results_table<-function(D){

	####
	#This function produces tex tables with SR function results
	####
  if(sum(!is.na(D$MCMC))){
    tab<-data.frame(Parameter=D$param_names,
                      MLE=D$MLE,
                      Median=c(apply(D$MCMC,2,function(x) quantile(x, .5)),D$other[2,]),
                      Lower=c(apply(D$MCMC,2,function(x) quantile(x, .025)),D$other[1,]),
                      Upper=c(apply(D$MCMC,2,function(x) quantile(x, .975)),D$other[3,]))
  }else{
    tab<-data.frame(Parameter=D$param_names,
                      MLE=D$MLE,
                      Median=c(D$other[2,]),
                      Lower=c(D$other[1,]),
                      Upper=c(D$other[3,]))
  }

    setwd(D$DIR)
    tabtmp<-xtable(tab, digits=D$digits,caption = D$caption,label=D$labs)
    digits(tabtmp)<-D$digits

    print(tabtmp,sanitize.text.function = function(x) {x},
      include.rownames = FALSE, 
  file=D$filename,caption.placement = "top")

}



model_pred_plot<-function(M, salvar=FALSE,DIR="",filename=""){

  sumfit<-apply(M$predBayes,2,function(x) quantile(x, probs=c(0.025,.5,0.975)))
  fitdf<-as.data.frame(t(sumfit))
  names(fitdf)<-c("lower","estimate","upper")
  fitdf<-cbind(fitdf,M$orig_data)
  fitdf$Type<-"Bayesian"

  fitdf1<-fitdf
  fitdf1$estimate<-M$predFreq
  fitdf1$lower<-NA
  fitdf1$upper<-NA
  fitdf1$Type<-"MLE"

  fitdf<-fitdf[order(SR$S_adj),]
  fitdf1<-fitdf1[order(SR$S_adj),]

  fitdf<-rbind(fitdf1,fitdf)

  p <- ggplot(fitdf)
  p <- p + geom_line(aes(x=S_adj,y=estimate, col=Type),size=1.5)
  p <- p + geom_ribbon(aes(x=S_adj,ymin=lower,ymax=upper, fill=Type),alpha=0.4)
  p <- p + geom_text(aes(x=S_adj,y=R,label=BroodYear ),hjust=0, vjust=0)
  p <- p + theme_bw(16) + scale_fill_viridis_d(end = 0.8,option = "B")
  p <- p + scale_color_viridis_d(end = 0.8,option = "B")
  p <- p + ylab("Recruits") + xlab("Spawners")
  print(p)

    if(salvar){
      setwd(DIR)
      ggsave(filename, plot=p, width=10,height=7)
    }

}
