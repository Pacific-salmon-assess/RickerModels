#=======================================================
# Calculate the avg and var quantities used to transform recursive bayes into 
#VRAP parameters
#



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

avgresult<-list(Rprdb0=Rprdb0,
Rprdb1 = Rprdb1,
Rprdb2 = Rprdb2, 
Rprdb3 = Rprdb3,
meanR0 = meanR0,
meanR1 = meanR1, 
meanR2 = meanR2,
meanR3 = meanR3,
varR0 = varR0,
varR1 = varR1,
varR2 = varR2,
varR3 = varR3,
avga0 = avga0,
avga1 = avga1,
avga2 = avga2,
avga3 = avga3,
avgu0 = avgu0,
avgu1 = avgu1,
avgu2 = avgu2,
avgu3 = avgu3

)
saveRDS(avgresult, file = "../data/vrap_RB_results.rds")
