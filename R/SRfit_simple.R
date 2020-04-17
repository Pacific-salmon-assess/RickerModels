#=================================================================
#Simple fit exploration of SR data
#Author: Catarina Wor
#Date: April 10th 2018
#=================================================================




#load in required packages and libraries
library(ggplot2)

#load in directories.R
source("C:/Users/worc/Documents/HarrisonSR/R/directories.R")

#read in simple data set
setwd(data_dir)
SR<-read.csv("Harrison_simples_Apr18.csv")

summary(SR)

#remove year 2004 because there is no estimate of survival for that year. 
Si <- SR$S_adj[!is.na(SR$Age2Surv)]
Ri <- SR$R[!is.na(SR$Age2Surv)]
log_RS <- log(SR$R/Si)[!is.na(SR$Age2Surv)]
Surv<-SR$Age2Surv[!is.na(SR$Age2Surv)]
PHatch<-SR$Prop_Hatchery[!is.na(SR$Age2Surv)]
byr<-SR$BroodYear[!is.na(SR$Age2Surv)]

plot(byr,Surv, type="b",lwd=2)
##### Ricker S/R model: log(R/S)=a+bS #####

lm(log(SR$R/SR$S_adj)~ SR$S_adj)

mod1<-lm(log_RS~Si)
mod2<-lm(log_RS~Si+Surv)
mod2l<-lm(log_RS~Si+log(Surv))

mod3<-lm(log_RS~Si+PHatch)
mod4<-lm(log_RS~Si+Surv+PHatch)
mod4l<-lm(log_RS~Si+log(Surv)+PHatch)


modaic<-AIC(mod1,mod2,mod2l,mod3,mod4,mod4l)

deltaAIC<-modaic$AIC-min(modaic$AIC)


mod1_a<-mod1$coefficients[1]
mod1_b<-mod1$coefficients[2]

mod2_a<-mod2$coefficients[1]
mod2_b<-mod2$coefficients[2]
mod2_b1<-mod2l$coefficients[3]

mod3_a<-mod3$coefficients[1]
mod3_b<-mod3$coefficients[2]
mod3_b1<-mod3$coefficients[3]


mod4_a<-mod4$coefficients[1]
mod4_b<-mod4$coefficients[2]
mod4_b1<-mod4$coefficients[3]
mod4_b2<-mod4$coefficients[4]


predR1<- Si[!is.na(Surv)]*exp(mod1_a+mod1_b*Si[!is.na(Surv)])
predR2<- Si[!is.na(Surv)]*exp(mod2_a+mod2_b*Si[!is.na(Surv)]+mod2_b1*Surv[!is.na(Surv)])
predR3<- Si[!is.na(Surv)]*exp(mod3_a+mod3_b*Si[!is.na(Surv)]+mod3_b1*PHatch[!is.na(Surv)])
predR4<- Si[!is.na(Surv)]*exp(mod4_a+mod4_b*Si[!is.na(Surv)]+mod4_b1*Surv[!is.na(Surv)]+mod4_b2*PHatch[!is.na(Surv)])
pred<-c(predR1,predR2,predR3,predR4)
model<-as.factor(rep(1:4, each=length(Si[!is.na(Surv)])))

SR<-SR[!is.na(Surv),]

dim(SR)
length(predR1)
df<-cbind(rbind(SR,SR,SR,SR),pred,model)
df<-df[order(sort(df$S_adj)),]

p <- ggplot(df)
#p <- p + geom_point(aes(x=S_adj,y=R))
p <- p + geom_text(aes(x=S_adj,y=R,label=BroodYear ),hjust=0, vjust=0)
p <- p + geom_line(aes(x=S_adj,y=pred, color=model), size=1.5)
p <- p + theme_bw(16)
p <- p + ylab("Recruits") + xlab("Spawners")

p




plot(SR$BroodYear,SR$Age2Surv, type="b", lwd=2)





















