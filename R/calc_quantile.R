#==============================================
#Title:calc_quantile
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to calculate 95% quantiles for plotting
#==============================================


calc_quantile<-function(x,ci){

 z=quantile(x,probs = c(.025,0.05,.5,.95,.975),na.rm = TRUE)

 return(z)

}