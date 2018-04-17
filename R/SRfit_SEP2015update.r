##Script translated from S+ originally written by D. Chen
##for Harrison S-R analysis  (Brown et al. unpub.)
##December 2013

rm(list = ls(all=TRUE)); #Remove all the objects in the memory
root.dir <- getwd(); #root directory for the data
oldpar <- par(no.readonly=TRUE)

filename <- c("HRdata_84-96.csv")           ## CHANGE INPUT FILE NAME HERE
data.file <- paste(root.dir,"/",filename, sep="");
HRdat84.96 <- read.csv(file=data.file, header=TRUE, sep=",")

require(boot)

boot.ricker <- function(dat)
{

	st  <- dat$S
	rtt  <- dat$R
		
	log.recruit2spawner  <- log(rtt/st)
	
	ricker.m  <- lm(log.recruit2spawner~ st)
	s2  <- sum(ricker.m$residuals**2)/(length(st)-2)
    # The above value is actually s-squared or 'residual mean square error'
    # or as called in S+ output, the 'residual standard error'

    a  <- ricker.m$coef[1]
    b  <- -ricker.m$coef[2] # b coerced to a positive
    a.cor <- a+0.5*s2
    sopt  <- a.cor*(0.5-0.07*a.cor)/b
    uopt  <- sopt*b
    ropt  <- sopt*exp(a.cor-b*sopt)   # Use if b supplied as positive
#   ropt  <- sopt*exp(a.cor+b*sopt)   # Use if b supplied as negative

    sopt.tr <- getSoptR.t(a,b,s2)$x
    uopt.tr <- sopt.tr*b
    ropt.tr <- sopt.tr*exp(a.cor-b*sopt.tr)
#   ropt.tr <- sopt.tr*exp(a.cor+b*sopt.tr)

	c(a,b,sopt.tr,uopt.tr,sopt,uopt,ropt.tr,ropt,s2)

}


boot.SR <- function(dat, do.plot=F)
{
	st  <- dat$S
	rtt  <- dat$R
	sr  <- data.standard(dat$SR)      #untransformed, standardized

	log.recruit2spawner  <- log(rtt/st)
	
	lm.m  <- lm(log.recruit2spawner~ st+sr)
	s2 <- sum(lm.m$residuals^2)/(length(st)-3)
    # The above formula gets the 'biased-corrected' 
    # residual squared error (because of subtracting 3
    # which equals the number of parameters in the
    # regression equation

    a  <- lm.m$coef[1]
    b  <- -lm.m$coef[2]   # b coerced to a positive
    g  <- lm.m$coef[3]
    a.cor <- a+0.5*s2
    sopt  <- a.cor*(0.5-0.07*a.cor)/b
    uopt  <- sopt*b
    ropt  <- sopt*exp(a.cor-b*sopt) # Use if b supplied as positive
#   ropt  <- sopt*exp(a.cor+b*sopt)

    srparam <- g*mean(sr)
    sopt.tr <- getSoptSR.t(a,b,s2,srparam)$estimate
    uopt.tr <- sopt.tr*b
    ropt.tr <- sopt.tr*exp(a.cor-b*sopt.tr)
#   ropt.tr <- sopt.tr*exp(a.cor+b*sopt.tr)

	c(a,b,sopt.tr,uopt.tr,sopt,uopt,ropt.tr,ropt,s2,g)

}


boot.logSR <- function(dat)
{

	st  <- dat$S
	rtt  <- dat$R
	sr  <- data.standard(log(dat$SR))          #ln-transformed, standardized
		
	log.recruit2spawner  <- log(rtt/st)
	
	loglm.m  <- lm(log.recruit2spawner~ st+sr)
	s2 <- sum(loglm.m$residuals^2)/(length(st)-3)
    # The above formula gets the 'biased-corrected' residual squared error 
    # (by subtracting 3, which equals the number of parameters in the
    # regression equation)

    a  <- loglm.m$coef[1]
    b  <- -loglm.m$coef[2]  # b coerced to a positive
    g  <- loglm.m$coef[3]
    a.cor <- a+0.5*s2
    sopt  <- a.cor*(0.5-0.07*a.cor)/b
    uopt  <- sopt*b
    ropt  <- sopt*exp(a.cor-b*sopt) # Use if b supplied as positive
#   ropt  <- sopt*exp(a.cor+b*sopt) # Use if b supplied as negative

    srparam <- g*mean(sr)
    sopt.tr <- getSoptSR.t(a,b,s2,srparam)$minimum
    uopt.tr <- sopt.tr*b
    ropt.tr <- sopt.tr*exp(a.cor-b*sopt.tr)
#   ropt.tr <- sopt.tr*exp(a.cor+b*sopt.tr)

	c(a,b,sopt.tr,uopt.tr,sopt,uopt,ropt.tr,ropt,s2,g)

}

data.standard  <- function(x){
	(x-mean(x)) /sqrt(var(x))
}

Bootresids1 <- function(dat, i)
{
	# This function performs bootstrapping using the residuals
	# from a regression equation
   d <- dat
   d$y <- d$fitted + d$resids[i]
   lm.bs <- lm(d$y ~ d$s)
   a <- lm.bs$coef[1]
   b <- -lm.bs$coef[2]  # b coerced to a positive
   s <- sum(lm.bs$residuals^2)/lm.bs$df.residual
   a.cor <- a+0.5*s
   
   sopt.ha <- a.cor*(0.5-0.07*a.cor)/b
   uopt.ha <- sopt.ha*b
   ropt.ha <- sopt.ha*exp(a.cor-b*sopt.ha) # Use if b supplied as positive
#  ropt.ha <- sopt.ha*exp(a.cor+b*sopt.ha)

   sopt.tr <- getSoptR.t(a,b,s)$x
   uopt.tr <- sopt.tr*b
   ropt.tr <- sopt.tr*exp(a.cor-b*sopt.tr)
#  ropt.tr <- sopt.tr*exp(a.cor+b*sopt.tr)

   c(a,b,sopt.tr,uopt.tr,sopt.ha,uopt.ha,ropt.tr,ropt.ha,s)

} # End of the Bootresids1 function


Bootresids2 <- function(dat, i)
{
	# This function performs bootstrapping using the residuals
	# from a regression equation
   d <- dat
   d$y <- d$fitted + d$resids[i]
   lm.bs <- lm(d$y ~ d$s + d$sr)
   a <- lm.bs$coef[1]
   b <- -lm.bs$coef[2]  # b coerced to a positive
   g <- lm.bs$coef[3]
   s <- sum(lm.bs$residuals^2)/lm.bs$df.residual
   a.cor <- a+0.5*s
   sopt.ha <- a.cor*(0.5-0.07*a.cor)/b
   uopt.ha <- sopt.ha*b
   ropt.ha <- sopt.ha*exp(a.cor-b*sopt.ha) # Use if b supplied as positive
#  ropt.ha <- sopt.ha*exp(a.cor+b*sopt.ha)

   srparam <- g*mean(d$sr)
   sopt.tr <- getSoptSR.t(a,b,s,srparam)$estimate
   uopt.tr <- sopt.tr*b
   ropt.tr <- sopt.tr*exp(a.cor-b*sopt.tr)
#  ropt.tr <- sopt.tr*exp(a.cor+b*sopt.tr)

   # Add these new values to return string
   c(a,b,sopt.tr,uopt.tr,sopt.ha,uopt.ha,ropt.tr,ropt.ha,s,g)

} # End of the Bootresids2 function


getSoptR.t <- function(a, b, s)
{
# Assumes the Ricker a value is in natural logs
# Ricker b given as a positive value here

	co <- c(a,b,s)
	assign("co", co)
	func <- function(x)
	{
		abs( (1 - co[2]*x)*exp(co[1] - (co[2]*x) + (co[3]*0.5)) - 1)
	}
	optimize(f=func, a=co[2],bmax.iter = 20)
} # End of the getSoptR.t function


getSoptSR.t <- function(a,b,s,srparam)
{ 
 #function for calculating Smsy (sopt) based on transcendental equation.
 #i.e. find the value of sopt in a Ricker curve where the slope = 1.
 # Assumes the Ricker a value is in natural logs
 # Value for beta supplied as a positive
 
    co <- c(a,b,s,srparam)
    assign("co",co)
    func <- function(x)
    {
       abs( (1 - co[2]*x)*exp(co[1] + co[4] - (co[2]*x) + (co[3]*0.5)) - 1)
      # the statement below is equivalent to that above
      # abs(((1-(co[2]*x))*exp(co[1])*exp(co[4])*exp(-co[2]*x)*exp(co[3]*0.5))-1)
    }

    optimize(f=func, interval=c(0,sopt),maximum=F)

} # End of the getSoptSR.t function


LBesc <- function(a,b,y,p,m)
{ 
 # function for calculating spawners at a specified proportion of the
 # yield at Smsy (i.e. MSY).
 # Assumes the Ricker a value is in natural logs
 # Currently, the Ricker 'a' corrected with sigma-squared is being passed to this
 # function
 # Value for beta supplied as a positive
 # y = MSY
 # p = prop. of MSY
 # print(paste("Alpha = ", a, "  Beta = ", b, "  MSY = ", y, "  Prop. of MSY = ", p))
    
    co <- c(a,b,y,p,m)
    assign("co",co)
    func <- function(x)
    {
        if(co[5]==1) {
           # This function solves for a proportion of MSY
           abs(x*exp(co[1] - co[2]*x) - x - (co[3]*co[4])) 
           # The form below crashes
           # -(x*exp(co[1] - co[2]*x) - x - (co[3]*co[4])) 
           }
        else if(co[5]==2) {
           # This function solves for a proportion of production at MSY
           abs(x*exp(co[1] - co[2]*x) - (co[3]*co[4])) 
           }
    }

    nlmin(func, 1, print.level=0, max.iter=20)

} # End of the LBesc function


SumBSdat <- function(dat)
{
# This function gets summary data (means, etc) without any selection of data.
# The results are based on the entire data set, i.e., whatever is passed to the
# function.

SampN <- nrow(dat)

ci90  <-  quantile(dat[,1], c(0.050, 0.950))
med  <-  quantile(dat[,1], probs=0.500)
print(paste("Alpha:  Mean: ", round(mean(dat[,1]),5), "   Median: ", round(med,5), "   SD: ", round(stdev(dat[,1]),6), "   N: ", SampN ))
print(paste("Alpha:  90% CI: (", round(ci90[1],5),",", round(ci90[2],5), ")   Range: ", round(diff(ci90),5) ))

ci90  <-  quantile(dat[,2], c(0.050, 0.950))
med  <-  quantile(dat[,2], probs=0.500)
print(paste("Beta:   Mean: ", round(mean(dat[,2]),8), "   Median: ", round(med,8), "   SD: ", round(stdev(dat[,2]),9), "   N: ", SampN ))
print(paste("Beta:   90% CI: (", round(ci90[1],8),",", round(ci90[2],8), ")   Range: ", round(diff(ci90),8) ))

cat("\n")
print("Results below are based on solving for Smsy in the linear form of the Ricker function")
cat("\n")

ci90  <-  quantile(dat[,3], c(0.050, 0.950))
med  <-  quantile(dat[,3], probs=0.500)
print(paste("Smsy (transcendental):   Mean: ", round(mean(dat[,3]),0), "    Median: ", round(med,0), "   SD: ", round(stdev(dat[,3]),0), "   N: ", SampN ))
print(paste("Smsy (transcendental):   90% CI: (", round(ci90[1],0),",", round(ci90[2],0), ")   Range: ", round(diff(ci90),0) ))

ci90  <-  quantile(dat[,4], c(0.050, 0.950))
med  <-  quantile(dat[,4], probs=0.500)
print(paste("Umsy (transcendental):   Mean: ", round(mean(dat[,4]),4), "    Median: ", round(med,4), "  SD: ", round(stdev(dat[,4]),5), "   N: ", SampN ))
print(paste("Umsy (transcendental):   90% CI: (", round(ci90[1],4),",", round(ci90[2],4), ")   Range: ", round(diff(ci90),4) ))

ci90  <-  quantile(dat[,7], c(0.050, 0.950))
med  <-  quantile(dat[,7], probs=0.500)
print(paste("Rmsy (transcendental):   Mean: ", round(mean(dat[,7]),0), "    Median: ", round(med,0), "   SD: ", round(stdev(dat[,7]),0), "   N: ", SampN ))
print(paste("Rmsy (transcendental):   90% CI: (", round(ci90[1],0),",", round(ci90[2],0), ")   Range: ", round(diff(ci90),0) ))

cat("\n")
print("Results below are based on using the Hilborn Approximation equation to obtain Smsy in the linear form of the Ricker function")
cat("\n")

ci90  <-  quantile(dat[,5], c(0.050, 0.950))
med  <-  quantile(dat[,5], probs=0.500)
print(paste("Smsy (Hilborn approx.):  Mean: ", round(mean(dat[,5]),0), "   Median: ", round(med,0), "   SD: ", round(stdev(dat[,5]),0), "   N: ", SampN ))
print(paste("Smsy (Hilborn approx.):  90% CI: (", round(ci90[1],0),",", round(ci90[2],0), ")   Range: ", round(diff(ci90),0) ))

ci90  <-  quantile(dat[,6], c(0.050, 0.950))
med  <-  quantile(dat[,6], probs=0.500)
print(paste("Umsy (Hilborn approx.):  Mean: ", round(mean(dat[,6]),4), "   Median: ", round(med,4), "   SD: ", round(stdev(dat[,6]),5), "   N: ", SampN ))
print(paste("Umsy (Hilborn approx.):  90% CI: (", round(ci90[1],4),",", round(ci90[2],4), ")   Range: ", round(diff(ci90),4) ))

ci90  <-  quantile(dat[,8], c(0.050, 0.950))
med  <-  quantile(dat[,8], probs=0.500)
print(paste("Rmsy (Hilborn approx.):  Mean: ", round(mean(dat[,8]),0), "   Median: ", round(med,0), "   SD: ", round(stdev(dat[,8]),0), "   N: ", SampN ))
print(paste("Rmsy (Hilborn approx.):  90% CI: (", round(ci90[1],0),",", round(ci90[2],0), ")   Range: ", round(diff(ci90),0) ))

n <- ncol(dat)
if(n==13 || n==14) {
   ci90  <-  quantile(dat[,n-1], c(0.050, 0.950))
   med  <-  quantile(dat[,n-1], probs=0.500)
   print(paste("LBesc (transcendental):   Mean: ", round(mean(dat[,n-1]),0), "    Median: ", round(med,0), "   SD: ", round(stdev(dat[,n-1]),0) ))
   print(paste("LBesc (transcendental):   90% CI: (", round(ci90[1],0),",", round(ci90[2],0), ")   Range: ", round(diff(ci90),0) ))

   ci90  <-  quantile(dat[,n], c(0.050, 0.950))
   med  <-  quantile(dat[,n], probs=0.500)
   print(paste("LBesc (Hilborn approx.):  Mean: ", round(mean(dat[,n]),0), "   Median: ", round(med,0), "   SD: ", round(stdev(dat[,n]),0) ))
   print(paste("LBesc (Hilborn approx.):  90% CI: (", round(ci90[1],0),",", round(ci90[2],0), ")   Range: ", round(diff(ci90),0) ))
   }

# This function doesn't return any data so the following statement is necessary
invisible(NULL)

} # End of SumBSdat function


plot.BSout <- function(x, obs, tit)
{
# Histogram of resample values with mean and observed values marked.
# Includes smooth density estimate.
	
# Function to calculate bandwidth from Venables and Ripley p. 138.
	bandwidth.nrd <- function(x)
	{
		r <- quantile(x, c(0.25, 0.75))
		h <- (r[2] - r[1])/1.34
		4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-0.2)
	}
# Function to calculate nclass from Venables and Ripley p. 126.
	nclass.FD <- function(x)
	{
		r <- quantile(x, c(0.25, 0.75))
		names(r) <- NULL
		h <- 2 * (r[2] - r[1]) * length(x)^{-1/3}
		ceiling(diff(range(x))/h)
	}
	

		xi <- x[!is.na(x)]
		
		hist.vals <- hist(xi, nclass = nclass.FD(xi), 
			probability = T, plot = F, include.lowest = T)
		dens <- density(xi, width = bandwidth.nrd(xi), 
			from = hist.vals$breaks[1], 
			to   = hist.vals$breaks[length(hist.vals$breaks)], na.rm = T)
		ymax <- max(c(hist.vals$counts, dens$y))
		barplot(hist.vals$counts, width = hist.vals$breaks, histo = T, col=100, 
			ylim = c(0, ymax), xlab = tit, ylab = "Density")
		line.vals <- approx(dens, xout = c(obs, mean(x) ))
		lines(rep(line.vals$x[1], 2), c(0, line.vals$y[1]))
		lines(rep(line.vals$x[2], 2), c(0, line.vals$y[2]), lty = 2)
		lines(dens)

		### add the arrows

#		ci95  <-  quantile(xi, c(0.025, 0.975))
#print(paste(tit,":", " 95% CI=(", round(ci95[1],3),",", round(ci95[2],3), ") with Range=", round(diff(ci95),3) ,sep=" " ))
		ci90  <-  quantile(xi, c(0.050, 0.950))
#print(paste(tit,":", " 90% CI=(", round(ci90[1],3),",", round(ci90[2],3), ") with Range=", round(diff(ci90),3) ,sep=" " ))

		#arrows(mean(xi), .150, ci95[1], .150)
		#arrows(mean(xi), .150, ci95[2], .150)
		#arrows(mean(xi),  .070, ci95[1], .07)
		#arrows(mean(xi), .070, ci95[2], .07)

		std  <- sqrt(var(xi))
#print(paste("mean=",mean(xi), "std=",round(std,3),sep=" " ))
		#arrows(mean(xi), .06, mean(xi)-1.96*std, .06, open=T)
		#arrows(mean(xi), .06, mean(xi)+1.96*std, .06, open=T)

		
		box()
		#if(is.na(match("main", names(list(...)))))
		#	title(names(x$observed)[i])
	
	invisible(NULL)
}


data.standard  <- function(x){
	(x-mean(x)) /sqrt(var(x))
	}
 

SR.fit <- function(dat, num.boot, do.plot, ricker=F, smooth=F, climatic=F, logclimatic=F, bsresids=F, propMSY=1, lowerbound=F, method=2) 
{
	print("The following is the command line used to generate all the output from function SR.fit()")
	print(paste("Numboot =",num.boot, "DoPlot =",do.plot, "Ricker =",ricker, 
    "Climatic =",climatic, "LogClimatic =",logclimatic, "BSresids =",bsresids, "propMSY =",propMSY, sep=" "))
	cat("\n")

   # The line below attaches a special library of additional bootstrap functions
   # The function 'boot', used to bootstrap with residuals, is part of this library
   # The function 'bootstrap', used to bootstrap the observed data, is always attached
	if(bsresids) { library(boot) }

	### read in the data ###

	st  <- dat$S
	rt  <- dat$R
	sr  <- as.numeric(dat$SR) # SR values must be supplied as proportions
 	year<- dat$BY
	#er <- data.standard(dat$ER)
	log.recruit2spawner  <- log(rt/st)
	len.dat  <- length(st)

if (ricker){
	cat("\n")
	print("##### Ricker S/R model: log(R/S)=a+bS #####")
		
	ricker.m  <- lm(log.recruit2spawner~ st)
	cat("\n")
	print("Below are the results of the regression analysis of the linear form of the standard Ricker function")
	cat("\n")
	print(summary(ricker.m))

        par(cex=1.5)
	par(mfrow=c(2,2), las=0)

	plot(ricker.m$fitted.values,ricker.m$residuals, xlab="Fitted", ylab="Residuals", pch=16)
	abline(h=0)
	plot(year,ricker.m$residuals, xlab="Year", ylab="Residuals", pch=16)
	abline(h=0)
	plot(acf(ricker.m$residuals, type="correlation", plot=F), main="")
	plot(acf(ricker.m$residuals, type="partial", plot=F), main="")

	#### do cal ###

	a <- ricker.m$coef[1]
	b <- -ricker.m$coef[2]  # b coerced to a positive value
	s2.r <- sum(ricker.m$residuals**2)/(length(st)-2)
	a.cor <- a+0.5*s2.r

   # Smsy by Hilborn approximation method
	sopt.ha1  <-  a.cor*(0.5-0.07*a.cor)/b
	uopt.ha1  <-  a.cor*(0.5-0.07*a.cor)
   	ropt.ha1  <-  sopt.ha1*exp(a.cor-b*sopt.ha1) # Use when b supplied as positive
	rrep.ha1  <-  a.cor/b
#  	ropt.ha1  <-  sopt.ha1*exp(a.cor+b*sopt.ha1)
   	yopt.ha1  <-  ropt.ha1-sopt.ha1

# NOTE: CODE RELATED TO THE CALCULATION OF A LOWER BOUND FOR ESCAPEMENT
# HAS MAINLY BEEN COMMENTED OUT

   # 'Lower Bound' escapement based on a proportion of MSY
   # The proportion of MSY is supplied as a parameter in the command line
#   if(method==1) LB.s.ha1   <-  LBesc(a.cor,b,yopt.ha1,propMSY,method)$x
#   else if(method==2) LB.s.ha1   <-  LBesc(a.cor,b,ropt.ha1,propMSY,method)$x
#   LB.r.ha1   <-  LB.s.ha1*exp(a.cor-b*LB.s.ha1) # Use when b supplied as positive
##  LB.r.ha1   <-  LB.s.ha1*exp(a.cor+b*LB.s.ha1)
#   LB.y.ha1   <-  LB.r.ha1-LB.s.ha1#


	# Smsy by solving transcendental equation method
	sopt.tr1  <-  getSoptR.t(a,b,s2.r)$x
	uopt.tr1  <-  b*sopt.tr1
 	ropt.tr1  <-  sopt.tr1*exp(a.cor-b*sopt.tr1) # Use when b supplied as positive
	rrep.tr1  <-  a.cor/b
# 	ropt.tr1  <-  sopt.tr1*exp(a.cor+b*sopt.tr1)
   	yopt.tr1  <-  ropt.tr1-sopt.tr1

   # 'Lower Bound' escapement based on a % of MSY; the neg. sign for 'b'
   # to force it to a positive is necessary
#   if(method==1) LB.s.tr1   <-  LBesc(a.cor,b,yopt.tr1,propMSY,method)$x
#   else if(method==2) LB.s.tr1   <-  LBesc(a.cor,b,ropt.tr1,propMSY,method)$x
#   LB.r.tr1   <-  LB.s.tr1*exp(a.cor-b*LB.s.tr1) # Use when b supplied as positive
##  LB.r.tr1   <-  LB.s.tr1*exp(a.cor+b*LB.s.tr1) 
#   LB.y.tr1   <-  LB.r.tr1-LB.s.tr1

   cat("\n")
   print("Parameter estimates obtained from the standard Ricker stock-recruit function")
	print("Note: The calculation of Smsy includes the bias correction to alpha involving the") 
	print("residual mean square error from the regression results based on the linear form")
	print("of the Ricker function")
	cat("\n")
   print(paste("Alpha= ", round(a,4) ))
   print(paste("Beta= ", signif(b,4) ))
   print(paste("Sigma^2= ", round(s2.r,4) ))
   cat("\n")
   print("Estimates of optimal spawners obtained by the Hilborn approximation method")
	cat("\n")
   print(paste("Smsy : ", round(sopt.ha1,0), "     Umsy : ", round(uopt.ha1,4) ))
   print(paste("Rmsy : ", round(ropt.ha1,0), "    Ymsy : ", round(yopt.ha1,0) ))
	print(paste("Rrep : ", round(rrep.ha1,0) ))
#   print(paste("Lower Bound Escapement : ", round(LB.s.ha1,0), "     at ", propMSY*100, "% of MSY"))
#   print(paste("Lower Bound Production : ", round(LB.r.ha1,0), "    at ", propMSY*100, "% of MSY"))
#   print(paste("Lower Bound Yield      : ", round(LB.y.ha1,0), "    at ", propMSY*100, "% of MSY"))

   cat("\n")
   print("Estimates of optimal spawners obtained by solving the transcendental equation")
	cat("\n")
   print(paste("Smsy : ", round(sopt.tr1,0), "     Umsy : ", round(uopt.tr1,4)))
   print(paste("Rmsy : ", round(ropt.tr1,0), "    Ymsy : ", round(yopt.tr1,0) ))
	print(paste("Rrep : ", round(rrep.tr1,0) ))
#   print(paste("Lower Bound Escapement : ", round(LB.s.tr1,0), "    at ", propMSY*100, "% of MSY"))
#   print(paste("Lower Bound Production : ", round(LB.r.tr1,0), "   at ", propMSY*100, "% of MSY"))
#   print(paste("Lower Bound Yield      : ", round(LB.y.tr1,0), "   at ", propMSY*100, "% of MSY"))

	### do bootstrapping ##

 	if(num.boot!=0){
 	set.seed(100)
	cat("\n")
 	print("##### Doing Bootstrapping for Ricker model ###")
	if(bsresids) { # Do bootstrapping using residuals
		print("Performing bootstrapping using residuals")
    	SRdat <- data.frame(y=log.recruit2spawner, s=st, fitted=fitted(ricker.m), resids=resid(ricker.m))
       do.boot.ricker <- boot(SRdat, Bootresids1, R=num.boot)
       print("Summary of bootstrap results:")
		print("Rows indicate the following in sequence:")
		print("1. alpha; 2. beta; 3. Smsy(TR); 4. Umsy(TR); 5. Smsy(HA); 6. Umsy(HA); 7. Rmsy(TR); 8. Rmsy(HA); 9. sigma-squared")
       print.boot(do.boot.ricker)
       cat("\n")
#      plot.boot(do.boot.ricker)
       a.boot     <- do.boot.ricker$t[,1]
       b.boot     <- do.boot.ricker$t[,2]  # b is returned as a positive here
       sopt.tr.boot     <- do.boot.ricker$t[,3]
       uopt.tr.boot     <- do.boot.ricker$t[,4]
       sopt.ha.boot     <- do.boot.ricker$t[,5]
       uopt.ha.boot     <- do.boot.ricker$t[,6]
       ropt.tr.boot     <- do.boot.ricker$t[,7]
       ropt.ha.boot     <- do.boot.ricker$t[,8]
       s.boot     <- do.boot.ricker$t[,9]
       }
	else {  # Do bootstrapping based on original data set
		print("Performing bootstrapping using the observed data")
		print("Rows indicate the following in sequence:")
		print("1. alpha; 2. beta; 3. Smsy(TR); 4. Umsy(TR); 5. Smsy(HA); 6. Umsy(HA); 7. Rmsy(TR); 8. Rmsy(HA); 9. sigma-squared")
		cat("\n")
    	do.boot.ricker <- bootstrap(dat, boot.ricker(dat), B=num.boot)
       print("Summary of the results from bootstrapping:")
       print(do.boot.ricker)
       cat("\n")
#	    print("Now plotting do.boot.ricker")
#	    plot(do.boot.ricker)
    	
	    a.boot     <- do.boot.ricker$replicates[,1]
	    b.boot     <- do.boot.ricker$replicates[,2] # b is returned as a positive here
       sopt.tr.boot     <- do.boot.ricker$replicates[,3]
       uopt.tr.boot     <- do.boot.ricker$replicates[,4]
       sopt.ha.boot     <- do.boot.ricker$replicates[,5]
       uopt.ha.boot     <- do.boot.ricker$replicates[,6]
       ropt.tr.boot     <- do.boot.ricker$replicates[,7]
       ropt.ha.boot     <- do.boot.ricker$replicates[,8]
       s.boot     <- do.boot.ricker$replicates[,9]

       } # End of the code generating bootstrap data from the observed data set

       yopt.tr.boot <- (ropt.tr.boot-sopt.tr.boot)
       yopt.ha.boot <- (ropt.ha.boot-sopt.ha.boot)
       if(lowerbound) {
          # These new vectors need to be initialized with at least one value in order
          # for S+ to create them in the upcoming 'for' loop
          lb.ha.boot <- 0
          lb.tr.boot <- 0
          for(i in 1:num.boot) {
              if(method==1) lb.ha.boot[i] <- LBesc(a.boot[i],b.boot[i],yopt.ha.boot[i],propMSY,method)$x
              else if(method==2) lb.ha.boot[i] <- LBesc(a.boot[i],b.boot[i],ropt.ha.boot[i],propMSY,method)$x
              if(sopt.tr.boot[i]>2) {
                 if(method==1) lb.tr.boot[i] <- LBesc(a.boot[i],b.boot[i],yopt.tr.boot[i],propMSY,method)$x
                 else if(method==2) lb.tr.boot[i] <- LBesc(a.boot[i],b.boot[i],ropt.tr.boot[i],propMSY,method)$x
                 }
              else
                 lb.tr.boot[i] <- 0
              }
          all.boot.ricker <- cbind(a.boot,b.boot,sopt.tr.boot,uopt.tr.boot,sopt.ha.boot,uopt.ha.boot,ropt.tr.boot,ropt.ha.boot,s.boot,yopt.tr.boot,yopt.ha.boot,lb.tr.boot,lb.ha.boot)
          dimnames(all.boot.ricker) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "Ymsy.tr", "Ymsy.ha", "LB.tr", "LB.ha"))
          }
        else {
          all.boot.ricker <- cbind(a.boot,b.boot,sopt.tr.boot,uopt.tr.boot,sopt.ha.boot,uopt.ha.boot,ropt.tr.boot,ropt.ha.boot,s.boot,yopt.tr.boot,yopt.ha.boot)
          dimnames(all.boot.ricker) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "Ymsy.tr", "Ymsy.ha"))
          }

       if (num.boot!=1000) {
          print(paste("The number of bootstrap replicates in the final sample = ", nrow(all.boot.ricker) ))
			print("Please Note: The final sample includes only replicates where a solution was obtained for Smsy")
          print("*** Summary of the bootstrap simulation results ***")
			cat("\n")
          SumBSdat(all.boot.ricker)

          par(mfrow=c(1,2))
          plot.BSout(all.boot.ricker[,3],sopt.tr1,"MSY Escapement (By Optimization)")
          plot.BSout(all.boot.ricker[,4],uopt.tr1,"MSY Exploitation Rate (By Optimization)")

          plot.BSout(all.boot.ricker[,5],sopt.ha1,"MSY Escapement (Hilborn Approximation)")
          plot.BSout(all.boot.ricker[,6],uopt.ha1,"MSY Exploitation Rate (Hilborn Approximation)")
          }
        else {
          sort.boot.ricker <- sort.col(all.boot.ricker, "<ALL>", 5, T)
          #trunc1.boot.ricker <- sort.boot.ricker[26:975,]
          trunc.boot.ricker <- select.rows(sort.boot.ricker, sort.boot.ricker[,3]>2)
          dimnames(trunc.boot.ricker) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "Ymsy.tr", "Ymsy.ha"))

          print(paste("The number of bootstrap replicates in the final sample = ", nrow(trunc.boot.ricker) ))
			print("The final sample includes only replicates where a solution was obtained for Smsy.")
          print("*** Summary of the bootstrap simulation results ***")
			cat("\n")
          SumBSdat(trunc.boot.ricker)

          par(mfrow=c(1,2))
    	   plot.BSout(trunc.boot.ricker[,3],sopt.tr1,"MSY Escapement (By Optimization)")
          plot.BSout(trunc.boot.ricker[,4],uopt.tr1,"MSY Exploitation Rate (By Optimization)")

          plot.BSout(trunc.boot.ricker[,5],sopt.ha1,"MSY Escapement (Hilborn Approximation)")
          plot.BSout(trunc.boot.ricker[,6],uopt.ha1,"MSY Exploitation Rate (Hilborn Approximation)")
          }

	} # end of bootstrapping
} # end of ricker


if(climatic){

	cat("\n")
	print("#### Ricker Model: log(R/S)=a+bS+SurvRate (SurvRate is untransformed) ###") 


	print("#### First check the relationship between R/S and SurvRate ####")
	par(mfrow=c(1,1))
	plot(sr, log.recruit2spawner, xlab="Survival Rate", ylab="Log (R/S)", pch=16)
   xpips <- par("xaxp")
   ypips <- par("yaxp")
   axis(side=1, ticks=T, labels=F, xaxp=c(xpips[1], xpips[2], xpips[3]*5), tck=-0.0075)
   axis(side=2, ticks=T, labels=F, yaxp=c(ypips[1], ypips[2], ypips[3]*5), tck=-0.0075)

	sr.m  <- lm(log.recruit2spawner~ sr)
	lines(sr, sr.m$fitted.values)  
	
    print("Results from regression analysis of Log recruits per spawner vs survival rate")
	print(summary(sr.m))

	sr  <- data.standard(as.numeric(sr))
	print("Here are the standardized SR values")
	print(sr)

	lm.m  <- lm(log.recruit2spawner~ st+sr) 
	cat("\n")
	print("Below are the results of the regression analysis of the linear form of the Ricker function with SurvRate")
	cat("\n")
 
	print(summary(lm.m))

	par(mfrow=c(2,2))
	plot(lm.m$fitted.values, lm.m$residuals, xlab="Fitted", ylab="Residuals", pch=16)
	abline(h=0)
	plot(year, lm.m$residuals, xlab="Year", ylab="Residuals", pch=16)
	abline(h=0)
	plot(acf(lm.m$residuals, type="correlation", plot=F), main="")
	plot(acf(lm.m$residuals, type="partial", plot=F), main="")

    a <- lm.m$coef[1]
    b <- -lm.m$coef[2]  # b is coerced to a positive here
    g <- lm.m$coef[3]
    s2.sr <- sum(lm.m$residuals**2)/(length(st)-3)
    #s <- sum(lm.m$residuals^2)/lm.m$df.residual # This produces the same result for sigma
    a.cor <- a+0.5*s2.sr
    
    # Smsy by Hilborn approximation method
    sopt.ha2  <- a.cor*(0.5-0.07*a.cor)/b
    uopt.ha2  <- a.cor*(0.5-0.07*a.cor)
    ropt.ha2  <- sopt.ha2*exp(a.cor-b*sopt.ha2)  # Use this when b supplied as positive
    rrep.ha2  <- a.cor/b
#   ropt.ha2  <- sopt.ha2*exp(a.cor+b*sopt.ha2)
    yopt.ha2  <- ropt.ha2-sopt.ha2

    # 'Lower Bound' escapement based on a % of MSY; the neg. sign for 'b'
    # to force it to a positive is necessary
#    if(method==1) LB.s.ha2  <- LBesc(a.cor,b,yopt.ha2,propMSY,method)$x
#    else if(method==2) LB.s.ha2  <- LBesc(a.cor,b,ropt.ha2,propMSY,method)$x
#    LB.r.ha2  <- LB.s.ha2*exp(a.cor-b*LB.s.ha2)  # Use this when b supplied as positive
##   LB.r.ha2  <- LB.s.ha2*exp(a.cor+b*LB.s.ha2)
#    LB.y.ha2  <- LB.r.ha2-LB.s.ha2
	
	# Smsy by solving transcendental equation method
    srparam   <- g*mean(sr)
    sopt.tr2  <- getSoptSR.t(a,b,s2.sr,srparam)$x
    uopt.tr2  <- b*sopt.tr2
    ropt.tr2  <- sopt.tr2*exp(a.cor-b*sopt.tr2)  # Use this when b supplied as positive
    rrep.tr2  <- a.cor/b
#   ropt.tr2  <- sopt.tr2*exp(a.cor+b*sopt.tr2)
    yopt.tr2  <- ropt.tr2-sopt.tr2

    # 'Lower Bound' escapement based on a % of MSY; the neg. sign for 'b'
    # is necessary to force it to a positive 
#    if(method==1) LB.s.tr2  <- LBesc(a.cor,b,yopt.tr2,propMSY,method)$x
#    else if(method==2) LB.s.tr2  <- LBesc(a.cor,b,ropt.tr2,propMSY,method)$x
#    LB.r.tr2  <- LB.s.tr2*exp(a.cor-b*LB.s.tr2)  # Use this when b supplied as positive
##   LB.r.tr2  <- LB.s.tr2*exp(a.cor+b*LB.s.tr2)
#    LB.y.tr2  <- LB.r.tr2-LB.s.tr2

    cat("\n")
    print("Parameter estimates obtained from a Ricker stock-recruit function including survival rate")
	print("Note: The calculation of Smsy includes the bias correction to alpha involving the") 
	print("residual mean square error from the regression results based on the linear form")
	print("of the Ricker function")
	cat("\n")
    print(paste("Alpha=", round(a,4) ))
    print(paste("Beta=", signif(b,4) ))
    print(paste("Gamma=", round(g,4) ))
    print(paste("Sigma^2=", round(s2.sr,4) ))
    print(paste("Gamma x mean SR=", signif(srparam,3) ))
 
    cat("\n")
    print("Estimates of optimal spawners obtained by the Hilborn approximation method")
	cat("\n")
    print(paste("Smsy : ", round(sopt.ha2, 0), "     Umsy : ", round(uopt.ha2,4) ))
    print(paste("Rmsy : ", round(ropt.ha2,0), "    Ymsy : ", round(yopt.ha2,0) ))
 	print(paste("Rrep : ", round(rrep.ha2,0) ))
#    print(paste("Lower Bound Escapement : ", round(LB.s.ha2,0), "    at ", propMSY*100, "% of MSY"))
#    print(paste("Lower Bound Production : ", round(LB.r.ha2,0), "   at ", propMSY*100, "% of MSY"))
#    print(paste("Lower Bound Yield      : ", round(LB.y.ha2,0), "   at ", propMSY*100, "% of MSY"))

    cat("\n")
    print("Estimates of optimal spawners obtained by solving the transcendental equation")
	cat("\n")
    print(paste("Smsy : ", round(sopt.tr2, 0), "     Umsy : ", round(uopt.tr2,4) ))
    print(paste("Rmsy : ", round(ropt.tr2,0), "    Ymsy : ", round(yopt.tr2,0) ))
	print(paste("Rrep : ", round(rrep.tr2,0) ))
#    print(paste("Lower Bound Escapement : ", round(LB.s.tr2,0), "    at ", propMSY*100, "% of MSY"))
#    print(paste("Lower Bound Production : ", round(LB.r.tr2,0), "   at ", propMSY*100, "% of MSY"))
#    print(paste("Lower Bound Yield      : ", round(LB.y.tr2,0), "   at ", propMSY*100, "% of MSY"))

    if(num.boot!=0){
       set.seed(100)
		cat("\n")
       print("##### Doing Bootstrapping for Ricker+SurvRate model ####")

       if(bsresids) { # Do bootstrapping using residuals
          print("Performing bootstrapping using residuals")
          SRdat <- data.frame(y=log.recruit2spawner, s=st, sr=sr, fitted=fitted(lm.m), resids=resid(lm.m))
          do.boot.SR <- boot(SRdat, Bootresids2, R=num.boot)
          print("Summary of bootstrap results:")
			print("Rows indicate the following in sequence:")
			print("1. alpha; 2. beta; 3. Smsy(TR); 4. Umsy(TR); 5. Smsy(HA); 6. Umsy(HA); 7. Rmsy(TR); 8. Rmsy(HA); 9. sigma-squared; 10. gamma")
          print.boot(do.boot.SR)
          cat("\n")
#         plot.boot(do.boot.SR)
          a.boot     <- do.boot.SR$t[,1]
          b.boot     <- do.boot.SR$t[,2]  # b is returned as a positive here
          sopt.tr.boot     <- do.boot.SR$t[,3]
          uopt.tr.boot     <- do.boot.SR$t[,4]
          sopt.ha.boot     <- do.boot.SR$t[,5]
          uopt.ha.boot     <- do.boot.SR$t[,6]
          ropt.tr.boot     <- do.boot.SR$t[,7]
          ropt.ha.boot     <- do.boot.SR$t[,8]
          s.boot     <- do.boot.SR$t[,9]
          g.boot     <- do.boot.SR$t[,10]
          }
       else {  # Do bootstrapping based on original data set
          print("Performing bootstrapping using the observed data")
          do.boot.SR <- bootstrap(dat, boot.SR(dat), B=num.boot)
          print("Summary of the results from bootstrapping:")
          print(do.boot.SR)
          cat("\n")
#         print("Now plotting do.boot.SR")
#         plot(do.boot.SR)
    	
          a.boot     <- do.boot.SR$replicates[,1]
          b.boot     <- do.boot.SR$replicates[,2]   # b is returned as a positive here
          sopt.tr.boot     <- do.boot.SR$replicates[,3]
          uopt.tr.boot     <- do.boot.SR$replicates[,4]
          sopt.ha.boot     <- do.boot.SR$replicates[,5]
          uopt.ha.boot     <- do.boot.SR$replicates[,6]
          ropt.tr.boot     <- do.boot.SR$replicates[,7]
          ropt.ha.boot     <- do.boot.SR$replicates[,8]
          s.boot     <- do.boot.SR$replicates[,9]
          g.boot     <- do.boot.SR$replicates[,10]
          } # End of the code generating bootstrap data from the observed data set

       yopt.tr.boot <- (ropt.tr.boot-sopt.tr.boot)
       yopt.ha.boot <- (ropt.ha.boot-sopt.ha.boot)
       if(lowerbound) {
          # These new vectors need to be initialized with at least one value in order
          # for S+ to create them in the upcoming 'for' loop
          lb.ha.boot <- 0
          lb.tr.boot <- 0
          for(i in 1:num.boot) {
              if(method==1) lb.ha.boot[i] <- LBesc(a.boot[i],b.boot[i],yopt.ha.boot[i],propMSY,method)$x
              else if(method==2) lb.ha.boot[i] <- LBesc(a.boot[i],b.boot[i],ropt.ha.boot[i],propMSY,method)$x
              if(sopt.tr.boot[i]>2) {
                 if(method==1) lb.tr.boot[i] <- LBesc(a.boot[i],b.boot[i],yopt.tr.boot[i],propMSY,method)$x
                 else if(method==2) lb.tr.boot[i] <- LBesc(a.boot[i],b.boot[i],ropt.tr.boot[i],propMSY,method)$x
                 }
              else
                 lb.tr.boot[i] <- 0
              }
          all.boot.SR <- cbind(a.boot,b.boot,sopt.tr.boot,uopt.tr.boot,sopt.ha.boot,uopt.ha.boot,ropt.tr.boot,ropt.ha.boot,s.boot,g.boot,yopt.tr.boot,yopt.ha.boot,lb.tr.boot,lb.ha.boot)
          dimnames(all.boot.SR) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "gamma","Ymsy.tr", "Ymsy.ha", "LB.tr", "LB.ha"))
          }
        else {
          all.boot.SR <- cbind(a.boot,b.boot,sopt.tr.boot,uopt.tr.boot,sopt.ha.boot,uopt.ha.boot,ropt.tr.boot,ropt.ha.boot,s.boot,g.boot,yopt.tr.boot,yopt.ha.boot)
          dimnames(all.boot.SR) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "gamma", "Ymsy.tr", "Ymsy.ha"))
          }

       if (num.boot!=1000) {
          print(paste("The number of bootstrap replicates in the final sample = ", nrow(all.boot.SR) ))
          print("*** Summary of the bootstrap simulation results ***")
			cat("\n")
          SumBSdat(all.boot.SR)

          par(mfrow=c(1,2))
          plot.BSout(all.boot.SR[,3],sopt.tr2,"MSY Escapement (By Optimization)")
          plot.BSout(all.boot.SR[,4],uopt.tr2,"MSY Exploitation Rate (By Optimization)")

          plot.BSout(all.boot.SR[,5],sopt.ha2,"MSY Escapement (Hilborn Approximation)")
          plot.BSout(all.boot.SR[,6],uopt.ha2,"MSY Exploitation Rate (Hilborn Approximation)")
          }
       else {
          sort.boot.SR <- sort.col(all.boot.SR, "<ALL>", 5, T)
          #trunc1.boot.SR <- sort.boot.SR[26:975,]
          trunc.boot.SR <- select.rows(sort.boot.SR, sort.boot.SR[,3]>2)
          dimnames(trunc.boot.SR) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "gamma", "Ymsy.tr", "Ymsy.ha"))

          print(paste("The number of bootstrap replicates in the final sample = ", nrow(trunc.boot.SR) ))
			print("Please Note: The final sample includes only replicates where a solution was obtained for Smsy")
          print("*** Summary of the bootstrap simulation results ***")
			cat("\n")
          SumBSdat(trunc.boot.SR)

          par(mfrow=c(1,2))
          plot.BSout(trunc.boot.SR[,3],sopt.tr2,"MSY Escapement (By Optimization)")
          plot.BSout(trunc.boot.SR[,4],uopt.tr2,"MSY Exploitation Rate (By Optimization)")

          plot.BSout(trunc.boot.SR[,5],sopt.ha2,"MSY Escapement (Hilborn Approximation)")
          plot.BSout(trunc.boot.SR[,6],uopt.ha2,"MSY Exploitation Rate (Hilborn Approximation)")
          }

	} # end of bootstrap

    print("At end of climatic")
} # end of climatic

if(logclimatic){

	cat("\n")
	print("#### Ricker Model: log(R/S)=a+bS+ln(SurvRate) (SurvRate is transformed) ###") 


	print("#### First check the relationship between R/S and ln(SurvRate) ####")
	par(mfrow=c(1,1))
	logsr  <- log(dat$SR) 
#	print(cbind(sr,logsr))
	plot(logsr, log.recruit2spawner, xlab="Ln Survival Rate", ylab="Ln (R/S)", pch=16)
   xpips <- par("xaxp")
   ypips <- par("yaxp")
   axis(side=1, labels=F, xaxp=c(xpips[1], xpips[2], xpips[3]*5), tck=-0.0075)
   axis(side=2, labels=F, yaxp=c(ypips[1], ypips[2], ypips[3]*5), tck=-0.0075)
	logsr.m  <- lm(log.recruit2spawner~ logsr)
	lines(logsr, logsr.m$fitted.values)  
	cat("\n")
    print("Results from regression analysis of Ln recruits per spawner vs Ln survival rate")
	cat("\n")
    print(summary(logsr.m))

	logsr  <- data.standard(log(dat$SR))
	print("Here are the standardized ln SR values")
	print(logsr)
	
	loglm.m  <- lm(log.recruit2spawner~ st+logsr)

	cat("\n")
	print("Below are the results of the regression analysis of the linear form of the Ricker function with log SurvRate")
	cat("\n")
	print(summary(loglm.m))
	print("Residuals from the linear regression")
	print(loglm.m$residuals)

	par(mfrow=c(2,2))
	plot(loglm.m$fitted.values, loglm.m$residuals, xlab="Fitted", ylab="Residuals", pch=16)
	abline(h=0)
	plot(year, loglm.m$residuals, xlab="Year", ylab="Residuals", pch=16)
	abline(h=0)
	plot(acf(loglm.m$residuals, type="correlation", plot=F), main="")
	plot(acf(loglm.m$residuals, type="partial", plot=F), main="")

    a <- loglm.m$coef[1]
    b <- -loglm.m$coef[2]  # b is coerced to a positive here
    g <- loglm.m$coef[3]
    s2.logsr <- sum(loglm.m$residuals**2)/(length(st)-3)
	#s <- sum(loglm.m$residuals^2)/loglm.m$df.residual # This produces the same result for sigma
    a.cor <- a+0.5*s2.logsr
    
    # Smsy by Hilborn approximation method
    sopt.ha3  <- a.cor*(0.5-0.07*a.cor)/b
    uopt.ha3  <- a.cor*(0.5-0.07*a.cor)
    ropt.ha3  <- sopt.ha3*exp(a.cor-b*sopt.ha3)  # Use when b supplied as positive
    rrep.ha3  <- a.cor/b
#   ropt.ha3  <- sopt.ha3*exp(a.cor+b*sopt.ha3)
    yopt.ha3  <- ropt.ha3-sopt.ha3

    # 'Lower Bound' escapement based on a proportion of MSY
    # The proportion is set as a parameter on the command line
#    if(method==1) LB.s.ha3  <- LBesc(a.cor,b,yopt.ha3,propMSY,method)$x
#    else if(method==2) LB.s.ha3  <- LBesc(a.cor,b,ropt.ha3,propMSY,method)$x
#    LB.r.ha3  <- LB.s.ha3*exp(a.cor-b*LB.s.ha3)  # Use when b supplied as positive
##   LB.r.ha3  <- LB.s.ha3*exp(a.cor+b*LB.s.ha3)
#    LB.y.ha3  <- LB.r.ha3-LB.s.ha3
	
	# Smsy by solving transcendental equation method
    srparam   <- g*mean(logsr)
    sopt.tr3  <- getSoptSR.t(a,b,s2.logsr,srparam)$x
    uopt.tr3  <- b*sopt.tr3
    ropt.tr3  <- sopt.tr3*exp(a.cor-b*sopt.tr3)  # Use when b supplied as positive
    rrep.tr3  <- a.cor/b
#   ropt.tr3  <- sopt.tr3*exp(a.cor+b*sopt.tr3)
    yopt.tr3  <- ropt.tr3-sopt.tr3

    # 'Lower Bound' escapement based on a proportion of MSY
    # The proportion is supplied as a parameter on the command line
#    if(method==1) LB.s.tr3  <- LBesc(a.cor,b,yopt.tr3,propMSY,method)$x
#    else if(method==2) LB.s.tr3  <- LBesc(a.cor,b,ropt.tr3,propMSY,method)$x
#    LB.r.tr3  <- LB.s.tr3*exp(a.cor-b*LB.s.tr3)  # Use when b supplied as positive
##   LB.r.tr3  <- LB.s.tr3*exp(a.cor+b*LB.s.tr3)
#    LB.y.tr3  <- LB.r.tr3-LB.s.tr3

    cat("\n")
    print("Parameter estimates obtained from a Ricker stock-recruit function with log (SR)")
	print("Note: The calculation of Smsy includes the bias correction to alpha involving the") 
	print("residual mean square error from the regression results based on the linear form")
	print("of the Ricker function")
	cat("\n")
    print(paste("Alpha=", round(a,4) ))
    print(paste("Beta=", signif(b,4) ))
    print(paste("Gamma=", round(g,4) ))
    print(paste("Sigma^2=", round(s2.logsr,4) ))
    print(paste("Gamma x mean logSR=", signif(srparam,3) ))
 
    cat("\n")
    print("Estimates of optimal spawners obtained by the Hilborn approximation method")
	cat("\n")
    print(paste("Smsy : ", round(sopt.ha3, 0), "     Umsy : ", round(uopt.ha3,4) ))
    print(paste("Rmsy : ", round(ropt.ha3,0), "    Ymsy : ", round(yopt.ha3,0) ))
	print(paste("Rrep : ", round(rrep.ha3,0) ))
#    print(paste("Lower Bound Escapement : ", round(LB.s.ha3,0), "    at ", propMSY*100, "% of MSY"))
#    print(paste("Lower Bound Production : ", round(LB.r.ha3,0), "   at ", propMSY*100, "% of MSY"))
#    print(paste("Lower Bound Yield      : ", round(LB.y.ha3,0), "   at ", propMSY*100, "% of MSY"))

    cat("\n")
    print("Estimates of optimal spawners obtained by solving the transcendental equation")
	cat("\n")
    print(paste("Smsy : ", round(sopt.tr3, 0), "     Umsy : ", round(uopt.tr3,4) ))
    print(paste("Rmsy : ", round(ropt.tr3,0), "    Ymsy : ", round(yopt.tr3,0) ))
	print(paste("Rrep : ", round(rrep.tr3,0) ))
#    print(paste("Lower Bound Escapement : ", round(LB.s.tr3,0), "    at ", propMSY*100, "% of MSY"))
#    print(paste("Lower Bound Production : ", round(LB.r.tr3,0), "   at ", propMSY*100, "% of MSY"))
#    print(paste("Lower Bound Yield      : ", round(LB.y.tr3,0), "   at ", propMSY*100, "% of MSY"))

    if(num.boot!=0){
       set.seed(100)
		cat("\n")
       print("##### Doing Bootstrapping for Ricker+log(SurvRate) model ###")

       if(bsresids) { # Do bootstrapping using residuals
          print("Performing bootstrapping using residuals")
    	   SRdat <- data.frame(y=log.recruit2spawner, s=st, sr=logsr, fitted=fitted(loglm.m), resids=resid(loglm.m))
          do.boot.logSR <- boot(SRdat, Bootresids2, R=num.boot)
          print("Summary of bootstrap results:")
			print("Rows indicate the following in sequence:")
			print("1. alpha; 2. beta; 3. Smsy(TR); 4. Umsy(TR); 5. Smsy(HA); 6. Umsy(HA); 7. Rmsy(TR); 8. Rmsy(HA); 9. sigma-squared; 10. gamma")
			cat("\n")
          print.boot(do.boot.logSR)
          cat("\n")
#         plot.boot(do.boot.logSR)
          a.boot     <- do.boot.logSR$t[,1]
          b.boot     <- do.boot.logSR$t[,2]  # b is returned as a positive here
          sopt.tr.boot     <- do.boot.logSR$t[,3]
          uopt.tr.boot     <- do.boot.logSR$t[,4]
          sopt.ha.boot     <- do.boot.logSR$t[,5]
          uopt.ha.boot     <- do.boot.logSR$t[,6]
          ropt.tr.boot     <- do.boot.logSR$t[,7]
          ropt.ha.boot     <- do.boot.logSR$t[,8]
          s.boot     <- do.boot.logSR$t[,9]
          g.boot     <- do.boot.logSR$t[,10]
          }
	    else {  # Do bootstrapping based on original data set
          print("Performing bootstrapping using the observed data")
          do.boot.logSR <- bootstrap(dat, boot.logSR(dat), B=num.boot)
          print("Summary of the results from bootstrapping:")
			cat("\n")
          print(do.boot.logSR)
          cat("\n")
 #        print("Now plotting do.boot.logSR")
 #        plot(do.boot.logSR)

          a.boot     <- do.boot.logSR$replicates[,1]
          b.boot     <- do.boot.logSR$replicates[,2] # b is returned as a positive here
          sopt.tr.boot     <- do.boot.logSR$replicates[,3]
          uopt.tr.boot     <- do.boot.logSR$replicates[,4]
          sopt.ha.boot     <- do.boot.logSR$replicates[,5]
          uopt.ha.boot     <- do.boot.logSR$replicates[,6]
          ropt.tr.boot     <- do.boot.logSR$replicates[,7]
          ropt.ha.boot     <- do.boot.logSR$replicates[,8]
          s.boot     <- do.boot.logSR$replicates[,9]
          g.boot     <- do.boot.logSR$replicates[,10]
          } # End of the code generating bootstrap data from the observed data set

       yopt.tr.boot <- (ropt.tr.boot-sopt.tr.boot)
       yopt.ha.boot <- (ropt.ha.boot-sopt.ha.boot)
       if(lowerbound) {
          # These new vectors need to be initialized with at least one value in order
          # for S+ to create them in the upcoming 'for' loop
          lb.ha.boot <- 0
          lb.tr.boot <- 0
          for(i in 1:num.boot) {
              if(method==1) lb.ha.boot[i] <- LBesc(a.boot[i],b.boot[i],yopt.ha.boot[i],propMSY,method)$x
              else if(method==2) lb.ha.boot[i] <- LBesc(a.boot[i],b.boot[i],ropt.ha.boot[i],propMSY,method)$x
              if(sopt.tr.boot[i]>2) {
                 if(method==1) lb.tr.boot[i] <- LBesc(a.boot[i],b.boot[i],yopt.tr.boot[i],propMSY,method)$x
                 else if(method==2) lb.tr.boot[i] <- LBesc(a.boot[i],b.boot[i],ropt.tr.boot[i],propMSY,method)$x
                 }
              else
                 lb.tr.boot[i] <- 0
              }
          all.boot.logSR <- cbind(a.boot,b.boot,sopt.tr.boot,uopt.tr.boot,sopt.ha.boot,uopt.ha.boot,ropt.tr.boot,ropt.ha.boot,s.boot,g.boot,yopt.tr.boot,yopt.ha.boot,lb.tr.boot,lb.ha.boot)
          dimnames(all.boot.logSR) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "gamma","Ymsy.tr", "Ymsy.ha", "LB.tr", "LB.ha"))
          }
        else {
          all.boot.logSR <- cbind(a.boot,b.boot,sopt.tr.boot,uopt.tr.boot,sopt.ha.boot,uopt.ha.boot,ropt.tr.boot,ropt.ha.boot,s.boot,g.boot,yopt.tr.boot,yopt.ha.boot)
          dimnames(all.boot.logSR) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "gamma", "Ymsy.tr", "Ymsy.ha"))
          }

      if(num.boot!=1000) {  
         print(paste("The number of bootstrap replicates in the final sample = ", nrow(all.boot.logSR) ))

         print("*** Summary of the bootstrap simulation results ***")
		  cat("\n")
         SumBSdat(all.boot.logSR)

         par(mfrow=c(1,2))
         plot.BSout(all.boot.logSR[,3],sopt.tr3,"MSY Escapement (By Optimization)")
         plot.BSout(all.boot.logSR[,4],uopt.tr3,"MSY Exploitation Rate (By Optimization)")

         plot.BSout(all.boot.logSR[,5],sopt.ha3,"MSY Escapement (Hilborn Approximation)")
         plot.BSout(all.boot.logSR[,6],uopt.ha3,"MSY Exploitation Rate (Hilborn Approximation)")
         }
      else {
         sort.boot.logSR <- sort.col(all.boot.logSR, "<ALL>", 5, T)
         #trunc1.boot.logSR <- sort.boot.logSR[26:975,]
         trunc.boot.logSR <- select.rows(sort.boot.logSR, sort.boot.logSR[,3]>2)
         dimnames(trunc.boot.logSR) <- list(NULL, c("alpha", "beta", "Smsy.tr", "Umsy.tr", "Smsy.ha", "Umsy.ha", "Rmsy.tr", "Rmsy.ha", "sigma", "gamma", "Ymsy.tr", "Ymsy.ha"))

         print(paste("The number of bootstrap replicates in the final sample = ", nrow(trunc.boot.logSR) ))
 		  print("Please Note: The final sample includes only replicates where a solution was obtained for Smsy")
         print("*** Summary of the bootstrap simulation results ***")
		  cat("\n")
         SumBSdat(trunc.boot.logSR)

         par(mfrow=c(1,2))
         plot.BSout(trunc.boot.logSR[,3],sopt.tr3,"MSY Escapement")
         plot.BSout(trunc.boot.logSR[,4],uopt.tr3,"MSY Exploitation Rate")
#        plot.BSout(trunc.boot.logSR[,3],sopt.tr3,"MSY Escapement (Transcendental Equation)")
#        plot.BSout(trunc.boot.logSR[,4],uopt.tr3,"MSY Exploitation Rate (Transcendental Equation)")

         plot.BSout(trunc.boot.logSR[,5],sopt.ha3,"MSY Escapement (Hilborn Approximation)")
         plot.BSout(trunc.boot.logSR[,6],uopt.ha3,"MSY Exploitation Rate (Hilborn Approximation)")
         }

	} # end of bootstrap
   print("At end of logclimatic")

} # end of logclimatic
	
    
if(do.plot)	{
	par(mfrow=c(1,1))
	x.max1 <- c(0,0,0)
	if(ricker) {
		smax <- -(ricker.m$coef[1]+s2.r/2)/ricker.m$coef[2]
		x.max1[1] <- smax+smax*0.04
#    	print(paste("Value of x.max1[1] :  ", x.max1[1]))
	}
   if(climatic) {
		smax <- -(lm.m$coef[1]+s2.sr/2)/lm.m$coef[2]
		x.max1[2] <- smax+smax*0.04
#		print(paste("Value of x.max1[2] :  ", x.max1[2]))
	}
	if(logclimatic) {
		smax <- -(loglm.m$coef[1]+s2.logsr/2)/loglm.m$coef[2]
		x.max1[3] <- smax+smax*0.04
#    	print(paste("Value of x.max1[3] :  ", x.max1[3]))
	}
   x.max <- max(x.max1)
#	print(paste("x.max = ", x.max))
 	plot(st,rt, xlim=c(0, x.max), ylim=c(0, max(rt)), xlab="Spawners", ylab="Recruits", pch=16)
   xpips <- par("xaxp")
   ypips <- par("yaxp")
   axis(side=1, ticks=T, labels=F, xaxp=c(xpips[1], xpips[2], xpips[3]*5), tck=-0.0075)
   axis(side=2, ticks=T, labels=F, yaxp=c(ypips[1], ypips[2], ypips[3]*5), tck=-0.0075)
   #print(paste("newmin = ", xpips[1], "   newmax = ", xpips[2], "  n ints = ", xpips[3]))
   #Below is a slightly different way to obtain pip marks on a graph
   #pips <- seq(from=new[1], to=new[2], length=new[3]*2+1)
   #print(pips)
   #axis(side=1, at=pips, ticks=T, labels=F, tck=-0.0075)
	lines(c(0,x.max), c(0,x.max	))
	
	x  <- seq(0, x.max, length=150)

if(ricker){
    # The equation below is correct because b is a negative
 	r1  <- exp(ricker.m$coef[1]+(s2.r/2))*x*exp(ricker.m$coef[2]*x)
# 	r1  <- exp(ricker.m$coef[1])*x*exp(ricker.m$coef[2]*x)
	lines(x, r1, lty=1)
    segments(x1=sopt.tr1, y1=0, x2=sopt.tr1, y2=ropt.tr1, lty=1)
#	abline(v=sopt.ha1, lty=1)  # Smsy - Hilborn Approx.
#	abline(v=sopt.tr1, lty=1)  # Smsy - transcendental
    print("Ricker curve plotted")
	} # end if ricker

if(climatic){
   # The equation below is correct because b is a negative
 	r2  <- exp(lm.m$coef[1]+(s2.sr/2))*x*exp(lm.m$coef[2]*x)
#	r2  <- exp(lm.m$coef[1])*x*exp(lm.m$coef[2]*x)
    if (!ricker && !logclimatic) {
       lines(x, r2, lty=1)
       segments(x1=sopt.tr2, y1=0, x2=sopt.tr2, y2=ropt.tr2, lty=1)
#      abline(v=sopt.ha2,lty=1)   # Smsy - Hilborn Approx.
#      abline(v=sopt.tr2,lty=1)   # Smsy - transcendental
       }
    else {
       lines(x, r2, lty=8)
       segments(x1=sopt.tr2, y1=0, x2=sopt.tr2, y2=ropt.tr2, lty=8)
#      abline(v=sopt.ha2,lty=8)   # Smsy - Hilborn Approx.
#      abline(v=sopt.tr2,lty=8)   # Smsy - transcendental
       }
    print("Climatic curve plotted")
	} # end of if climatic

if(logclimatic){
    # The equation below is correct because b is a negative
    # The calculation of production below includes the correction for
    # process and sampling error
    r3  <- exp(loglm.m$coef[1]+(s2.logsr/2))*x*exp(loglm.m$coef[2]*x)
    # The calculation of production below does not include a correction for
    # process and sampling error
#   r3  <- exp(loglm.m$coef[1])*x*exp(loglm.m$coef[2]*x)
    if (!ricker && !climatic) {
       lines(x, r3, lty=1)
       segments(x1=sopt.tr3, y1=0, x2=sopt.tr3, y2=ropt.tr3, lty=1)
#      segments(x1=0, y1=ropt.tr3, x2=sopt.tr3, y2=ropt.tr3)
#      abline(v=sopt.ha3,lty=1)  # Smsy - Hilborn Approx.
#      abline(v=sopt.tr3,lty=1)  # Smsy - transcendental
       }
    else {
       lines(x, r3, lty=4)
       segments(x1=sopt.tr3, y1=0, x2=sopt.tr3, y2=ropt.tr3, lty=4)
#      abline(v=sopt.ha3,lty=4)  # Smsy - Hilborn Approx.
#      abline(v=sopt.tr3,lty=4)  # Smsy - transcendental
       }
    print("LogClimatic curve plotted")
	} # end of if log climatic

if(climatic){
	par(mfrow=c(1,1))
 	plot(st,rt, xlim=c(0, x.max), ylim=c(0, max(rt)), xlab="Spawners", ylab="Recruits", pch=16)
   axis(side=1, ticks=T, labels=F, xaxp=c(xpips[1], xpips[2], xpips[3]*5), tck=-0.0075)
   axis(side=2, ticks=T, labels=F, yaxp=c(ypips[1], ypips[2], ypips[3]*5), tck=-0.0075)
	lines(c(0,x.max), c(0,x.max))

    # The equation below is correct because b is a negative
	r2.adj  <- rt*exp(-lm.m$coef[3]*sr)
	
	print("Here are the R values adjusted for the co-variate")
	print(r2.adj)

	lines(x, r2, lty=1)
    segments(x1=sopt.tr2, y1=0, x2=sopt.tr2, y2=ropt.tr2, lty=1)
#	abline(v=sopt.ha2,lty=1)  # Smsy - Hilborn Approx.
#	abline(v=sopt.tr2,lty=1)  # Smsy - transcendental
	points(st, r2.adj, pch=0)
    # The following three lines are what produce the vertical line segments
    # extending between the real 'R's' and 'R's' adjusted for survival rate
    for (i in 1:len.dat){
        segments(st[i], rt[i], st[i], r2.adj[i])
        }

	# Graph Number 7 - Shows S-R curve with adjusted R values as points
	par(mfrow=c(1,1))
 	plot(st,r2.adj, xlim=c(0, x.max), ylim=c(0, max(rt)), xlab="Spawners", ylab="Recruits", pch=0)
   axis(side=1, ticks=T, labels=F, xaxp=c(xpips[1], xpips[2], xpips[3]*5), tck=-0.0075)
   axis(side=2, ticks=T, labels=F, yaxp=c(ypips[1], ypips[2], ypips[3]*5), tck=-0.0075)
	lines(c(0,x.max), c(0,x.max))
	lines(x, r2, lty=1)
    segments(x1=sopt.tr2, y1=0, x2=sopt.tr2, y2=ropt.tr2, lty=1)
#	abline(v=sopt.ha2,lty=1)  # Smsy - Hilborn Approx.
#	abline(v=sopt.tr2,lty=1)  # Smsy - transcendental
#	points(st, r2.adj, pch=0)

	} # end of if climatic

if(logclimatic){
	par(mfrow=c(1,1))
 	plot(st,rt, xlim=c(0, x.max), ylim=c(0, max(rt)), xlab="Spawners", ylab="Recruits", pch=16)
   axis(side=1, ticks=T, labels=F, xaxp=c(xpips[1], xpips[2], xpips[3]*5), tck=-0.0075)
   axis(side=2, ticks=T, labels=F, yaxp=c(ypips[1], ypips[2], ypips[3]*5), tck=-0.0075)
	lines(c(0,x.max), c(0,x.max))

    # The equation below is correct because b is a negative
	r3.adj  <- rt*exp(-loglm.m$coef[3]*logsr)  #logsr = index values

	print("Here are the R values adjusted for the co-variate")
	print(r3.adj)

	lines(x, r3, lty=1)
   segments(x1=sopt.tr3, y1=0, x2=sopt.tr3, y2=ropt.tr3, lty=1)
#	abline(v=sopt.ha3,lty=1)  # Smsy - Hilborn Approx.
# 	abline(v=sopt.tr3,lty=1)  # Smsy - transcendental
	points(st, r3.adj, pch=0)
    # The following three lines are what produce the vertical line segments
    # extending between the real 'R's' and 'R's' adjusted for survival rate
    for (i in 1:len.dat){
       segments(st[i], rt[i], st[i], r3.adj[i])
       }

	# Graph Number 7 - Shows S-R curve with adjusted R values as points
	par(mfrow=c(1,1))
 	plot(st,r3.adj, xlim=c(0, x.max), ylim=c(0, max(rt)), xlab="Spawners", ylab="Recruits", pch=0)
   axis(side=1, ticks=T, labels=F, xaxp=c(xpips[1], xpips[2], xpips[3]*5), tck=-0.0075)
   axis(side=2, ticks=T, labels=F, yaxp=c(ypips[1], ypips[2], ypips[3]*5), tck=-0.0075)
	lines(c(0,x.max), c(0,x.max))
	lines(x, r3, lty=1)
   segments(x1=sopt.tr3, y1=0, x2=sopt.tr3, y2=ropt.tr3, lty=1)
#	abline(v=sopt.ha3,lty=1)  # Smsy - Hilborn Approx.
# 	abline(v=sopt.tr3,lty=1)  # Smsy - transcendental
#	points(st, r3.adj, pch=0)

	} # end of if logclimatic

} # end of if plot

if(num.boot!=0) # save the boostrapping results
{
	if(num.boot!=1000) {
      if(ricker)
         outdat <- list(boot.R=do.boot.ricker, all.boot.Ricker=all.boot.ricker, SRdat=dat)
      if(climatic)
         outdat <- list(boot.SR=do.boot.SR, all.boot.SR=all.boot.SR, SRdat=dat)
      if(logclimatic)
         outdat <- list(boot.logSR=do.boot.logSR, all.boot.logSR=all.boot.logSR, SRdat=dat)
      }
    else {
      if(ricker)
         outdat <- list(boot.R=do.boot.ricker, all.BSreps.Ricker=sort.boot.ricker, trunc.BSreps.Ricker=trunc.boot.ricker, SRdat=dat)
      if(climatic)
         outdat <- list(boot.SR=do.boot.SR, all.BSreps.SR=sort.boot.SR, trunc.boot.SR=trunc.BSreps.SR, SRdat=dat, Radj=r2.adj)
      if(logclimatic)
         outdat <- list(boot.logSR=do.boot.logSR, all.BSreps.logSR=sort.boot.logSR, trunc.BSreps.logSR=trunc.boot.logSR, SRdat=dat, Radj=r3.adj)
      }

# 	if(ricker & climatic & logclimatic)
# 	list(boot.ricker=do.boot.ricker, boot.SR=do.boot.SR, boot.logSR=do.boot.logSR) 
# 	if(ricker & logclimatic)
# 	list(boot.ricker=do.boot.ricker, boot.logSR=do.boot.logSR) 

    # This return statement seems to be necessary in order for data to
    # saved at the end of processing.
    return(outdat)

    } # end of saving output results from bootstrap analysis

}  # end of SR.fit function

######## end all Subs  ############

######## now it is time to run the data #########

# command line to process data
#library(boot)



