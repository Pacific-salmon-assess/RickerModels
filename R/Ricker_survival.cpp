#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

 // dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs_R);   // observed log recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  DATA_VECTOR(obs_survival);    // observed  survival up to age 2
  
  PARAMETER(alpha);
  //PARAMETER(atwosurv);
  PARAMETER(logbeta);
  //PARAMETER(rho);
  PARAMETER(logSigObs);
 // PARAMETER(logSigtheta);
  //PARAMETER_VECTOR(logsurv)
  
  int timeSteps=obs_R.size();

  Type beta=exp(logbeta);
  Type SigObs     = exp(logSigObs);
  //Type Sigtheta     = exp(logSigtheta);
  Type Smax  = Type(1.0)/beta;

  
  Type tau     = Type(1.0)/(SigObs*SigObs);
  
  vector<Type> pred_logRS(timeSteps), obs_logRS(timeSteps), residuals(timeSteps), alphat(timeSteps),Smsy(timeSteps), umsy(timeSteps);

 

  //priors on precision and variance ratio
  //Type ans= -dbeta(rho,Type(3.0),Type(3.0),true);  
  //Type ans= -dnorm(logSigObs,Type(0.0),Type(5.0),true);   
    Type ans= Type(0);
 
  

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_survival(i))){
      obs_logRS(i) = log(obs_R(i)/obs_S(i));
      pred_logRS(i) = alpha + log(obs_survival(i)) - beta * obs_S(i) ;
      alphat(i) = alpha + log(obs_survival(i));
      //pred_logR(i) = logRS(i) + log(obs_S(i));
      residuals(i) = obs_logRS(i) - pred_logRS(i);
      umsy(i)     = Type(.5) * alphat(i) - Type(0.07) * (alphat(i) * alphat(i));
      Smsy(i)     =  alphat(i)/beta * (Type(0.5) -Type(0.07) * alphat(i));  
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),SigObs,true);
    }
  
  }

 


  REPORT(pred_logRS)
  REPORT(alpha)
  REPORT(alphat)
  REPORT(residuals)
  //REPORT(atwosurv)
  REPORT(tau)
  REPORT(beta)
  REPORT(SigObs)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  return ans;
}

