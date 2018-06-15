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
  DATA_VECTOR(obs_logR);   // observed log recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  PARAMETER(alpha);
  PARAMETER(logbeta);
  //PARAMETER(rho);
  PARAMETER(logSigObs);
  
  int timeSteps=obs_logR.size();

  Type beta=exp(logbeta);
  Type SigObs     = exp(logSigObs);
  Type Smax  = Type(1.0)/beta;
  
  Type tau     = Type(1.0)/(SigObs*SigObs);
  
  vector<Type> pred_logR(timeSteps), logRS(timeSteps);

 

  //priors on precision and variance ratio
  //Type ans= -dbeta(rho,Type(3.0),Type(3.0),true);  
  //Type ans= -dnorm(logSigObs,Type(0.0),Type(5.0),true);   
    Type ans= Type(0);
  

  
  

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logR(i))){
      logRS(i) = alpha - beta * obs_S(i) ;
      pred_logR(i) = logRS(i) + log(obs_S(i)); 
      ans+=-dnorm(obs_logR(i),pred_logR(i),SigObs,true);
    }
  
  }

  REPORT(pred_logR)
  REPORT(alpha)
  REPORT(tau)
  REPORT(beta)
  REPORT(SigObs)
  REPORT(Smax)
  return ans;
}

