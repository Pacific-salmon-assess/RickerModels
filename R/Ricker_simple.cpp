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
  
  PARAMETER(alpha);
  PARAMETER(logbeta);
  //PARAMETER(rho);
  PARAMETER(logSigObs);
  
  int timeSteps=obs_R.size();

  Type beta = exp(logbeta);
  Type SigObs = exp(logSigObs);
  Type Smax  = Type(1.0)/beta;
  
  Type tau     = Type(1.0)/(SigObs*SigObs);
  
  vector<Type> pred_logRS(timeSteps), obs_logRS(timeSteps), residuals(timeSteps); 
  
  Type ans= Type(0);
  

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_R(i))){
      obs_logRS(i) = log(obs_R(i)/obs_S(i)) ;
      pred_logRS(i) = alpha - beta * obs_S(i) ; 
      residuals(i) = obs_logRS(i) - pred_logRS(i);
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),SigObs,true);
    }
  
  }
  Type umsy     = Type(.5) * alpha - Type(0.07) * (alpha * alpha);

  REPORT(pred_logRS)
  REPORT(residuals)
  REPORT(alpha)
  REPORT(tau)
  REPORT(beta)
  REPORT(SigObs)
  REPORT(Smax)
  REPORT(umsy)
  return ans;
}

