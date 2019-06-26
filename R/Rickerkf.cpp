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
  
  PARAMETER(alphao);
  PARAMETER(logbeta);
 
  PARAMETER(logSigalpha);
  PARAMETER(logSigObs);

  PARAMETER_VECTOR(alpha);
  
  int timeSteps=obs_logR.size();

  Type beta=exp(logbeta);
  
  Type Sigalpha=exp(logSigalpha);
  Type SigObs=exp(logSigObs);
  Type Smax  = Type(1.0)/beta;

  //Type tauR = Type(1.0)/SigObs;
  //Type taualpha = Type(1.0)/Sigalpha;


  vector<Type> pred_logR(timeSteps), logRS(timeSteps);

  //Type ans= -dgamma(tauR,Type(0.01),Type(0.001),true);   
  //ans+= -dgamma(taualpha,Type(0.01),Type(0.001),true); 

  //Type ans= -dgamma(tauR,Type(0.01),Type(0.00001),true);  
  //ans+= -dgamma(taualpha,Type(0.01),Type(0.00001),true); 

  Type ans= -dnorm(logSigalpha,Type(0.0),Type(1.0),true);  
  ans+= -dnorm(logSigObs,Type(0.0),Type(1.0),true); 

  ans+= -dnorm(alpha(0),alphao,Sigalpha,true); 
  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm(alpha(i),alpha(i-1),Sigalpha,true);
  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logR(i))){
      logRS(i) = alpha(i) - beta * obs_S(i) ;
      pred_logR(i) = logRS(i) + log(obs_S(i)); 
      ans+=-dnorm(obs_logR(i),pred_logR(i),SigObs,true);
    }
  
  }

  REPORT(pred_logR)
  REPORT( alpha)
  REPORT(SigObs)
  REPORT(Sigalpha)
  REPORT(beta)
  REPORT(alphao)
  REPORT(Smax)
  return ans;
}

