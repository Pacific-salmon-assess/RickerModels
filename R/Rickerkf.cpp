#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs_logR);   // observed log recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  PARAMETER(alphao);
  PARAMETER(Smax);
 
  PARAMETER(logSigalpha);
  PARAMETER(logSigObs);

  PARAMETER_VECTOR(alpha);
  
  int timeSteps=obs_logR.size();
  
  Type Sigalpha=exp(logSigalpha);
  Type SigObs=exp(logSigObs);
  Type beta = Type(1.0)/Smax;


  vector<Type> pred_logR(timeSteps), logRS(timeSteps);


  Type ans= -dnorm(alpha(0),alphao,Sigalpha,true); 
  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm(alpha(i),alpha(i-1),Sigalpha,true);
  
  }

  for(int i=0;i<timeSteps;i++){
    
    logRS(i) = alpha(i) - beta * obs_S(i) ;
    pred_logR(i) = logRS(i) + log(obs_S(i)); 
    ans+=-dnorm(obs_logR(i),pred_logR(i),SigObs,true);
  
  }

  REPORT(pred_logR)
  REPORT( alpha)
  REPORT(SigObs)
  REPORT(Sigalpha)
  REPORT(beta )
  REPORT(alphao)
  return ans;
}

