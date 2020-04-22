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
  
  DATA_SCALAR(prbeta1);
  DATA_SCALAR(prbeta2);
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value
  //rho         -> Proportion of total variance associated with obs error.
  //varphi      -> Total precision
  //alpha       -> Time-varying alpha

  PARAMETER(alphao);
  PARAMETER(logbeta);
  PARAMETER(rho);
  PARAMETER(logvarphi);

  PARAMETER_VECTOR(alpha);
  

  
  int timeSteps=obs_logR.size();

  Type beta = exp(logbeta);
  Type Smax  = Type(1.0)/beta;
  
  //theta       -> total standard deviation
  //sig         -> obs error std
  //tau         -> proc error (alpha) std
  
  Type varphi     = exp(logvarphi);
  Type theta     = sqrt(Type(1.0)/varphi);
  Type sig       = sqrt(rho) * theta;
  Type tau        = sqrt(Type(1.0)-rho) * theta ;


  vector<Type> pred_logR(timeSteps), logRS(timeSteps),umsy(timeSteps), Smsy(timeSteps);

  

  //priors on precision and variance ratio
  Type ans= -dbeta(rho,prbeta1,prbeta2,true);  
  
  ans+= -dnorm(alpha(0),alphao,tau,true);
  umsy(0)     = Type(.5) * alpha(0) - Type(0.07) * (alpha(0) * alpha(0));
  Smsy(0)     =  alpha(0)/beta * (Type(0.5) -Type(0.07) * alpha(0));  

  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm(alpha(i),alpha(i-1),tau,true);
  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logR(i))){
      logRS(i) = alpha(i) - beta * obs_S(i) ;
      pred_logR(i) = logRS(i) + log(obs_S(i));
      umsy(i) = Type(.5) * alpha(i) - Type(0.07) * (alpha(i) * alpha(i)); 
      Smsy(i) =  alpha(i)/beta * (Type(0.5) -Type(0.07) * alpha(i));
      ans+=-dnorm(obs_logR(i),pred_logR(i),sig,true);
    }
  
  }

  REPORT(pred_logR)
  REPORT(alpha)
  REPORT(sig)
  REPORT(tau)
  REPORT(rho)
  REPORT(beta)
  REPORT(varphi)
  REPORT(alphao)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  return ans;
}

