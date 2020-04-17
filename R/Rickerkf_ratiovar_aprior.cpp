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
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(obs_logR);   // observed log recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  //DATA_SCALAR(eps);
  DATA_SCALAR(prbeta1);
  DATA_SCALAR(prbeta2);
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value
  //rho         -> Proportion of total variance associated with obs error.
  //varphi      -> Total precision
  //alpha       -> Time-varying alpha

  PARAMETER(logalphao);
  PARAMETER(logbeta);
  PARAMETER(rho);
  PARAMETER(logvarphi);

  PARAMETER_VECTOR(logalpha);
  

  
  int timeSteps=obs_logR.size();

  Type beta=exp(logbeta);
  Type Smax  = Type(1.0)/beta;
  
  //theta       -> total standard deviation
  //sig         -> obs error std
  //tau         -> proc error (alpha) std
  
  Type varphi     = exp(logvarphi);
  Type theta     = sqrt(Type(1.0)/varphi);
  Type sig       = sqrt(rho) * theta;
  Type tau        = sqrt(Type(1.0)-rho) * theta ;

  Type pen;
  pen=0;
  Type eps;
  eps=1e-3;


  vector<Type> pred_logR(timeSteps), logRS(timeSteps),umsy(timeSteps);

  

  //priors on precision and variance ratio
  Type ans= -dbeta(rho,prbeta1,prbeta2,true);  
   
  
  ans+= -dnorm(exp(logalpha(0)),exp(logalphao),tau,true);
  umsy(0)     = Type(.5) * exp(logalpha(0)) - Type(0.07) * (exp(2*logalpha(0)) ); 

  //Type ans= -dnorm(alpha(0),alphao,tau,true); 
  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm(exp(logalpha(i)),exp(logalpha(i-1)),tau,true);
  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logR(i))){
      logRS(i) = exp(logalpha(i)) - beta * obs_S(i) ;
      pred_logR(i) = logRS(i) + log(obs_S(i));
      umsy(i)     = Type(.5) * exp(logalpha(i)) - Type(0.07) * (exp(2*logalpha(i) )); 


      
      ans+=-dnorm(obs_logR(i),pred_logR(i),sig,true);


    }
  
  }

  REPORT(pred_logR)
  REPORT(logalpha)
  REPORT(sig)
  REPORT(tau)
  REPORT(rho)
  REPORT(beta)
  REPORT(varphi)
  REPORT(logalphao)
  REPORT(Smax)
  REPORT(umsy)
  return ans;
}

