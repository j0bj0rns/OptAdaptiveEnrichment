## PublicHealthUtility
## Contains implementations of the public health utility function.
## Versions for both the summation and integration approach to computing the stage 2 utility are defined here.

##
# PHUtility. Calculates the utility for the adaptive enrichment design based on all data available after the trial.
# This version should be used with the summation approach.
##

## Top level generic function.
## There are a total of six combinations for possible stage 1 and stage 2 designs,
## since the Null, FullEnrichment and PartialEnrichment designs are available at the first stage,
## and the FullEnrichment and PartialEnrichment designs are available at the second stage.
## Note that there is no possibility of early stopping after the first stage.
## Implementing this would require a generalisation of the testing procedure, since stage-wise multiplicity would be introduced.
PHUtility=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("PHUtility",design1);
}

## Generic function for Null first stage design.
PHUtility.Null=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("PHUtility.Null",design2);
}

## Method for Null first stage design and FullEnrichment second stage design (i.e., a one-stage FullEnrichment design).
## Only H_S is tested, so no alpha split needs to be done for this design.
PHUtility.Null.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS2){
    psiS=PsiS(zS1=0,zS2,wS=0,z_ahalf=qnorm(1-params$alpha));
    
    C=TrialCost(n1=0,n2=design2$nS,lambda1=0.5,lambda2=1,params);       
    U=params$r*params$N*params$lambda*psiS*(deltaS-params$muS)-C;    
    return(U);
}

## Method for Null first stage design and PartialEnrichment second stage design.
PHUtility.Null.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS2,zSC2){
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;
    
    deltaF=params$lambda*deltaS+(1-params$lambda)*deltaSC;
    psiF=PsiF(zS1=0,zS2,zSC1=0,zSC2,lambda1=0.5,lambda2,wS=0,wSC=0,wF=0,params$lambda,params$z_ahalf,params$z_etaS,params$z_etaSC);
    psiS=PsiS(zS1=0,zS2,wS=0,params$z_ahalf);
    
    C=TrialCost(n1=0,n2,lambda1=0.5,lambda2,params);       
    U=params$r*params$N*(psiF*(deltaF-params$muF)+params$lambda*psiS*(1-psiF)*(deltaS-params$muS))-C;    
    return(U);
}

## Generic function for FullEnrichment first stage design.
PHUtility.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("PHUtility.FullEnrichment",design2);
}

## Method for FullEnrichment first stage design and FullEnrichment second stage design.
## Only H_S is tested, so no alpha split needs to be done for this design.
PHUtility.FullEnrichment.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zS2){
    psiS=PsiS(zS1,zS2,params$wS,z_ahalf=qnorm(1-params$alpha));

    C=TrialCost(n1=design1$nS,n2=design2$nS,lambda1=1,lambda2=1,params);       
    U=params$r*params$N*params$lambda*psiS*(deltaS-params$muS)-C;    
    return(U);
}

## Method for FullEnrichment first stage design and PartialEnrichment second stage design.
## For this case, it is assumed that only the second stage data can be used for deciding on approval in F.
PHUtility.FullEnrichment.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zS2,zSC2){
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;
    
    deltaF=params$lambda*deltaS+(1-params$lambda)*deltaSC;    
    zF2=ZF(zS2,zSC2,params$lambda,lambda2);
  
    psiF=InvNormCombTest(z1=0,z2=zF2,w=0,params$z_ahalf) & InvNormCombTest(zS1,zS2,params$wS,params$z_etaS) & InvNormCombTest(z1=0,z2=zSC2,w=0,params$z_etaSC);        
    psiS=PsiS(zS1,zS2,params$wS,params$z_ahalf);
    
    C=TrialCost(n1=design1$nS,n2,lambda1=1,lambda2,params);       
    U=params$r*params$N*(psiF*(deltaF-params$muF)+params$lambda*psiS*(1-psiF)*(deltaS-params$muS))-C;    
    return(U);
}


## Generic function for PartialEnrichment first stage design.
PHUtility.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("PHUtility.PartialEnrichment",design2);
}

## Method for PartialEnrichment first stage design and FullEnrichment second stage design.
## It is assumed that H_F is retained, that is, there can only be approval in S given a
## first stage partial enrichment design followed by a full enrichment design.
## Only H_S is tested, so no alpha split needs to be done for this design.
PHUtility.PartialEnrichment.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zSC1,zS2){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
         
    psiS=PsiS(zS1,zS2,params$wS,z_ahalf=qnorm(1-params$alpha));
    
    C=TrialCost(n1,n2=design2$nS,lambda1,lambda2=1,params);       
    U=params$r*params$N*params$lambda*psiS*(deltaS-params$muS)-C;    
    return(U);
}

## Method for PartialEnrichment first stage design and PartialEnrichment second stage design.
PHUtility.PartialEnrichment.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zSC1,zS2,zSC2){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;
    
    deltaF=params$lambda*deltaS+(1-params$lambda)*deltaSC;
    psiF=PsiF(zS1,zS2,zSC1,zSC2,lambda1,lambda2,params$wS,params$wSC,params$wF,params$lambda,params$z_ahalf,params$z_etaS,params$z_etaSC);
    psiS=PsiS(zS1,zS2,params$wS,params$z_ahalf);
    
    C=TrialCost(n1,n2,lambda1,lambda2,params);       
    U=params$r*params$N*(psiF*(deltaF-params$muF)+params$lambda*psiS*(1-psiF)*(deltaS-params$muS))-C;    
    return(U);   
}

##########
# Code for calculating the conditional expected utility from the public health view given a first stage observation zS1 and zSC1. 
# Notation follows the manuscript "Numeric Evaluations for "Opimizing Adaptive Enrichment Designs" " (Version 17_06_2016)
# Note that sigma^2 is the same as sigma_T^2+sigma_C^2 in the manuscript.
##########

########################################
## Power for rejecting HF             ##
########################################
Integrand_IRF=function(zS2,cF,xi2,lam,lam2,mS,mSC,b){
  g=g(zS2,cF,xi2,lam,lam2);
  return(dnorm(zS2-mS)*(1-pnorm(pmax(b,g)-mSC)));
}

IRF=function(a,cF,xi2,lam,lam2,mS,mSC,b){
  return(integrate(Integrand_IRF,lower=a,upper=Inf,rel.tol=.Machine$double.eps^0.5,cF=cF,xi2=xi2,lam=lam,lam2=lam2,mS=mS,mSC=mSC,b=b)$value)
}

########################################
## Power for rejecting exactly HS     ##
########################################

#Note that the integrand for the Integral IRSRF is the same as the integrand for IRF! Therefore no separate Integrand_IRSRF function is needed.
IRSRF=function(a,cF,xi2,lam,lam2,mS,mSC,b,cS){
  m=pmax(cS,a);
  return(integrate(Integrand_IRF,lower=m,upper=Inf,rel.tol=.Machine$double.eps^0.5,cF=cF,xi2=xi2,lam=lam,lam2=lam2,mS=mS,mSC=mSC,b=b)$value)
}

#Note that the Integral IRS can be evaluated directly, allowing to calculate IRFCRS=IRS-IRSRF with only one call of integrate()
IRFCRS=function(a,cF,xi2,lam,lam2,mS,mSC,b,cS){
  IRS=1-pnorm(cS-mS);
  return(IRS-IRSRF(a,cF,xi2,lam,lam2,mS,mSC,b,cS));
}

##
# PHUtilityInt. Calculation of Public Health Utility.
# This is the version to be used for the integration approach.
##

## Top level generic function.
PHUtilityInt=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("PHUtilityInt",design1);
}

## Generic function for Null first stage design.
PHUtilityInt.Null=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("PHUtilityInt.Null",design2);
}

## Method for Null first stage design and FullEnrichment second stage design.
## Only H_S is tested, so no alpha split needs to be done for this design.
PHUtilityInt.Null.FullEnrichment=function(deltaS,deltaSC,design1,design2,params){
    mS=mS(deltaS,lam2=1,n2=design2$nS,params$sigma);
    cS=cS(z_ahalf=qnorm(1-params$alpha),wS=0,zS1=0);  
    IRS=1-pnorm(cS-mS);
        
    C=TrialCost(n1=0,n2=design2$nS,lambda1=0.5,lambda2=1,params);    
    U=params$r*params$N*params$lambda*IRS*(deltaS-params$muS)-C;
    return(U);
}

## Method for Null first stage design and PartialEnrichment second stage design.
PHUtilityInt.Null.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params){
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;

    mS=mS(deltaS,lambda2,n2,params$sigma);
    mSC=mSC(deltaSC,lambda2,n2,params$sigma);
    deltaF=params$lambda*deltaS+(1-params$lambda)*deltaSC;   
    a=a(params$z_etaS,wS=0,zS1=0);        
    b=b(params$z_etaSC,wSC=0,zSC1=0);
    cS=cS(params$z_ahalf,wS=0,zS1=0);
    cF=cF(params$z_ahalf,wF=0,zF1=0);
    xi2=Xi(params$lambda,lambda2);
    
    C=TrialCost(n1=0,n2,lambda1=0.5,lambda2,params);    
    
    U=params$r*params$N*(IRF(a,cF,xi2,params$lambda,lambda2,mS,mSC,b)*(deltaF-params$muF) + params$lambda*IRFCRS(a,cF,xi2,params$lambda,lambda2,mS,mSC,b,cS)*(deltaS-params$muS))-C;
    return(U);   
}

## Generic function for FullEnrichment first stage design.
PHUtilityInt.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("PHUtilityInt.FullEnrichment",design2);
}

## Method for FullEnrichment first stage design and FullEnrichment second stage design
## Only H_S is tested, so no alpha split needs to be done for this design.
PHUtilityInt.FullEnrichment.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1){
    mS=mS(deltaS,lam2=1,n2=design2$nS,params$sigma);
    cS=cS(z_ahalf=qnorm(1-params$alpha),params$wS,zS1);
  
    IRS=1-pnorm(cS-mS);

    C=TrialCost(n1=design1$nS,n2=design2$nS,lambda1=1,lambda2=1,params);    
    U=params$r*params$N*params$lambda*IRS*(deltaS-params$muS)-C;
    return(U);
}

## Method for FullEnrichment first stage design and PartialEnrichment second stage design.
## It is assumed that approval in the full population may only be based on zF2.
PHUtilityInt.FullEnrichment.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1){
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;

    mS=mS(deltaS,lambda2,n2,params$sigma);
    mSC=mSC(deltaSC,lambda2,n2,params$sigma);
    deltaF=params$lambda*deltaS+(1-params$lambda)*deltaSC;
    a=a(params$z_etaS,params$wS,zS1);        
    b=b(params$z_etaSC,wSC=0,zSC1=0);
    cS=cS(params$z_ahalf,params$wS,zS1);
    cF=cF(params$z_ahalf,wF=0,zF1=0);    
    xi2=Xi(params$lambda,lambda2);
        
    C=TrialCost(n1=design1$nS,n2,lambda1=1,lambda2,params);       
    U=params$r*params$N*(IRF(a,cF,xi2,params$lambda,lambda2,mS,mSC,b)*(deltaF-params$muF) + params$lambda*IRFCRS(a,cF,xi2,params$lambda,lambda2,mS,mSC,b,cS)*(deltaS-params$muS))-C;           
    return(U);
}

## Generic function for PartialEnrichment first stage design.
PHUtilityInt.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("PHUtilityInt.PartialEnrichment",design2);
}

## Method for PartialEnrichment first stage design and FullEnrichment second stage design.
## It is assumed that H_F is retained, that is, there can only be approval in S given a
## first stage partial enrichment design followed by a full enrichment design.
## Only H_S is tested, so no alpha split needs to be done for this design.
PHUtilityInt.PartialEnrichment.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zSC1){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;       

    mS=mS(deltaS,lam2=1,n2=design2$nS,params$sigma);
    cS=cS(z_ahalf=qnorm(1-params$alpha),params$wS,zS1);
    IRS=1-pnorm(cS-mS);
  
    C=TrialCost(n1,n2=design2$nS,lambda1,lambda2=1,params);
    U=params$r*params$N*params$lambda*IRS*(deltaS-params$muS)-C;    
    return(U);
}

## Method for PartialEnrichment first stage design and PartialEnrichment second stage design.
PHUtilityInt.PartialEnrichment.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zSC1){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;
    
    mS=mS(deltaS,lambda2,n2,params$sigma);
    mSC=mSC(deltaSC,lambda2,n2,params$sigma);
    deltaF=params$lambda*deltaS+(1-params$lambda)*deltaSC;
  
    zF1=ZF(zS1,zSC1,params$lambda,lambda1);
    
    a=a(params$z_etaS,params$wS,zS1);        
    b=b(params$z_etaSC,params$wSC,zSC1);
    cS=cS(params$z_ahalf,params$wS,zS1);
    cF=cF(params$z_ahalf,params$wF,zF1);
    
    xi2=Xi(params$lambda,lambda2);
      
    C=TrialCost(n1,n2,lambda1,lambda2,params);    
    U=params$r*params$N*(IRF(a,cF,xi2,params$lambda,lambda2,mS,mSC,b)*(deltaF-params$muF) + params$lambda*IRFCRS(a,cF,xi2,params$lambda,lambda2,mS,mSC,b,cS)*(deltaS-params$muS))-C;    
    return(U);
}
