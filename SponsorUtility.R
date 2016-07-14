## SponsorUtility
## Contains implementations of the sponsor utility function.
## Versions for both the summation and integration approach to computing the stage 2 utility are defined here.

##
# DeltaSHatOneStage. Compute the estimate of the effect in the subgroup over one stage.
##
DeltaSHatOneStage=function(design,zS,sigma){
    return(zS*sqrt(sigma^2/design$nS));
}

##
# DeltaSHatTwoStages. Compute the estimate of the effect in the subgroup over both stages.
##
DeltaSHatTwoStages=function(design1,design2,zS1,zS2,sigma){  
    deltaSHat1=DeltaSHatOneStage(design1,zS1,sigma);
    deltaSHat2=DeltaSHatOneStage(design2,zS2,sigma);
    return((design1$nS*deltaSHat1+design2$nS*deltaSHat2)/(design1$nS+design2$nS));
}

##
# DeltaFHatOneStage. Compute the estimate of the effect in the full population over one stage.
##
DeltaFHatOneStage=function(design,zS,zSC,lambda,sigma){
    return(lambda*sqrt(sigma^2/design$nS)*zS+(1-lambda)*sqrt(sigma^2/design$nSC)*zSC);
}

##
# DeltaFHatTwoStages. Compute the estimate of the effect in the full population over both stages.
##
DeltaFHatTwoStages=function(design1,design2,zS1,zSC1,zS2,zSC2,lambda,sigma){
    n1=design1$nS+design1$nSC;
    n2=design2$nS+design2$nSC;    
    deltaFHat1=DeltaFHatOneStage(design1,zS1,zSC1,lambda,sigma);
    deltaFHat2=DeltaFHatOneStage(design2,zS2,zSC2,lambda,sigma);
    return((n1*deltaFHat1+n2*deltaFHat2)/(n1+n2));
}

##
# SUtility. Calculates the utility (sponsor) for the adaptive enrichment design based on all data available after the trial.
# This is the version to be used with the summation approach.
##

## Top level generic function.
## There are a total of nine combinations for possible stage 1 and stage 2 designs,
## since the Null, FullEnrichment and PartialEnrichment designs are available at each stage.
SUtility=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("SUtility",design1);
}

## Generic function for Null first stage design.
SUtility.Null=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("SUtility.Null",design2);
}

## Method for Null first stage design and FullEnrichment second stage design.
SUtility.Null.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS2){        
    deltaSHat=DeltaSHatOneStage(design2,zS2,params$sigma);    
    psiS=PsiS(zS1=0,zS2,wS=0,z_ahalf=qnorm(1-params$alpha));
   
    C=TrialCost(n1=0,n2=design2$nS,lambda1=0.5,lambda2=1,params);       
    U=params$r*params$N*params$lambda*psiS*pmax(deltaSHat-params$muS,0)-C;    
    return(U);
}

## Method for Null first stage design and PartialEnrichment second stage design.
SUtility.Null.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS2,zSC2){
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;
    
    deltaSHat=DeltaSHatOneStage(design2,zS2,params$sigma);
    deltaFHat=DeltaFHatOneStage(design2,zS2,zSC2,params$lambda,params$sigma);

    psiF=PsiF(zS1=0,zS2,zSC1=0,zSC2,lambda1=0.5,lambda2,wS=0,wSC=0,wF=0,params$lambda,params$z_ahalf,params$z_etaS,params$z_etaSC);
    psiS=PsiS(zS1=0,zS2,wS=0,params$z_ahalf);
    
    C=TrialCost(n1=0,n2,lambda1=0.5,lambda2,params);       
    U=params$r*params$N*(psiF*pmax(deltaFHat-params$muF,0)+params$lambda*psiS*(1-psiF)*pmax(deltaSHat-params$muS,0))-C;   
    return(U);
}

## Generic function for FullEnrichment first stage design.
SUtility.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("SUtility.FullEnrichment",design2);
}

## Method for FullEnrichment first stage design and FullEnrichment second stage design.
SUtility.FullEnrichment.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zS2){
    deltaSHat=DeltaSHatTwoStages(design1,design2,zS1,zS2,params$sigma);
    psiS=PsiS(zS1,zS2,params$wS,z_ahalf=qnorm(1-params$alpha));

    C=TrialCost(n1=design1$nS,n2=design2$nS,lambda1=1,lambda2=1,params);       
    U=params$r*params$N*params$lambda*psiS*pmax(deltaSHat-params$muS,0)-C;
    return(U);
}

## Method for FullEnrichment first stage design and PartialEnrichment second stage design.
## It is assumed that approval (and estimate) in the full population may only be based on zF2.
SUtility.FullEnrichment.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zS2,zSC2){
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;

    deltaSHat=DeltaSHatTwoStages(design1,design2,zS1,zS2,params$sigma);
    deltaFHat=DeltaFHatOneStage(design2,zS2,zSC2,params$lambda,params$sigma);
    
    zF2=ZF(zS2,zSC2,params$lambda,lambda2);
    psiF=InvNormCombTest(z1=0,z2=zF2,w=0,params$z_ahalf) & InvNormCombTest(zS1,zS2,params$wS,params$z_etaS) & InvNormCombTest(z1=0,z2=zSC2,w=0,params$z_etaSC);        
    psiS=PsiS(zS1,zS2,params$wS,params$z_ahalf);
    
    C=TrialCost(n1=design1$nS,n2,lambda1=1,lambda2,params);       
    U=params$r*params$N*(psiF*pmax(deltaFHat-params$muF,0)+params$lambda*psiS*(1-psiF)*pmax(deltaSHat-params$muS,0))-C;

    return(U);
}

## Generic function for PartialEnrichment first stage design.
SUtility.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("SUtility.PartialEnrichment",design2);
}

## Method for PartialEnrichment first stage design and FullEnrichment second stage design.
## It is assumed that H_F is retained, that is, there can only be approval in S given a
## first stage partial enrichment design followed by a full enrichment design.
## Only H_S is tested, so no alpha split needs to be done for this design.
SUtility.PartialEnrichment.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zSC1,zS2){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;

    deltaSHat=DeltaSHatTwoStages(design1,design2,zS1,zS2,params$sigma);     
    psiS=PsiS(zS1,zS2,params$wS,z_ahalf=qnorm(1-params$alpha));
    
    C=TrialCost(n1,n2=design2$nS,lambda1,lambda2=1,params);       
    U=params$r*params$N*params$lambda*psiS*pmax(deltaSHat-params$muS,0)-C;
    return(U);
}

## Method for PartialEnrichment first stage design and PartialEnrichment second stage design.
SUtility.PartialEnrichment.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zSC1,zS2,zSC2){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;
    
    deltaSHat=DeltaSHatTwoStages(design1,design2,zS1,zS2,params$sigma);
    deltaFHat=DeltaFHatTwoStages(design1,design2,zS1,zSC1,zS2,zSC2,params$lambda,params$sigma);   
    
    psiF=PsiF(zS1,zS2,zSC1,zSC2,lambda1,lambda2,params$wS,params$wSC,params$wF,params$lambda,params$z_ahalf,params$z_etaS,params$z_etaSC);
    psiS=PsiS(zS1,zS2,params$wS,params$z_ahalf);
    C=TrialCost(n1,n2,lambda1,lambda2,params);

    U=params$r*params$N*(psiF*pmax(deltaFHat-params$muF,0)+params$lambda*psiS*(1-psiF)*pmax(deltaSHat-params$muS,0))-C;   
    return(U);
}

############################################
##Integrals for the HF part of the utility##
############################################

IntegrandS3RF=function(zSC2,a,mS,mSC,xi2,DeltaF1,n1,n2,lam,lam2,cF,muF,sigma){
  invg=invg(zSC2,cF,xi2,lam,lam2);
  invk=invk(zSC2,DeltaF1,n1,n2,lam,lam2,muF,sigma);
  return(zSC2*dnorm(zSC2-mSC)*(1-pnorm(pmax(a,invg,invk)-mS)));
}

IntegrandS1RFplusS2RF=function(zS2,b,mS,mSC,xi2,DeltaF1,n1,n2,lam,lam2,cF,muF,sigma){
  g=g(zS2,cF,xi2,lam,lam2);
  k=k(zS2,DeltaF1,n1,n2,lam,lam2,muF,sigma);
  return(dnorm(zS2-mS)*(1-pnorm(pmax(b,g,k)-mSC))*(DeltaF1-muF+n2/(n1+n2)*zS2*zetaS2(lam,lam2,n2,sigma)));
}

S1RFplusS2RF=function(a,b,mS,mSC,xi2,DeltaF1,n1,n2,lam,lam2,cF,muF,sigma){
  return(integrate(IntegrandS1RFplusS2RF,a,Inf,b=b,mS=mS,mSC=mSC,xi2=xi2,DeltaF1=DeltaF1,n1=n1,n2=n2,lam=lam,lam2=lam2,cF=cF,muF=muF,sigma=sigma)$value)
}

S3RF=function(a,b,mS,mSC,xi2,DeltaF1,n1,n2,lam,lam2,cF,muF,sigma){
  return(n2/(n1+n2)*zetaSC2(lam,lam2,n2,sigma)*integrate(IntegrandS3RF,b,Inf,a=a,mS=mS,mSC=mSC,xi2=xi2,DeltaF1=DeltaF1,n1=n1,n2=n2,lam=lam,lam2=lam2,cF=cF,muF=muF,sigma=sigma)$value)
}

SRF=function(a,b,mS,mSC,xi2,DeltaF1,n1,n2,lam,lam2,cF,muF,sigma){
  return(S1RFplusS2RF(a,b,mS,mSC,xi2,DeltaF1,n1,n2,lam,lam2,cF,muF,sigma)+S3RF(a,b,mS,mSC,xi2,DeltaF1,n1,n2,lam,lam2,cF,muF,sigma))
}

############################################
##Integrals for the HS part of the utility##
############################################

#The evaluation strategy follows the manuscript Numeric Evaluations for Opimizint Adaptive Enrichment designs.
#This mean that the Utility for rejecting only HS i.e. via the Formula SRFCRS=SRS-SRFRS 

#Note: The following evaluation of can be done directly via calculation of an antiderivative
#Should muS really be included in the integrand?
IntegrandSRS=function(zS2,mS,DeltaS1,n1,n2,lam1,lam2,cS,muS,sigma){
  return(dnorm(zS2-mS)*(DeltaS1-muS+lam2*n2/(lam1*n1+lam2*n2)*sqrt(sigma^2/(lam2*n2))*zS2));
}

SRS=function(mS,kappa,DeltaS1,n1,n2,lam1,lam2,cS,muS,sigma){
  maxcSkappa=pmax(kappa,cS)
  return(integrate(IntegrandSRS,maxcSkappa,Inf,mS=mS,DeltaS1=DeltaS1,n1=n1,n2=n2,lam1=lam1,lam2=lam2,cS=cS,muS=muS,sigma=sigma)$value);
}

IntegrandSRFRS=function(zS2,mS,mSC,DeltaS1,cF,xi2,lam,lam1,lam2,n1,n2,b,muS,sigma){
  g=g(zS2,cF,xi2,lam,lam2);
  return(dnorm(zS2-mS)*(DeltaS1-muS+lam2*n2/(lam1*n1+lam2*n2)*sqrt(sigma^2/(lam2*n2))*zS2)*(1-pnorm(pmax(b,g)-mSC)));
}

SRFRS=function(a,kappa,cS,mS,mSC,DeltaS1,cF,xi2,lam,lam1,lam2,n1,n2,b,muS,sigma){
  M=pmax(a,kappa,cS);
  return(integrate(IntegrandSRFRS,M,Inf,mS=mS,mSC=mSC,DeltaS1=DeltaS1,cF=cF,xi2=xi2,lam=lam,lam1=lam1,lam2=lam2,n1=n1,n2=n2,b=b,muS=muS,sigma=sigma)$value)
}

SRFCRS=function(a,kappa,cS,mS,mSC,DeltaS1,cF,xi2,lam,lam1,lam2,n1,n2,b,muS,sigma){
  return(SRS(mS,kappa,DeltaS1,n1,n2,lam1,lam2,cS,muS,sigma)-SRFRS(a,kappa,cS,mS,mSC,DeltaS1,cF,xi2,lam,lam1,lam2,n1,n2,b,muS,sigma))
}

##
# SUtilityInt. Calculation of Sponsor Utility.
# This is the version to be used for the integration approach.
##

## Top level generic function.
SUtilityInt=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("SUtilityInt",design1);
}

## Generic function for Null first stage design.
SUtilityInt.Null=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("SUtilityInt.Null",design2);
}

## Method for Null first stage design and FullEnrichment second stage design.
SUtilityInt.Null.FullEnrichment=function(deltaS,deltaSC,design1,design2,params){
    n1=0;
    lambda1=0.5;
    n2=design2$nS;
    lambda2=1;

    mS=mS(deltaS,lambda2,n2,params$sigma);
    cS=cS(z_ahalf=qnorm(1-params$alpha),wS=0,zS1=0);    
    kappa=kappa(DeltaS1=0,n1,n2,lambda1,lambda2,params$muS,params$sigma);

    IRS=SRS(mS,kappa,DeltaS1=0,n1,n2,lambda1,lambda2,cS,params$muS,params$sigma);
    
    C=TrialCost(n1,n2,lambda1,lambda2,params);
    U=params$r*params$N*params$lambda*IRS-C;
    return(U);
}

## Method for Null first stage design and PartialEnrichment second stage design.
SUtilityInt.Null.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS2,zSC2){
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;

    mS=mS(deltaS,lambda2,n2,params$sigma);
    mSC=mSC(deltaSC,lambda2,n2,params$sigma);
    
    a=a(params$z_etaS,wS=0,zS1=0);        
    b=b(params$z_etaSC,wSC=0,zSC1=0);
    cS=cS(params$z_ahalf,wS=0,zS1=0);
    cF=cF(params$z_ahalf,wF=0,zF1=0); 

    xi2=Xi(params$lambda,lambda2);
    kappa=kappa(DeltaS1=0,n1=0,n2,lam1=0.5,lambda2,params$muS,params$sigma);
                
    C=TrialCost(n1=0,n2,lambda1=0.5,lambda2,params);
    U=params$r*params$N*(SRF(a,b,mS,mSC,xi2,DeltaF1=0,n1=0,n2,params$lambda,lambda2,cF,params$muF,params$sigma)+params$lambda*SRFCRS(a,kappa,cS,mS,mSC,DeltaS1=0,cF,xi2,params$lambda,lam1=0.5,lambda2,n1=0,n2,b,params$muS,params$sigma))-C;
    return(U);
}

## Generic function for FullEnrichment first stage design.
SUtilityInt.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("SUtilityInt.FullEnrichment",design2);
}

## Method for FullEnrichment first stage design and FullEnrichment second stage design
## Only H_S is tested, so no alpha split needs to be done for this design.
SUtilityInt.FullEnrichment.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1){
    n1=design1$nS;
    lambda1=1;
    n2=design2$nS;
    lambda2=1;

    mS=mS(deltaS,lambda2,n2,params$sigma);
    cS=cS(z_ahalf=qnorm(1-params$alpha),params$wS,zS1);
    DeltaS1=DeltaS1(zS1,n1,n2,lambda1,lambda2,params$sigma);
    kappa=kappa(DeltaS1,n1,n2,lambda1,lambda2,params$muS,params$sigma);

    IRS=SRS(mS,kappa,DeltaS1,n1,n2,lambda1,lambda2,cS,params$muS,params$sigma);
    
    C=TrialCost(n1,n2,lambda1,lambda2,params);
    U=params$r*params$N*params$lambda*IRS-C;
    return(U);
}

## Method for FullEnrichment first stage design and PartialEnrichment second stage design.
## It is assumed that approval in the full population may only be based on zF2.
SUtilityInt.FullEnrichment.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1){  
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;
    
    mS=mS(deltaS,lambda2,n2,params$sigma);
    mSC=mSC(deltaSC,lambda2,n2,params$sigma);   
    
    a=a(params$z_etaS,params$wS,zS1);        
    b=b(params$z_etaSC,wSC=0,zSC1=0);
    
    cS=cS(params$z_ahalf,params$wS,zS1);
    cF=cF(params$z_ahalf,wF=0,zF1=0);

    DeltaS1=DeltaS1(zS1,n1=design1$nS,n2,lam1=1,lambda2,params$sigma);

    xi2=Xi(params$lambda,lambda2);
    kappa=kappa(DeltaS1,n1=design1$nS,n2,lam1=1,lambda2,params$muS,params$sigma);
                
    C=TrialCost(n1=design1$nS,n2,lambda1=1,lambda2,params);
    U=params$r*params$N*(SRF(a,b,mS,mSC,xi2,DeltaF1=0,n1=0,n2,params$lambda,lambda2,cF,params$muF,params$sigma)+params$lambda*SRFCRS(a,kappa,cS,mS,mSC,DeltaS1,cF,xi2,params$lambda,lam1=1,lambda2,n1=design1$nS,n2,b,params$muS,params$sigma))-C;
    return(U);
}

## Generic function for PartialEnrichment first stage design.
SUtilityInt.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,...){
    UseMethod("SUtilityInt.PartialEnrichment",design2);
}

## Method for PartialEnrichment first stage design and FullEnrichment second stage design.
## It is assumed that H_F is retained, that is, there can only be approval in S given a
## first stage partial enrichment design followed by a full enrichment design.
## Only H_S is tested, so no alpha split needs to be done for this design.
SUtilityInt.PartialEnrichment.FullEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zSC1){ 
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
    n2=design2$nS;
    lambda2=1;

    mS=mS(deltaS,lambda2,n2,params$sigma);
    cS=cS(z_ahalf=qnorm(1-params$alpha),params$wS,zS1);
    DeltaS1=DeltaS1(zS1,n1,n2,lambda1,lambda2,params$sigma);
    kappa=kappa(DeltaS1,n1,n2,lambda1,lambda2,params$muS,params$sigma);

    IRS=SRS(mS,kappa,DeltaS1,n1,n2,lambda1,lambda2,cS,params$muS,params$sigma);
    
    C=TrialCost(n1,n2,lambda1,lambda2,params);
    U=params$r*params$N*params$lambda*IRS-C;
    
    return(U);   
}

## Method for PartialEnrichment first stage design and PartialEnrichment second stage design.
SUtilityInt.PartialEnrichment.PartialEnrichment=function(deltaS,deltaSC,design1,design2,params,zS1,zSC1){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;    

    mS=mS(deltaS,lambda2,n2,params$sigma);
    mSC=mSC(deltaSC,lambda2,n2,params$sigma);

    zF1=ZF(zS1,zSC1,params$lambda,lambda1);
    
    a=a(params$z_etaS,params$wS,zS1);        
    b=b(params$z_etaSC,params$wSC,zSC1);
    cS=cS(params$z_ahalf,params$wS,zS1);
    cF=cF(params$z_ahalf,params$wF,zF1);

    DeltaS1=DeltaS1(zS1,n1,n2,lambda1,lambda2,params$sigma);
    DeltaF1=DeltaF1(zS1,zSC1,n1,n2,lambda1,lambda2,params$lambda,params$sigma);

    xi2=Xi(params$lambda,lambda2);
    kappa=kappa(DeltaS1,n1,n2,lambda1,lambda2,params$muS,params$sigma);
                
    C=TrialCost(n1,n2,lambda1,lambda2,params);
    U=params$r*params$N*(SRF(a,b,mS,mSC,xi2,DeltaF1,n1,n2,params$lambda,lambda2,cF,params$muF,params$sigma)+params$lambda*SRFCRS(a,kappa,cS,mS,mSC,DeltaS1,cF,xi2,params$lambda,lambda1,lambda2,n1,n2,b,params$muS,params$sigma))-C;
    return(U);
}
