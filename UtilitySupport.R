## UtilitySupport

## Contains some support functions for defining the utility functions for the
## sponsor and public health view.

## Notation

# z_etaS...z-value for the reduced level etaS, i.e. z_etaS=qnorm(1-etaS)
# z_etaSC...z-value for the reduced level etaSC, i.e. z_etaS=qnorm(1-etaSC)
# z_sig_level...z-value for corresponding to a given significance level
# wS...weight for the first stage data for HS
# wSC...weight for the first stage complement data (data from first stage not belonging to S)
# wF...weight for the first stage data for HF
# zS1...the first stage z-value for the subgroup S
# zS2...the second stage z-value for the subgroup S
# zSC1...the first stage z-value for the subgroup SC
# zSC2...the second stage z-value for the subgroup SC
# lam...the prevalence of the subgroup
# lam1...first stage trial prevalence
# lam2...second stage trial prevalence
# cF,cS...can be interpreted as a significance levels adjusted for the the interim observations (for details see manuscript)
# n1...first stage sample size
# n2...second stage sample size
# sigma^2...equals the sum sigma_T^2+sigma_C^2, which are the variances of
#           the treatment and the control (these must not depend on the biomarker status of the patient)


# InvNormComb. Inverse normal combination function.
InvNormComb=function(x,y,w){
    return(sqrt(w)*x+sqrt(1-w)*y);
}

##
# InvNormCombTest. Inverse normal combination function test.
# z1 = Stage 1 Z-statistic.
# z2 = Stage 2 Z-statistic.
# w = Weight for inverse normal combination.
# q = Z-value corresponding to some one-sided significance level.
##
InvNormCombTest=function(z1,z2,w,q){
    return(InvNormComb(z1,z2,w)>q);
}

##
# Xi. The factor used in the definition of the enrichment adjusted test statistic zF.
# lambda = Prevalence in the target population.
# trialLambda = Prevalence in the enriched trial (either in stage 1 or stage 2).
##
Xi=function(lambda,trialLambda){
    return(1/sqrt(lambda^2/trialLambda+(1-lambda)^2/(1-trialLambda)));
}

# ZF. Compute zF given zS and zSC.
ZF=function(zS,zSC,lambda,trialLambda){
    return(Xi(lambda,trialLambda)*((lambda/sqrt(trialLambda))*zS+((1-lambda)/sqrt(1-trialLambda))*zSC));    
}

##
# PsiS. Test for HS, based on the first and second stage z-values zS1 and zS2.
# Returns TRUE if the test rejects and FALSE otherwise.
# wS = Weight for inverse normal combination.
# z_sig_level = Z-value corresponding to a one-sided significance level of alpha / 2.
##
PsiS=function(zS1,zS2,wS,z_sig_level){
    return(InvNormCombTest(zS1,zS2,wS,z_sig_level));
}

##
# PsiF. Test for HF, based on the first and second stage data.
# Returns TRUE if the test rejects and FALSE otherwise.
##
PsiF=function(zS1,zS2,zSC1,zSC2,lambda1,lambda2,wS,wSC,wF,lambda,z_sig_level,z_etaS,z_etaSC){    
    zF1=ZF(zS1, zSC1, lambda, lambda1);
    zF2=ZF(zS2, zSC2, lambda, lambda2);
    
    return(InvNormCombTest(zF1,zF2,wF,z_sig_level) &
           InvNormCombTest(zS1,zS2,wS,z_etaS) &
           InvNormCombTest(zSC1,zSC2,wSC,z_etaSC));  
}

# ScreeningCost. Screening costs for a single stage.
ScreeningCost=function(lambda,trialLambda,n,cScreening){
    return(2*cScreening*n*max(trialLambda/lambda,(1-trialLambda)/(1-lambda)));        
}

# TrialCost. Total cost over both stages of performing the trial.
TrialCost=function(n1,n2,lambda1,lambda2,params){
    return(params$cSetup+params$cBiomarker+2*params$cPerPatient*(n1+n2)+ScreeningCost(params$lambda,lambda1,n1,params$cScreening)+ScreeningCost(params$lambda,lambda2,n2,params$cScreening));
}

# Returns: The variable a (for definition see manuscript)
a=function(z_etaS,wS,zS1){
    return((z_etaS - sqrt(wS)*zS1)/(sqrt(1-wS)));
}

# Returns: The variable b (for definition see manuscript)
b=function(z_etaSC,wSC,zSC1){
    return((z_etaSC - sqrt(wSC)*zSC1)/(sqrt(1-wSC)));
}

# Returns: The variable cF (for definition see manuscript)
cF=function(z_sig_level,wF,zF1){
    return((z_sig_level-sqrt(wF)*zF1)/sqrt(1-wF));
}

# Returns: The variable cS (for definition see manuscript)
cS=function(z_sig_level,wS,zS1){
    return((z_sig_level-sqrt(wS)*zS1)/sqrt(1-wS));
}

# Returns: The value zSC2 needed to reject HF as a function of zS1 (for details see manuscript)
g=function(zS2,cF,xi2,lam,lam2){
    return((cF/xi2-lam/sqrt(lam2)*zS2)*sqrt(1-lam2)/(1-lam));
}

# Returns: The inverse of function g (as a function of zSC1, for details see manuscript)
invg=function(zSC2,cF,xi2,lam,lam2){
    return((cF/xi2-(1-lam)/sqrt(1-lam2)*zSC2)*sqrt(lam2)/lam);
}

# Returns: The value zSC2 needed to obtain an estimate for deltaF larger than muF (for details see manuscript)
k=function(zS2,DeltaF1,n1,n2,lam,lam2,muF,sigma){
    return(((muF-DeltaF1)*(n1+n2)/n2 - lam*sqrt(sigma^2/(lam2*n2))*zS2)*sqrt((1-lam2)*n2/(sigma^2*(1-lam)^2)));
}

# Returns: The inverse of function k (as a function of zSC1, for details see manuscript)
invk=function(zSC2,DeltaF1,n1,n2,lam,lam2,muF,sigma){
    ((muF-DeltaF1)*(n1+n2)/n2-(1-lam)*sqrt(sigma^2/((1-lam2)*n2))*zSC2)*sqrt(lam2*n2/(sigma^2*lam^2));
}

# Returns zetaS2 (for details see manuscript)
zetaS2=function(lam,lam2,n2,sigma){
    return(lam*sqrt(sigma^2/(lam2*n2)))
}

# Returns zetaSC2 (for details see manuscript)
zetaSC2=function(lam,lam2,n2,sigma){
    return((1-lam)*sqrt(sigma^2/((1-lam2)*n2)))
}

# Returns: The value zS2 needed to obtain an estimate for deltaS larger than muS
kappa=function(DeltaS1,n1,n2,lam1,lam2,muS,sigma){
    return((muS-DeltaS1)*(lam1*n1+lam2*n2)/(lam2*n2)*sqrt(lam2*n2/sigma^2));
}

# The mean of the random variable zS2 if the true effect equals deltaS and the trial parameters are n2 and lam2
mS=function(deltaS,lam2,n2,sigma){
    return(deltaS*sqrt(lam2*n2/sigma^2));
}

# The mean of the random variable zSC2 if the true effect equals deltaSC and the trial parameters are n2 and lam2
mSC=function(deltaSC,lam2,n2,sigma){
    return(deltaSC*sqrt((1-lam2)*n2/sigma^2));
}

# Returns: The variable DeltaS1, for details see manuscript
DeltaS1=function(zS1,n1,n2,lam1,lam2,sigma){
    return(lam1*n1/(lam1*n1+lam2*n2)*sqrt(sigma^2/(lam1*n1))*zS1); 
}

# Returns: The variable DeltaF1, for details see manuscript
DeltaF1=function(zS1,zSC1,n1,n2,lam1,lam2,lam,sigma){
    return(n1/(n1+n2)*( lam*sqrt(sigma^2/(lam1*n1))*zS1 + (1-lam)*sqrt(sigma^2/((1-lam1)*n1))*zSC1 ) );  
}


