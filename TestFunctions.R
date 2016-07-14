## TestFunctions

## This file contains a collection of test functions and other examples.

## Check the implementation of GridHalfWidth using simulation.
Test.GridHalfWidth=function() {
    epsilon=0.9;
    hw=GridHalfWidth(epsilon);

    nSims=1000000;
    zS=rnorm(nSims);
    zSC=rnorm(nSims);

    cat("epsilon =",epsilon,"\n");    
    cat("Simulation epsilon =", sum(zS < -hw | zS > hw | zSC < -hw | zSC > hw) / nSims, "\n");
}

##### Check stage 2 summation results against the stage 2 integration results.
## stage1Type = "Null", "FullEnrichment" or "PartialEnrichment".
## stage2Type = "FullEnrichment" or "PartialEnrichment".

Test.Stage2EU=function(stage1Type,stage2Type){
    ## Set the parameter values
    params=Parameters(
        lambda=0.5,
        alpha=0.025,etaS=0.3,etaSC=0.3,
        wF=0.5,wS=0.5,wSC=0.5,
        N=10^6,r=10^4,
        cSetup=10^6,cBiomarker=10^6,cPerPatient=10^5,cScreening=10^5,
        muS=0.1,muF=0.1,
        sigma=sqrt(2));

    ## Construct the weak biomarker prior with dS=0.3
    dS=0.3;
    prior=cbind(c(0,dS,dS,dS),c(0,0,dS/2,dS),c(0.2,0.2,0.3,0.3));       

    epsilon1=10^-4;
    stepSize1=0.1;
    epsilon2=10^-4;
    stepSize2=0.02;

    ## Define the fixed stage 1 design and stage 1 data
    design1=0;
    zS1=0;
    zSC1=0;
    if (stage1Type=="Null") {
        design1=Null();
    }   

    if (stage1Type=="FullEnrichment") {
        design1=FullEnrichment(nS=100);
        zS1=0.1;        
    }

    if (stage1Type=="PartialEnrichment") {
        design1=PartialEnrichment(nS=100,nSC=100);
        zS1=0.1;
        zSC1=-0.1;
    }
        
    stage1Designs=list(design1);

    ## Define the set of possible two-stage designs
    n2S=seq(100,500,20);
    n2SC=200;
    
    stage2Designs=list();
    if (stage2Type=="FullEnrichment") {
        stage2Designs=lapply(n2S,function(x) FullEnrichment(x));
    }    

    if (stage2Type=="PartialEnrichment") {
        stage2Designs=lapply(n2S,function(x) PartialEnrichment(x,n2SC));
    }
    
    ## Construct a decision problem for the Public health view (summation approach)
    dpPHSum=DecisionProblemSum(prior=prior,utility=PHUtility,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        epsilon1=epsilon1,stepSize1=stepSize1,
        epsilon2=epsilon2,stepSize2=stepSize2);

    ## Construct the same decision problem but use the Sponsor view (summation approach)
    dpSSum=DecisionProblemSum(prior=prior,utility=SUtility,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        epsilon1=epsilon1,stepSize1=stepSize1,
        epsilon2=epsilon2,stepSize2=stepSize2);
    
    ## Construct a decision problem for the Public health view (integration approach)
    dpPHInt=DecisionProblemInt(prior=prior,utility=PHUtilityInt,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        epsilon1=epsilon1,stepSize1=stepSize1);
    
    ## Construct the same decision problem but use the Sponsor view (integration approach)
    dpSInt=DecisionProblemInt(prior=prior,utility=SUtilityInt,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        epsilon1=epsilon1,stepSize1=stepSize1);
    
    ## Compute Posterior    
    post=if (stage1Type=="Null") {
        ComputePosterior(design1,prior);
    } else if (stage1Type=="FullEnrichment") {
        ComputePosterior(design1,prior,params$sigma,zS1);
    } else if (stage1Type=="PartialEnrichment") {
        ComputePosterior(design1,prior,params$sigma,zS1,zSC1);
    }

    ## Compute the 2nd stage expected utilities using the summation approach
    s2PHSumEUs=0;
    s2SSumEUs=0;

    if (stage2Type=="FullEnrichment") {
         ## Create grid for zS2 (For FullEnrichment in second stage)
        hw=GridHalfWidth(epsilon2); ## Half width of grid square
        Gz2=seq(-hw,hw,stepSize2);

        if (stage1Type=="Null") {
            s2PHSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpPHSum,params,post,Gz2));
            s2SSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpSSum,params,post,Gz2));
        }
        else if (stage1Type=="FullEnrichment") {
            s2PHSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpPHSum,params,post,Gz2,zS1));
            s2SSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpSSum,params,post,Gz2,zS1));
        }
        else if (stage1Type=="PartialEnrichment") {
            s2PHSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpPHSum,params,post,Gz2,zS1,zSC1));
            s2SSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpSSum,params,post,Gz2,zS1,zSC1));
        }
        
        
    } else if (stage2Type=="PartialEnrichment") {
        ## Create grid for (zS2,zSC2) (For PartialEnrichment in second stage)
        hw=GridHalfWidth(epsilon2); ## Half width of grid square    
        Gz2=as.matrix(expand.grid(seq(-hw,hw,stepSize2),seq(-hw,hw,stepSize2)));
        colnames(Gz2)=NULL;

        if (stage1Type=="Null") {
            s2PHSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpPHSum,params,post,Gz2));
            s2SSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpSSum,params,post,Gz2));
        }
        else if (stage1Type=="FullEnrichment") {
            s2PHSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpPHSum,params,post,Gz2,zS1));
            s2SSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpSSum,params,post,Gz2,zS1));
        }
        else if (stage1Type=="PartialEnrichment") {
            s2PHSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpPHSum,params,post,Gz2,zS1,zSC1));
            s2SSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpSSum,params,post,Gz2,zS1,zSC1));
        }               
    }
    
    ## Compute the 2nd stage expecteded utilities using the integration approach       
    s2PHIntEUs=0;
    s2SIntEUs=0;

    if (stage1Type=="Null") {
        s2PHIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpPHInt,params,post));
        s2SIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpSInt,params,post));
    }
    else if (stage1Type=="FullEnrichment") {
        s2PHIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpPHInt,params,post,zS1));
        s2SIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpSInt,params,post,zS1));
    }
    else if (stage1Type=="PartialEnrichment") {
        s2PHIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpPHInt,params,post,zS1,zSC1));
        s2SIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpSInt,params,post,zS1,zSC1));
    }
       
    plot(c(100, 500),c(1*10^9, -1*10^9),type="n",xlab="n2S",ylab="Stage 2 EU")
    lines(n2S,s2PHSumEUs,lty=1,col="black");
    lines(n2S,s2PHIntEUs,lty=2,col="black");
    lines(n2S,s2SSumEUs,lty=1,col="red");
    lines(n2S,s2SIntEUs,lty=2,col="red");
}

Test.OptStage2EU=function(){
    ## Set the parameter values
    params=Parameters(
        lambda=0.5,
        alpha=0.025,etaS=0.3,etaSC=0.3,
        wF=0.5,wS=0.5,wSC=0.5,
        N=10^6,r=10^4,
        cSetup=10^6,cBiomarker=10^6,cPerPatient=10^5,cScreening=0,
        muS=0.1,muF=0.1,
        sigma=1);

    ## Construct three point prior for (deltaS, deltaSC)
    dS=0.3;
    prior=cbind(c(dS,dS,dS),c(0,dS/2,dS),c(1/3,1/3,1/3));
  
    ## Construct a decision problem    
    dp=DecisionProblem(
        prior=prior,
        utility=PHUtility,
        n1s=c(100),lambda1s=c(0.5),epsilon1=10^-3,stepSize1=1,
        n2s=seq(300,600,10),lambda2s=seq(0.05,0.95,0.025),epsilon2=10^-3,stepSize2=0.2);
    
    ## Set stage 1 choices and data
    n1=100;
    lambda1=0.5;
    zS1=0;
    zSC1=0;           

    cat("Integration approach gives:\n")
    print(NumOptFastv1(v,n1,lambda1,zS1,zSC1,params$lambda,params$N,params$r,cFix=params$cSetup,params$cPerPatient,params$cScreening,params$wF,params$wS,params$wSC,params$sigma));
    
    cat("Summation approach gives:\n")
    print(OptStage2EU(dp,params,n1,lambda1,zS1,zSC1));    
}

## Test SolveDP
Test.SolveDP=function(){
    ## Set the parameter values
    params=Parameters(
        lambda=0.5,
        alpha=0.025,etaS=0.3,etaSC=0.3,
        wF=0.5,wS=0.5,wSC=0.5,
        N=10^6,r=10^4,
        cSetup=10^6,cBiomarker=10^6,cPerPatient=10^5,cScreening=10^3,
        muS=0.1,muF=0.1,
        sigma=1);

    ## Construct three point prior for (deltaS, deltaSC)
    dS=0.3;
    prior=cbind(c(dS,dS,dS),c(0,dS/2,dS),c(1/3,1/3,1/3));
  
    ## Construct a decision problem    
    dp=DecisionProblem(
        prior=prior,
        utility=PHUtility,
        n1s=c(100),lambda1s=c(0.5),epsilon1=10^-3,stepSize1=0.2,
        n2s=seq(100,1000,50),lambda2s=seq(0.1,0.9,0.05),epsilon2=10^-3,stepSize2=0.2);    

    policy=SolveDP(dp,params);

    cat("Stage 1 optimal EU =",policy$eu1,"\n")
    cat("Stage 1 optimal sample size =",policy$n1,"\n")
    cat("Stage 1 optimal prevalence =",policy$lambda1,"\n")      
    
    PlotPolicy(policy,zSC1s=list(-2,0,2));
}

## Use the symmetry over the two stages to check the backward induction algorithm.
## The numerical result of two decision problems are compared.
## In problem I, a sample size and prevalence are chosen in stage 2, while a fixed sample size and prevalence is set for stage 1.
## In problem II, a sample size and prevalence are chosen in stage 1, while a fixed sample size and prevalence is set for stage 2.
## By taking the prior for problem II to be the posterior of problem I, and redefining the utility function for problem II,
## the SolveDP method should give (approximately) the same results.
## NOTE: wF=0.5,wS=0.5,wSC=0.5 are necessary. Symmetry then follows by the form of PHUtility.
Test.BackwardInduction=function(){
    ## Set the parameter values
    params=Parameters(
        lambda=0.5,
        alpha=0.025,etaS=0.3,etaSC=0.3,
        wF=0.5,wS=0.5,wSC=0.5, ## all equal to 0.5 important for symmetry
        N=10^6,r=10^4,
        cSetup=10^6,cBiomarker=10^6,cPerPatient=10^5,cScreening=10^3,
        muS=0.1,muF=0.1,
        sigma=1);

    ## Construct three point prior for (deltaS, deltaSC)
    dS=0.3;
    prior=cbind(c(dS,dS,dS),c(0,dS/2,dS),c(1/3,1/3,1/3));

    ## Define fixed actions and observations in stage 1 for problem I
    n1=100;
    lambda1=0.7;  
    zS1=-0.5;
    zSC1=0;       
    
    ## Construct decision problem I
    dpI=DecisionProblem(
        prior=prior,
        utility=PHUtility,
        n1s=c(n1),lambda1s=c(lambda1),epsilon1=10^-4,stepSize1=1,
        n2s=seq(100,1000,100),lambda2s=seq(0.1,0.9,0.1),epsilon2=10^-4,stepSize2=0.2);

    outI=OptStage2EU(dpI,params,n1,lambda1,zS1,zSC1);
    
    ## Construct decision problem II

    ## Get posterior corresponding to fixed observation for stage 1
    post=ComputePosterior(prior,zS1,zSC1,n1,lambda1,params$lambda,params$sigma);

    PHUtilityII=function(zS2Fixed,zSC2Fixed){
        function(n1,lambda1,zS1,zSC1,n2,lambda2,zS2,zSC2,deltaS,deltaSC,params){
             PHUtility(n1,lambda1,zS1,zSC1,n2,lambda2,zS2=zS2Fixed,zSC2=zSC2Fixed,deltaS,deltaSC,params)
        }
    }   
    
    dpII=DecisionProblem(
        prior=post,
        utility=PHUtilityII(zS1,zSC1),
        n1s=seq(100,1000,100),lambda1s=seq(0.1,0.9,0.1),epsilon1=10^-4,stepSize1=0.2,
        n2s=c(n1),lambda2s=c(lambda1),epsilon2=10^-4,stepSize2=1);

    outII=SolveDP(dpII,params);

    cat("Problem I, opt. utility =",outI$eu2,"\n")
    cat("Problem I, opt. stage 2 sample size =",outI$n2,"\n")
    cat("Problem I, opt. stage 2 prevalence =",outI$lambda2,"\n")
    cat("Problem II opt. utility =",outII$eu1,"\n")
    cat("Problem II opt. stage 1 sample size =",outII$n1,"\n")
    cat("Problem II opt. stage 1 prevalence =",outII$lambda1,"\n")
}
