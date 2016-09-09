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

## Check stage 2 summation results against the stage 2 integration results.
## stage1Type = "Null" or "PartialEnrichment".
## stage2Type = "Null", "FullEnrichment" or "PartialEnrichment".
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
    
    epsilon2=10^-4;
    zGridPoints2=400;

    ## Define the fixed stage 1 design and stage 1 data
    design1=if (stage1Type=="Null") {
        design1=Null();
    } else if (stage1Type=="PartialEnrichment") {
        design1=PartialEnrichment(nS=100,nSC=100);
    }

    zS1=0.1;
    zSC1=-0.1;
    
    stage1Designs=list(design1);

    ## Define the set of possible two-stage designs
    n2S=seq(100,500,10);
    n2SC=200;
    
    stage2Designs=list();
    if (stage2Type=="Null") {
        stage2Designs=lapply(n2S,function(x) Null());
    } 
    if (stage2Type=="FullEnrichment") {
        stage2Designs=lapply(n2S,function(x) FullEnrichment(x));
    }    
    if (stage2Type=="PartialEnrichment") {
        stage2Designs=lapply(n2S,function(x) PartialEnrichment(x,n2SC));
    }

    ## Construct decision problems. Since only the stage 2 expectation is the only thing computed,
    ## many parameters may be set to NA. Only the utility function needs to be specified for these test objects.
    dpPHSum=DecisionProblem(prior=NA,utility=PHUtility,stage1Designs=NA,stage2Designs=NA,stage2Opt=NA,stage2Exp=NA,policyStepSize=NA,epsilon1=NA,zGridPoints1=NA);
    dpSSum=DecisionProblem(prior=NA,utility=SUtility,stage1Designs=NA,stage2Designs=NA,stage2Opt=NA,stage2Exp=NA,policyStepSize=NA,epsilon1=NA,zGridPoints1=NA);    
    dpPHInt=DecisionProblem(prior=NA,utility=PHUtilityInt,stage1Designs=NA,stage2Designs=NA,stage2Opt=NA,stage2Exp=NA,policyStepSize=NA,epsilon1=NA,zGridPoints1=NA);
    dpSInt=DecisionProblem(prior=NA,utility=SUtilityInt,stage1Designs=NA,stage2Designs=NA,stage2Opt=NA,stage2Exp=NA,policyStepSize=NA,epsilon1=NA,zGridPoints1=NA);
            
    ## Compute Posterior    
    post=if (stage1Type=="Null") {
        ComputePosterior(design1,prior,params$sigma);
    } else if (stage1Type=="PartialEnrichment") {
        ComputePosterior(design1,prior,params$sigma,zS1,zSC1);
    }

    ## Compute the 2nd stage expected utilities
    hw=0;
    Gz2=0;
    s2PHSumEUs=0;
    s2SSumEUs=0;
    s2PHIntEUs=0;
    s2SIntEUs=0;

    if (stage2Type=="FullEnrichment") {
        ## Create grid for zS2 (For FullEnrichment in second stage)
        hw=GridHalfWidth(epsilon2); ## Half width of grid square
        Gz2=seq(-hw,hw,length.out=zGridPoints2);               
    } else if (stage2Type=="PartialEnrichment") {
        ## Create grid for (zS2,zSC2) (For PartialEnrichment in second stage)
        hw=GridHalfWidth(epsilon2); ## Half width of grid square    
        Gz2=as.matrix(expand.grid(seq(-hw,hw,length.out=zGridPoints2),seq(-hw,hw,length.out=zGridPoints2)));
        colnames(Gz2)=NULL;
    }    

    if (stage1Type=="Null") {
        s2PHSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpPHSum,params,post,Gz2));
        s2SSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpSSum,params,post,Gz2));
        s2PHIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpPHInt,params,post));
        s2SIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpSInt,params,post));
    }  
    if (stage1Type=="PartialEnrichment") {
        s2PHSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpPHSum,params,post,Gz2,zS1,zSC1));
        s2SSumEUs=sapply(stage2Designs,function(d2) Stage2EUSum(design1,d2,dpSSum,params,post,Gz2,zS1,zSC1));
        s2PHIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpPHInt,params,post,zS1,zSC1));
        s2SIntEUs=sapply(stage2Designs,function(d2) Stage2EUInt(design1,d2,dpSInt,params,post,zS1,zSC1));
    }
    
    yMin=min(c(s2PHSumEUs,s2SSumEUs,s2PHIntEUs,s2SIntEUs));
    yMax=max(c(s2PHSumEUs,s2SSumEUs,s2PHIntEUs,s2SIntEUs));
    labels=c("Public health summation", "Public health integration","Sponsor summation", "Sponsor integration");
    col=c("black","black","red","red");
    
    plot(c(100, 500),c(yMax, yMin),type="n",xlab="n2S",ylab="Stage 2 EU")
    lines(n2S,s2PHSumEUs,lty=1,col="black");
    lines(n2S,s2PHIntEUs,lty=2,col="black");
    lines(n2S,s2SSumEUs,lty=1,col="red");
    lines(n2S,s2SIntEUs,lty=2,col="red");
    legend("topright",legend=labels,col=col,lty=c(1,2,1,2));
}

## Test OptStage2EU.
## decisionMaker = "Sponsor" or "PublicHealth".
## stage1Type = "Null" or "PartialEnrichment".
## stage2Opt = Optimization method = "Grid" or "Optim".
## stage2Exp = Expectation method = "Sum" or "Int".
Test.OptStage2EU=function(decisionMaker,stage1Type,stage2Opt,stage2Exp){
    
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

    ## Determine the fineness of the stage 2 summation grid
    epsilon2=10^-4;
    zGridPoints2=400;

    ## Define the fixed stage 1 design and stage 1 data
    design1=if (stage1Type=="Null") {
        Null();
    } 
    else if (stage1Type=="PartialEnrichment") {
        PartialEnrichment(nS=100,nSC=100);        
    }

    zS1=0.1;
    zSC1=1;
    
    stage1Designs=list(design1);

    ## Define the set of possible two-stage designs
    n2S=seq(100,500,20);
    n2SC=seq(100,500,20);
    
    stage2Designs=list(Null(),FullEnrichment(n2S),PartialEnrichment(n2S,n2SC));

    ## Construct decision problem
    
    utility=
        if (decisionMaker=="Sponsor" && stage2Exp=="Sum") SUtility
        else if (decisionMaker=="Sponsor" && stage2Exp=="Int") SUtilityInt
        else if (decisionMaker=="PublicHealth" && stage2Exp=="Sum") PHUtility
        else if (decisionMaker=="PublicHealth" && stage2Exp=="Int") PHUtilityInt
    
    dp=DecisionProblem(
        prior=prior,
        utility=utility,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        stage2Opt=stage2Opt,stage2Exp=stage2Exp,
        policyStepSize=NA,
        epsilon1=NA,zGridPoints1=NA,
        epsilon2=epsilon2,zGridPoints2=zGridPoints2);
    
    ## Optimise
    out=if (stage1Type=="Null") {
        OptStage2EU(design1,dp,params);          
    }
    else if (stage1Type=="PartialEnrichment") {
        OptStage2EU(design1,dp,params,zS1,zSC1);
    }

    ## Print results
    cat("Stage 2 optimization results:\n"); print(out);
}

## Use the symmetry over the two stages to check the backward induction algorithm.
## The numerical result of two decision problems are compared.
## In problem I, a sample size and prevalence are chosen in stage 2, while a fixed sample size and prevalence is set for stage 1.
## In problem II, a sample size and prevalence are chosen in stage 1, while a fixed sample size and prevalence is set for stage 2.
## By taking the prior for problem II to be the posterior of problem I, and redefining the utility function for problem II,
## the SolveDP method should give the same results.
## NOTE: wF=0.5,wS=0.5,wSC=0.5 are necessary. Symmetry then follows by the form of PHUtility.
Test.BackwardInduction=function(){
    ## Set the parameter values
    params=Parameters(
        lambda=0.5,
        alpha=0.025,etaS=0.3,etaSC=0.3,
        wF=0.5,wS=0.5,wSC=0.5, ## all equal to 0.5 important for symmetry
        N=10^6,r=10^4,
        cSetup=10^6,cBiomarker=10^6,cPerPatient=10^5,cScreening=10^5,
        muS=0.1,muF=0.1,
        sigma=sqrt(2));

    epsilon1=10^-4;
    zGridPoints1=21;
    epsilon2=10^-4;
    zGridPoints2=21;
    
    ## Construct the weak biomarker prior with dS=0.3
    dS=0.3;
    prior=cbind(c(0,dS,dS,dS),c(0,0,dS/2,dS),c(0.2,0.2,0.3,0.3));      

    ## Define fixed actions and observations in stage 1 for problem I
    n1S=10;
    n1SC=10;  
    zS1=0.1;
    zSC1=-0.2;       

    design1=PartialEnrichment(n1S,n1SC);
    stage1Designs=list(design1);
    stage2Designs=list(PartialEnrichment(seq(100,500,100),seq(100,500,100)));
    
    ## Construct and solve decision problem I
    dpI=DecisionProblem(
        prior=prior,utility=PHUtility,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        stage2Opt="Grid",stage2Exp="Sum",
        policyStepSize=1,
        epsilon1=epsilon1,zGridPoints1=zGridPoints1,epsilon2=epsilon2,zGridPoints2=zGridPoints2);    
    
    outI=OptStage2EU(design1,dpI,params,zS1,zSC1);
    
    ## Get posterior corresponding to fixed observation for stage 1
    post=ComputePosterior(design1,prior,params$sigma,zS1,zSC1);
        
    ## Construct and solve decision problem II
    PHUtilityII=function(zS2Fixed,zSC2Fixed){
        function(deltaS,deltaSC,design1,design2,params,zS1,zSC1,zS2,zSC2){
             PHUtility(deltaS,deltaSC,design1,design2,params,zS1,zSC1,zS2=zS2Fixed,zSC2=zSC2Fixed);
        }
    }   
    
    dpII=DecisionProblem(
        prior=post,utility=PHUtilityII(zS1,zSC1),
        stage1Designs=stage2Designs,stage2Designs=stage1Designs,
        stage2Opt="Grid",stage2Exp="Sum",
        policyStepSize=1,
        epsilon1=epsilon1,zGridPoints1=zGridPoints1,epsilon2=epsilon2,zGridPoints2=zGridPoints2); 

    outII=SolveDP(dpII,params);

    cat("Problem I, opt. utility =",outI$eu2,"\n")
    cat("Problem I, opt. stage 2 sample size (S)=",outI$design2$nS,"\n")
    cat("Problem I, opt. stage 2 sample size (SC)=",outI$design2$nSC,"\n")
    cat("Problem II opt. utility =",outII$eu1,"\n")
    cat("Problem II opt. stage 1 sample size (S)=",outII$design1$nS,"\n")
    cat("Problem II opt. stage 1 sample size (SC)=",outII$design1$nSC,"\n")
}

## Example of how to use SolveDP.
## In this example, the optimal policy is found for the public health view using the integration approach.
Example.SolveDP=function(){
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
    zGridPoints1=9;

    nS=seq(100,500,100);
    nSC=seq(100,500,100);
    
    stage1Designs=list(PartialEnrichment(nS,nSC));
    stage2Designs=list(PartialEnrichment(nS,nSC));
    
    ## Construct a decision problem
    dp=DecisionProblem(
        prior=prior,utility=PHUtilityInt,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        stage2Opt="Grid",stage2Exp="Int",
        policyStepSize=1,
        epsilon1=epsilon1,zGridPoints1=zGridPoints1);
   
    policy=SolveDP(dp,params);

    cat("Stage 1 optimal EU =",policy$eu1,"\n");
    cat("Stage 1 optimal design:\n"); print(policy$design1);
}
