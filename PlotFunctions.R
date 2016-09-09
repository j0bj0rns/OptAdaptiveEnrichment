## PlotFunctions
## Contains functions for plotting the results of the optimization.

## Simple plotting function adding a legend to 'image'
customImage=function(x,y,z,...){
    nLevels=20;

    old.par=par(mar=c(5,4,4,2+5)+0.1,mgp=c(2.5,1,0));
    
    colMap <- colorRampPalette(c("blue","red" ))(nLevels);
    image(x,y,z,col=colMap,...);

    #grconvertX(0.5, "device")y=grconvertY(1, "device"),
    legend(x="right",
           legend=format(c(max(z),(min(z)+max(z))/2,min(z)),digits=6),
           fill = colMap[c(nLevels, nLevels/2, 1)],
           xpd = NA,
           inset=c(-0.25,0));

    par(old.par);
}

## Plot the expected utility of stage 2, given a fixed stage 1 design.
PlotStage2EU=function(decisionMaker,design1,stage2Designs,prior,params,...){
    
    ## Construct the decision problem
    stage1Designs=list(design1);

    dp=DecisionProblem(
        prior=prior,
        utility=if (decisionMaker=="PublicHealth") PHUtilityInt else if (decisionMaker=="Sponsor") SUtilityInt,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        stage2Opt="Grid",stage2Exp="Int",
        policyStepSize=NA,
        epsilon1=NA,zGridPoints1=NA);   
    
    ## Compute Posterior
    post=ComputePosterior(design1,prior,params$sigma,...);   

    ## Compute the 2nd stage expected utilities
    design2=stage2Designs[[1]];
    if(class(design2)=="PartialEnrichment"){
        s2EU=matrix(0,nrow=length(design2$nS),ncol=length(design2$nSC)); 
        
        for (i in 1:length(design2$nS)) {
            for (j in 1:length(design2$nSC)) {
                s2EU[i,j]=Stage2EUInt(design1,PartialEnrichment(design2$nS[i],design2$nSC[j]),dp,params,post,...);
            }
        }

        ## Plot contour of stage 2 expected utility as function of nS and nSC
        customImage(design2$nS,design2$nSC,s2EU/10^6,
                    xlab=expression(n[S]^2),ylab=expression(n[SC]^2),
                    main=paste("Stage 2 EU, partial enrichment,",decisionMaker));
    }
    
    if(class(stage2Designs[[1]])=="FullEnrichment"){
        s2EU=vector(mode="numeric",length=length(design2$nS));
           
        for (i in 1:length(design2$nS)) {
            s2EU[i]=Stage2EUInt(design1,FullEnrichment(design2$nS[i]),dp,params,post,...);
        }
     
        ## Plot stage 2 expected utility as function of nS
        plot(design2$nS,s2EU/10^6,
             title=title(main=paste("Stage 2 EU, full enrichment,",decisionMaker)),
             xlab=expression(n[S]^2),ylab="EU");        
    }
}

## Example usage of PlotStage2EU
Example.PlotStage2EU=function(){
    ## Define decision maker
    decisionMaker="Sponsor";  
    
    ## Define stage 1 design and fixed results
    
    ## Null design
    design1=Null(); 
    
    ## PartialEnrichment design
    ##design1=PartialEnrichment(nS=100,nSC=100);
    ##zS1=0;
    ##zSC1=0;
    
    ## Define space of stage 2 designs
    n2S=seq(100,500,20);
    n2SC=seq(100,500,20);
    ##stage2Designs=list(FullEnrichment(n2S));
    stage2Designs=list(PartialEnrichment(n2S,n2SC));    

    ## Construct the weak biomarker prior with dS=0.3
    dS=0.3;
    prior=cbind(c(0,dS,dS,dS),c(0,0,dS/2,dS),c(0.2,0.2,0.3,0.3));   
    
    ## Set the parameter values
    params=Parameters(
        lambda=0.5,
        alpha=0.025,etaS=0.3,etaSC=0.3,
        wF=0.5,wS=0.5,wSC=0.5,
        N=10^6,r=10^4,
        cSetup=10^6,cBiomarker=10^6,cPerPatient=10^5,cScreening=10^5,
        muS=0.1,muF=0.1,
        sigma=sqrt(2));

    if (class(design1)=="Null"){
        PlotStage2EU(decisionMaker,design1,stage2Designs,prior,params);
    }
  
    if (class(design1)=="PartialEnrichment"){
        PlotStage2EU(decisionMaker,design1,stage2Designs,prior,params,zS1,zSC1);
    }
}

## Given a fixed partial enrichment design for stage 1,
## plot the optimal stage 2 EU and sample size choices as functions of the stage 1 outcome (zS, zSC).
PlotStage2Optimals=function(decisionMaker,design1,stage2Designs,prior,params){
    ## Construct the decision problem
    stage1Designs=list(design1);
    policyStepSize=0.5;
    epsilon1=10^-4;
    zGridPoints1=2; ## Smallest value allowed; it is policyStepSize which determines the step size in (zS,zSC)-space
    
    dp=DecisionProblem(
        prior=prior,
        utility=if (decisionMaker=="PublicHealth") PHUtilityInt else if (decisionMaker=="Sponsor") SUtilityInt,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        stage2Opt="Grid",stage2Exp="Int",
        policyStepSize=policyStepSize,
        epsilon1=epsilon1,zGridPoints1=zGridPoints1);    

    ## Solve the decision problem to obtain the optimal policy 
    policy=SolveDP(dp,params);
    
    ## Construct matrices for plotting using filled.contour
    zS=seq(policy$zSmin,policy$zSmax,policyStepSize);
    zSC=seq(policy$zSCmin,policy$zSCmax,policyStepSize);
    
    s2OptEU=matrix(0,nrow=length(zS),ncol=length(zSC));
    s2OptDesign=matrix(0,nrow=length(zS),ncol=length(zSC));
    s2OptnS=matrix(0,nrow=length(zS),ncol=length(zSC));
    s2OptnSC=matrix(0,nrow=length(zS),ncol=length(zSC));
        
    for (i in 1:length(zS)) {
        for (j in 1:length(zSC)) {
            s2OptEU[i,j]=policy$eu2(zS[i],zSC[j]);
            
            if (class(policy$design2(zS[i],zSC[j]))=="Null") {
                s2OptDesign[i,j]=0;
                s2OptnS[i,j]=0;
                s2OptnSC[i,j]=0;
            }

            if (class(policy$design2(zS[i],zSC[j]))=="FullEnrichment") {
                s2OptDesign[i,j]=1;
                s2OptnS[i,j]=policy$design2(zS[i],zSC[j])$nS;
                s2OptnSC[i,j]=0;
            }

            if (class(policy$design2(zS[i],zSC[j]))=="PartialEnrichment") {
                s2OptDesign[i,j]=2;
                s2OptnS[i,j]=policy$design2(zS[i],zSC[j])$nS;
                s2OptnSC[i,j]=policy$design2(zS[i],zSC[j])$nSC;
            }            
        }
    }

    ## Plot the optimal stage 2 utility as a function of (zS,zSC)
    x11();
    customImage(zS,zSC,s2OptEU/10^6,
                xlab=expression(z[S]^1),ylab=expression(z[SC]^1),
                main=paste("Stage 2 optimal EU,",decisionMaker));         
    
    ## Plot the optimal stage 2 design type as a function of (zS,zSC)
    x11();
    customImage(zS,zSC,s2OptDesign,
                xlab=expression(z[S]^1),ylab=expression(z[SC]^1),
                main=paste("Stage 2 optimal design,",decisionMaker)); 

    
    ## Plot the optimal stage 2 nS as a function of (zS,zSC)
    x11();
    customImage(zS,zSC,s2OptnS,
                xlab=expression(z[S]^1),ylab=expression(z[SC]^1),
                main=paste("Stage 2 optimal nS,",decisionMaker));
    
    ## Plot the optimal stage 2 nSC as a function of (zS,zSC)
    x11();
    customImage(zS,zSC,s2OptnSC,
                xlab=expression(z[S]^1),ylab=expression(z[SC]^1),
                main=paste("Stage 2 optimal nSC,",decisionMaker));
}

## Example usage of PlotStage2Optimals
Example.PlotStage2Optimals=function(){
    ## Define decision maker
    decisionMaker="Sponsor";  
    
    ## Define stage 1 design   
    design1=PartialEnrichment(nS=100,nSC=100); 
        
    ## Define space of stage 2 designs
    n2S=seq(100,500,100);
    n2SC=seq(100,500,100);
    stage2Designs=list(Null(),FullEnrichment(n2S),PartialEnrichment(n2S,n2SC));

    ## Construct the weak biomarker prior with dS=0.3
    dS=0.3;
    prior=cbind(c(0,dS,dS,dS),c(0,0,dS/2,dS),c(0.2,0.2,0.3,0.3));   
    
    ## Set the parameter values
    params=Parameters(
        lambda=0.5,
        alpha=0.025,etaS=0.3,etaSC=0.3,
        wF=0.5,wS=0.5,wSC=0.5,
        N=10^6,r=10^4,
        cSetup=10^6,cBiomarker=10^6,cPerPatient=10^5,cScreening=10^5,
        muS=0.1,muF=0.1,
        sigma=sqrt(2));
  
    PlotStage2Optimals(decisionMaker,design1,stage2Designs,prior,params);
}


## PlotStage1EU: Plot the expected stage 1 utility (given optimal continuation in stage 2) as a function of the sample size parameters for a stage 1 partial enrichment design.
PlotStage1EU=function(decisionMaker,design1,stage2Designs,prior,params){
    
    ## Construct the decision problem
    stage1Designs=list(design1);
    epsilon1=10^-2;
    zGridPoints1=5;
    
    dp=DecisionProblem(
        prior=prior,
        utility=if (decisionMaker=="PublicHealth") PHUtilityInt else if (decisionMaker=="Sponsor") SUtilityInt,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        stage2Opt="Grid",stage2Exp="Int",
        policyStepSize=NA,
        epsilon1=epsilon1,zGridPoints1=zGridPoints1);

    s1EU=matrix(0,nrow=length(design1$nS),ncol=length(design1$nSC)); 
    
    for (i in 1:length(design1$nS)) {
        for (j in 1:length(design1$nSC)) {
            s1EU[i,j]=Stage1EU(PartialEnrichment(design1$nS[i],design1$nSC[j]),dp,params);
        }
    }

    ## Plot stage 1 expected utility as function of nS and nSC
    customImage(design1$nS,design1$nSC,s1EU/10^6,
                xlab=expression(n[S]^2),ylab=expression(n[SC]^2),
                main=paste("Expected Utility, partial enrichment,",decisionMaker));    
}

## Example usage of PlotStage1EU
Example.PlotStage1EU=function(){
    ## Define decision maker
    decisionMaker="PublicHealth";  
    
    ## Define space of stage 1 designs
    design1=PartialEnrichment(nS=seq(100,500,100),nSC=seq(100,500,100)); 
        
    ## Define space of stage 2 designs
    n2S=seq(100,500,100);
    n2SC=seq(100,500,100);
    stage2Designs=list(Null(),FullEnrichment(n2S),PartialEnrichment(n2S,n2SC));

    ## Construct the weak biomarker prior with dS=0.3
    dS=0.3;
    prior=cbind(c(0,dS,dS,dS),c(0,0,dS/2,dS),c(0.2,0.2,0.3,0.3));   
    
    ## Set the parameter values
    params=Parameters(
        lambda=0.5,
        alpha=0.025,etaS=0.3,etaSC=0.3,
        wF=0.5,wS=0.5,wSC=0.5,
        N=10^6,r=10^4,
        cSetup=10^6,cBiomarker=10^6,cPerPatient=10^5,cScreening=10^5,
        muS=0.1,muF=0.1,
        sigma=sqrt(2));
  
    PlotStage1EU(decisionMaker,design1,stage2Designs,prior,params);
}
