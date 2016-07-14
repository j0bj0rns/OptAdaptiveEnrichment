## PlotFunctions
## Contains functions for plotting the results of the optimization.

PlotStage2EU=function(decisionMaker,design1,stage2Designs,prior,params,...){
    UseMethod("PlotStage2EU",design1);
}

## Generic function for Null first stage design.
PlotStage2EU.Null=function(deltaS,deltaSC,design1,design2,params){
    
}

## Plot the expected utility of stage 2, given a fixed stage 1 design.
PlotStage2EU=function(decisionMaker,design1,stage2Designs,prior,params,...){
    
    ## Construct a decision problem for the Public health view (integration approach)
    epsilon1=NA; ## Value not used in this function
    stepSize1=NA; ## Value not used in this function

    stage1Designs=list(design1);
    dp=if (decisionMaker=="PublicHealth"){
        DecisionProblemInt(prior=prior,utility=PHUtilityInt,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        epsilon1=epsilon1,stepSize1=stepSize1);
    } else if (decisionMaker=="Sponsor"){
        DecisionProblemInt(prior=prior,utility=SUtilityInt,
        stage1Designs=stage1Designs,stage2Designs=stage2Designs,
        epsilon1=epsilon1,stepSize1=stepSize1);
    }                      
    
    ## Compute Posterior
    post=ComputePosterior(design1,prior,params$sigma,...);   

    ## Compute the 2nd stage expected utilities
    design2=stage2Designs[[1]];
    if(class(design2)=="PartialEnrichment"){
        s2EUs=matrix(0,nrow=length(design2$nSC),ncol=length(design2$nS)); 
        
        for (j in 1:length(design2$nSC)) {
            for (i in 1:length(design2$nS)) {
                s2EUs[i,j]=Stage2EUInt(design1,PartialEnrichment(design2$nS[i],design2$nSC[j]),dp,params,post,...);
            }
        }

        ## Plot contour of stage 2 expected utility as function of nS and nSC                
        s2EUs=pmax(s2EUs,quantile(s2EUs,probs=0.8));
        filled.contour(design2$nS,design2$nSC,s2EUs/10^6,plot.title=title(main=paste("Expected Utility, partial enrichment,",decisionMaker)));
        mtext(expression(n[S]^2),side=1,line=2);
        mtext(expression(n[SC]^2),side=2,line=2);
    }
    
    if(class(stage2Designs[[1]])=="FullEnrichment"){
        s2EUs=vector(mode="numeric",length=length(design2$nS));
           
        for (i in 1:length(design2$nS)) {
            s2EUs[i]=Stage2EUInt(design1,FullEnrichment(design2$nS[i]),dp,params,post,...);
        }
     
        ## Plot stage 2 expected utility as function of nS
        s2EUs=pmax(s2EUs,quantile(s2EUs,probs=0.8));
        plot(design2$nS,s2EUs/10^6,
             title=title(main=paste("Expected Utility, full enrichment,",decisionMaker)),
             xlab=expression(n[S]^2),
             ylab="EU");
        
    }
}

## Example usage of PlotStage2EU
Example.PlotStage2EU=function(){
    ## Define decision maker
    decisionMaker="Sponsor";  
    
    ## Define stage 1 design and fixed results

    ## Null design
    design1=Null();

    ## FullEnrichment design
    ##design1=FullEnrichment(nS=100);
    ##zS1=0;
    
    ## PartialEnrichment design
    ##design1=PartialEnrichment(nS=100,nSC=100);
    ##zS1=0;
    ##zSC1=0;
    
    ## Define space of stage 2 designs
    n2SC=seq(50,500,20);
    n2S=seq(50,500,20);
##    stage2Designs=list(PartialEnrichment(n2S,n2SC));    
    stage2Designs=list(FullEnrichment(n2S));

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
    else if (class(design1)=="FullEnrichment"){
        PlotStage2EU(decisionMaker,design1,stage2Designs,prior,params,zS1);
    }
    else if (class(design1)=="PartialEnrichment"){
        PlotStage2EU(decisionMaker,design1,stage2Designs,prior,params,zS1,zSC1);
    }
}
