## Optimization
## Contains basic data structures and functions for optimizing over both stages.

###################
# Data structures #
###################

##
# Parameters. Construct a named list with the problem parameters.
##
Parameters=function(lambda,alpha,etaS,etaSC,wF,wS,wSC,N,r,cSetup,cBiomarker,cPerPatient,cScreening,muS,muF,sigma) {
    return(list(lambda=lambda,
                alpha=alpha,etaS=etaS,etaSC=etaSC,
                wF=wF,wS=wS,wSC=wSC,
                N=N,r=r,
                cSetup=cSetup,cBiomarker=cBiomarker,cPerPatient=cPerPatient,cScreening=cScreening,
                muS=muS,muF=muF,
                sigma=sigma,
                z_etaS=qnorm(1-etaS),z_etaSC=qnorm(1-etaSC)));
}

##
# DecisionProblem. Construct a named list with a specification of the decision problem.
# prior = Discrete prior for the effects sizes, represented using a matrix with three columns.
# utility = Utility function to be used (i.e., sponsor or public health view).
# stage1Designs = list containing the possible stage 1 designs.
# stage2Designs = list containing the possible stage 2 designs.

# stage2Opt = The optimization method for the second stage, equal to either "Grid" or "Optim".
# stage2Exp = The method used to compute the second stage expectation, equal to either "Sum" or "Int".
# policyStepSize = The step size used in the discrete approximation to the stage 2 policy function.
# epsilon1,epsilon2 = Probabilities determining the size of the region of probable Z-values.
# zGridPoints1,zGridPoints2 = Number of points in each dimension in the Z-grid when using the summation approach.
# epsilon2 and stepSize2 need not be specified if stage2Exp = "Int".
#
# In order to support the specification of different forms of the utility functions for different design types,
# the utility functions must be implemented as generic functions (via R's S3 object system).
# The generic function should take the arguments (deltaS,deltaSC,design1,design2,params,...).
# Method dispatch must be implemented for the two design arguments, design1 and design2.
# If stage2Exp="Sum", then the dots should be replaced by the stage 1 z-values followed by the stage 2 z-values,
# as appropriate for the designs used.
# If stage2Exp="Int", then the dots should be replaced by the stage 1 z-values only,
# since integration is done over the stage 2 z-values.
#
# Note that the stage 1 expected utility is always computed using the summation approach, optimizing over a grid.
##
DecisionProblem=function(prior,utility,stage1Designs,stage2Designs,stage2Opt,stage2Exp,policyStepSize,epsilon1,zGridPoints1,epsilon2=NA,zGridPoints2=NA) {
    return(list(prior=prior,utility=utility,
                stage1Designs=stage1Designs,stage2Designs=stage2Designs,
                stage2Opt=stage2Opt,stage2Exp=stage2Exp,
                policyStepSize=policyStepSize,
                epsilon1=epsilon1,zGridPoints1=zGridPoints1,
                epsilon2=epsilon2,zGridPoints2=zGridPoints2));
}

##
# Constructors for objects representing different stage 1 and stage 2 designs.
# nS = n * trialLambda (e.g., for stage 1 n1S = n1 * lambda1).
# nSC = n * (1 - trialLambda) (e.g., for stage 2 n2SC = n2 * (1 - lambda2).

# When used in DecisionProblem objects, nS and nSC are vectors containing the possible values for the subgroup sample sizes.
# If the summation approach is used in stage 2, optimisation is done over all the values.
# If the integration approach is used in stage 2, optimisation is done over the range defined by the min and max values.

# These design objects are used as inputs when specifying the choice space for a decision problem.
# They are also used as outputs, in which case nS and nSC consists of single points,
# namely the optimal values found during optimization.
##

## S3 class representing a single choice of no observations.
## Using a Null design in stage 1 corresponds to a reduction of the two-stage problem
## to a single stage optimization problem.
## Using a Null design in stage 2 corresponds to early stopping (for futility).
Null=function() {
    return(structure(list(),class="Null"));
}

## S3 class representing a full enrichment design (lambda2 = 1).
FullEnrichment=function(nS) {
    return(structure(list(nS=nS),class="FullEnrichment")); 
}

## S3 class representing a partial enrichment design.
PartialEnrichment=function(nS,nSC) {
    return(structure(list(nS=nS,nSC=nSC),class="PartialEnrichment"));   
}

##
# Policy. Represents an action policy for both stages.
# eu1 = Optimal stage 1 expected utility, a single number.
# design1 = Optimal stage 1 design.
# eu2 = Optimal stage 2 expected utility, either a function mapping the stage 1 outcome to the corresponding optimal stage 2 EU, or just a single number in case the stage 1 design is Null.
# design2 = Optimal stage 2 design, a function mapping the stage 1 outcome to the corresponding optimal stage 2 design, or just a single design in case the stage 1 design is Null.
# zSmin = Min value of zS1 for when used as an argument to the eu2 and design2 policy functions.
# zSmax = Max value of zS1 for when used as an argument to the eu2 and design2 policy functions.
# zSCmin = Min value of zSC1 for when used as an argument to the eu2 and design2 policy functions.
# zSCmax = Max value of zSC1 for when used as an argument to the eu2 and design2 policy functions.
##
Policy=function(eu1,design1,eu2,design2,zSmin,zSmax,zSCmin,zSCmax) {
    return(list(eu1=eu1,design1=design1,
                eu2=eu2,design2=design2,
                zSmin=zSmin,zSmax=zSmax,
                zSCmin=zSCmin,zSCmax=zSCmax));
}


###########################
# Some basic functions    #
###########################

##
# GridHalfWidth. Compute the half width of the square grid for (zS, zSC) that will give a probability
# of 1 - epsilon for an outcome to end up in the square.
# If w denotes the half width, then the coordinates of the top-right and bottom-left corner of the square
# are (w, w) and (-w,-w), respectively, given that the center of the square is at the origin.
# The formula which the algorithm is based on is easily shown using that zS and zSC are independent and
# standard normal under the assumption that (deltaS, deltaSC) = (0, 0).
##
GridHalfWidth=function(epsilon){   
    return(qnorm((1+sqrt(1-epsilon))/2));
}

##
# SampleDensity. The bivariate normal density for the Z-statistics in stage 1 and stage 2, given deltaS and deltaSC.
##
SampleDensity=function(deltaS,deltaSC,zS,zSC,n,trialLambda,sigma){
    return(dmvnorm(cbind(zS,zSC),
                   mean=c(deltaS*sqrt(trialLambda*n/sigma^2),deltaSC*sqrt((1-trialLambda)*n/sigma^2)),
                   sigma=diag(2)));
}

##
# SampleDensity1. The univariate normal density for the Z-statistic in stage 1 and stage 2, given deltaS.
##
SampleDensity1=function(deltaS,zS,n,sigma){
    return(dnorm(zS,mean=deltaS*sqrt(n/sigma^2),sd=1));
}

## Normalise (i.e., scale) a numeric vector so that the sum of its entries becomes 1.
normalise=function(x){
    return(x/sum(x));
}

##
# ComputePosterior. Computes the posterior for the trial effect sizes deltaS and deltaSC.
# The prior for the pair (deltaS, deltaSC) is represented using a matrix.
# Each row in the submatrix consisting of the first two columns represents a possible point.
# The third column in the same row assigns the prior probability for the point.
# The way in which the posterior must be computed depends on the stage 1 design, 
# and therefore ComputePosterior is implemented as a generic function.
##
ComputePosterior=function(design1,prior,sigma,...){
    UseMethod("ComputePosterior",design1)
}

ComputePosterior.Null=function(design1,prior,sigma){
    ## Just return the prior as the posterior, since nothing has been observed with a Null design
    return(prior);
}

ComputePosterior.PartialEnrichment=function(design1,prior,sigma,zS1,zSC1){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
    
    ## Compute the sampling probabilities
    priorProbs=prior[,3];
    sampleProbs=mapply(SampleDensity,prior[,1],prior[,2],
        MoreArgs=list(zS=zS1,zSC=zSC1,n=n1,trialLambda=lambda1,sigma=sigma));
    
    ## Return the updated matrix representing the posterior distribution
    return(cbind(prior[,1:2,drop=FALSE],normalise(priorProbs*sampleProbs)));
}

#######################################################
# Functions for optimization using backward induction #
#######################################################

##
# SolveDP. Solve a decision problem and return an optimal policy.
# dp = Decision problem object.
# params = Parameter object.
##
SolveDP=function(dp,params){
    optStage1EU=-Inf;
    optStage1Design=NA;

    for (d in dp$stage1Designs) {
        if (class(d)=="Null") {
            design1=Null();
            out=Stage1EU(design1,dp,params);
            
            if(out>optStage1EU){
                optStage1EU=out;
                optStage1Design=design1;
            }
        }        
        
        if (class(d)=="PartialEnrichment") {
            for (nS in d$nS) {
                for (nSC in d$nSC) {
                    design1=PartialEnrichment(nS,nSC);
                    out=Stage1EU(design1,dp,params);

                    if(out>optStage1EU){                      
                        optStage1EU=out;
                        optStage1Design=design1;
                    }                    
                }
            }
        }        
    }

    return(Stage2Policy(optStage1Design,optStage1EU,dp,params));
}

##
# Stage1EU. Compute the expected utility for a stage 1 design, given that stage 2 is performed optimally.
##
Stage1EU=function(design1,...) {
    UseMethod("Stage1EU",design1);
}

Stage1EU.Null=function(design1,dp,params){
    out=OptStage2EU(design1,dp,params);
    return(out$eu2);
}

Stage1EU.PartialEnrichment=function(design1,dp,params){
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
    
    ## Create grid for (zS1,zSC1)
    hw=GridHalfWidth(dp$epsilon1); ## Half width of grid square
    Gz1=as.matrix(expand.grid(seq(-hw,hw,length.out=dp$zGridPoints1),seq(-hw,hw,length.out=dp$zGridPoints1)));
    colnames(Gz1)=NULL;

    ## Compute discretized version of the sample density over a grid centered at (0,0)
    sd=cbind(Gz1[,1],Gz1[,2],normalise(SampleDensity(deltaS=0,deltaSC=0,zS=Gz1[,1],zSC=Gz1[,2],n=n1,trialLambda=lambda1,params$sigma)));
    
    Stage1CondEU=function(deltaS,deltaSC) {
        ## Shift the grid points and compute the utility for each grid point
        zS1=sd[,1]+deltaS*sqrt(design1$nS/params$sigma^2);
        zSC1=sd[,2]+deltaSC*sqrt(design1$nSC/params$sigma^2);          
        condEU1=mapply(function(x,y) OptStage2EU(design1,dp,params,x,y)$eu2,zS1,zSC1);        
        
        ## Return expected utility
        return(sum(sd[,3]*condEU1));
    }       
    
    ## Compute and return expected utility    
    return(sum(dp$prior[,3]*mapply(Stage1CondEU,dp$prior[,1],dp$prior[,2])));    
}
    
##
# OptStage2EU. Optimise the stage 2 expected utility.
##
OptStage2EU=function(design1,dp,params,...) {
    ## Compute posterior distribution given stage 1 data
    post=ComputePosterior(design1,dp$prior,params$sigma,...);

    optimals=function(eu2,stage2Design){
        return(list(eu2=eu2,stage2Design=stage2Design));
    }
    
    opts=optimals(-Inf,NA);
    out=0;
    
    for (d in dp$stage2Designs) {
        
        if (class(d)=="Null") {
            design2=Null();
            out=if (dp$stage2Exp=="Sum") {
                Stage2EUSum(design1,design2,dp,params,post,...);
            } else if (dp$stage2Exp=="Int"){
                Stage2EUInt(design1,design2,dp,params,post,...);
            }
            
            if(out>opts$eu2){
                opts=optimals(out,design2);
            }
        }
        
        if (class(d)=="FullEnrichment") {
            ## Create wrapper function to optimize
            f=if (dp$stage2Exp=="Sum") {
                ## Create grid for zS2
                hw=GridHalfWidth(dp$epsilon2); ## Half width of grid square
                Gz2=seq(-hw,hw,length.out=dp$zGridPoints2);
                
                function(x){
                    design2=FullEnrichment(x[1]);
                    Stage2EUSum(design1,design2,dp,params,post,Gz2,...)
                }
            } else if (dp$stage2Exp=="Int") {
                function(x){
                    design2=FullEnrichment(x[1]);
                    Stage2EUInt(design1,design2,dp,params,post,...)
                }
            }
            
            if (dp$stage2Opt=="Grid") {
                ## Compute stage 2 expected utility for each stage 2 design         
                for (nS in d$nS) {
                    out=f(nS);
                    if(out>opts$eu2){
                        opts=optimals(out,FullEnrichment(nS));
                    }  
                }    
            }

            if (dp$stage2Opt=="Optim") {
                ## Set start point for local optima search
                start=(min(d$nS)+max(d$nS))/2;
                
                ## Find optimal value of n2S and return the result
                opt=optim(start,fn=f,method="L-BFGS-B",
                    lower=min(d$nS),upper=max(d$nS),control=list(fnscale=-1));                 
                
                if(opt$value>opts$eu2){                
                    opts=optimals(opt$value,FullEnrichment(opt$par[1]));
                }
            }
        }
        
        if (class(d)=="PartialEnrichment") {
            ## Create wrapper function to optimize           
            f=if (dp$stage2Exp=="Sum") {
                ## Create grid for (zS2,zSC2)
                hw=GridHalfWidth(dp$epsilon2); ## Half width of grid square
                Gz2=as.matrix(expand.grid(seq(-hw,hw,length.out=dp$zGridPoints2),seq(-hw,hw,length.out=dp$zGridPoints2)));
                colnames(Gz2)=NULL;
                    
                function(x){
                    design2=PartialEnrichment(x[1],x[2]);
                    Stage2EUSum(design1,design2,dp,params,post,Gz2,...)
                }
            } else if (dp$stage2Exp=="Int") {
                function(x){
                    design2=PartialEnrichment(x[1],x[2]);
                    Stage2EUInt(design1,design2,dp,params,post,...)
                }
            }
            
            if (dp$stage2Opt=="Grid") {
                ## Compute stage 2 expected utility for each stage 2 design         
                for (nS in d$nS) {
                    for (nSC in d$nSC) {
                        out=f(c(nS,nSC));
                        if(out>opts$eu2){
                            opts=optimals(out,PartialEnrichment(nS,nSC));
                        }  
                    }
                }                
            }

            if (dp$stage2Opt=="Optim") {           
                ## Set start point for local optima search
                start=c((min(d$nS)+max(d$nS))/2,(min(d$nSC)+max(d$nSC))/2);
                
                ## Find optimal combination of n2S and n2SC and return the result
                opt=optim(start,fn=f,method="L-BFGS-B",
                    lower=c(min(d$nS),min(d$nSC)),upper=c(max(d$nS),max(d$nSC)),control=list(fnscale=-1));                 
                
                if(opt$value>opts$eu2){                
                    opts=optimals(opt$value,PartialEnrichment(opt$par[1],opt$par[2]));
                }
            }       
        }    
    }

    return(list(eu2=opts$eu2,design2=opts$stage2Design));
}

##
# MatFun. Convert a 2-dimensional matrix holding values at grid points into a function.
# M = the matrix.
# xMin = x-position corresponding to first row of M.
# yMin = y-position corresponding to first column of M.
# stepSize = step size between grid points (in both dimensions).
##
MatFun=function(M,xMin,yMin,stepSize){
    return(function(x,y){
        row=round((x-xMin)/stepSize)+1;
        col=round((y-yMin)/stepSize)+1;
        return(M[cbind(row,col)]);
    });
}

##
# MatDesignFun. Construct a function returning the optimal stage 2 design, which depends on the stage 1 data.
# designMat = Matrix holding the design type.
# nSMat = Optimal sample size for subgroup S, if applicable.
# nSCMat = Optimal sample size for subgroup SC, if applicable.
# xMin = x-position corresponding to first row of the matrices.
# yMin = y-position corresponding to first column of the matrices.
# stepSize = step size between grid points (in both dimensions).
##
MatDesignFun=function(designMat,nSMat,nSCMat,xMin,yMin,stepSize){
    return(function(x,y){
        row=round((x-xMin)/stepSize)+1;
        col=round((y-yMin)/stepSize)+1;
        
        d=designMat[cbind(row,col)];
        if (d=="Null") { Null(); }
        else if (d=="FullEnrichment") { FullEnrichment(nSMat[row,col]); }
        else if (d=="PartialEnrichment") { PartialEnrichment(nSMat[row,col],nSCMat[row,col]); }
    });
}

##
# Stage2Policy. Construct the policy of using design1 in stage 1 and then continue optimally in stage 2.
##
Stage2Policy=function(design1,...) {
    UseMethod("Stage2Policy",design1);
}

Stage2Policy.Null=function(design1,eu1,dp,params){
    out=OptStage2EU(design1,dp,params);
    return(Policy(eu1=eu1,design1=design1,eu2=out$eu2,design2=out$design2,
                  zSmin=NA,zSmax=NA,zSCmin=NA,zSCmax=NA));
}

Stage2Policy.PartialEnrichment=function(design1,eu1,dp,params){
    ## Construct a grid containing all sufficiently probable pairs (zS1, zSC1),
    ## as (deltaS, deltaSC) varies over the prior points.   
    hw=GridHalfWidth(dp$epsilon1);  

    priorMeansZS=dp$prior[,1]*sqrt(design1$nS/params$sigma^2);
    priorMeansZSC=dp$prior[,2]*sqrt(design1$nSC/params$sigma^2);
    
    minStepZS=floor((min(priorMeansZS)-hw)/dp$policyStepSize);
    maxStepZS=ceiling((max(priorMeansZS)+hw)/dp$policyStepSize);
    minStepZSC=floor((min(priorMeansZSC)-hw)/dp$policyStepSize);
    maxStepZSC=ceiling((max(priorMeansZSC)+hw)/dp$policyStepSize);
        
    ## Initialize matrices that will hold the results from the stage 2 optimization
    matRows=maxStepZS-minStepZS+1;
    matCols=maxStepZSC-minStepZSC+1;

    eu2Mat=matrix(nrow=matRows,ncol=matCols);
    designMat=matrix(nrow=matRows,ncol=matCols);
    nSMat=matrix(nrow=matRows,ncol=matCols);
    nSCMat=matrix(nrow=matRows,ncol=matCols);  

    ## Construct and return the optimal policy for stage 2
    for(i in 1:matRows) {
        zS1=(minStepZS+i-1)*dp$policyStepSize;
        
        for(j in 1:matCols) {
            zSC1=(minStepZSC+j-1)*dp$policyStepSize;

            out=OptStage2EU(design1,dp,params,zS1,zSC1);
            
            designType=class(out$design2);
            eu2Mat[i,j]=out$eu2;                                        
            designMat[i,j]=designType;
            
            if(designType=="FullEnrichment") {
                nSMat[i,j]=out$design2$nS;
            }
            
            if(designType=="PartialEnrichment") {
                nSMat[i,j]=out$design2$nS;
                nSCMat[i,j]=out$design2$nSC;
            }
        }                           
    }
    
    zSmin=minStepZS*dp$policyStepSize;
    zSmax=maxStepZS*dp$policyStepSize;
    zSCmin=minStepZSC*dp$policyStepSize;
    zSCmax=maxStepZSC*dp$policyStepSize;
    
    return(Policy(eu1=eu1,design1=design1,
                  eu2=MatFun(eu2Mat,zSmin,zSCmin,dp$policyStepSize),
                  design2=MatDesignFun(designMat,nSMat,nSCMat,zSmin,zSCmin,dp$policyStepSize),
                  zSmin=zSmin,zSmax=zSmax,zSCmin=zSCmin,zSCmax=zSCmax));
}

##
# Stage2EUSum. Compute the stage 2 expected utility given the stage 1 data using summation.
# sd = Discretised sample distribution for zS2 or (zS2, zSC2) with center 0 or (0,0).
# Each row in the submatrix consisting of the first (or first two) columns represents a possible point.
# The last column in the same row assigns the probability for the point.
##
Stage2EUSum=function(design1,design2,dp,params,post,...){
    UseMethod("Stage2EUSum",design2);
}

Stage2EUSum.Null=function(design1,design2,dp,params,post,...){
    ## Compute and return expected utility
    stage2CondEUs=mapply(Stage2CondEUSum,post[,1],post[,2],
        MoreArgs=list(design1,design2,dp,params,...));
    
    return(sum(post[,3]*stage2CondEUs));
}

Stage2EUSum.FullEnrichment=function(design1,design2,dp,params,post,Gz2,...){
    ## Compute discretized version of the sample density over a grid centered at 0
    sd=cbind(Gz2,normalise(SampleDensity1(deltaS=0,zS=Gz2,n=design2$nS,params$sigma))); 
    
    ## Compute and return expected utility
    stage2CondEUs=mapply(Stage2CondEUSum,post[,1],post[,2],
        MoreArgs=list(design1,design2,dp,params,sd,...));
    
    return(sum(post[,3]*stage2CondEUs));
}

Stage2EUSum.PartialEnrichment=function(design1,design2,dp,params,post,Gz2,...){
    n2=design2$nS+design2$nSC;
    lambda2=design2$nS/n2;

    ## Compute discretized version of the sample density over a grid centered at (0,0)
    sd=cbind(Gz2[,1],Gz2[,2],normalise(SampleDensity(deltaS=0,deltaSC=0,zS=Gz2[,1],zSC=Gz2[,2],n=n2,trialLambda=lambda2,params$sigma))); 
    
    ## Compute and return expected utility
    stage2CondEUs=mapply(Stage2CondEUSum,post[,1],post[,2],
        MoreArgs=list(design1,design2,dp,params,sd,...));

    return(sum(post[,3]*stage2CondEUs));
}

##
# Stage2CondEUSum. Compute the stage 2 expected utility given deltaS and deltaSC.
##
Stage2CondEUSum=function(deltaS,deltaSC,design1,design2,dp,params,...){
    UseMethod("Stage2CondEUSum",design2);
}

Stage2CondEUSum.Null=function(deltaS,deltaSC,design1,design2,dp,params,...){
    ## Return expected utility
    return(dp$utility(deltaS,deltaSC,design1,design2,params,...));
}

Stage2CondEUSum.FullEnrichment=function(deltaS,deltaSC,design1,design2,dp,params,sd,...){
    ## Shift the grid points and compute utility for each grid point
    zS2=sd[,1]+deltaS*sqrt(design2$nS/params$sigma^2);
    condEU2=dp$utility(deltaS,deltaSC,design1,design2,params,...,zS2);
    
    ## Return expected utility
    return(sum(sd[,2]*condEU2));
}

Stage2CondEUSum.PartialEnrichment=function(deltaS,deltaSC,design1,design2,dp,params,sd,...){    
    ## Shift the grid points and compute the utility for each grid point
    zS2=sd[,1]+deltaS*sqrt(design2$nS/params$sigma^2);
    zSC2=sd[,2]+deltaSC*sqrt(design2$nSC/params$sigma^2);    
    condEU2=dp$utility(deltaS,deltaSC,design1,design2,params,...,zS2,zSC2);    
    
    ## Return expected utility
    return(sum(sd[,3]*condEU2));
}

##
# Stage2EUInt. Compute the stage 2 expected utility given the stage 1 data using integration.
##
Stage2EUInt=function(design1,design2,dp,params,post,...){
    u=mapply(dp$utility,post[,1],post[,2],MoreArgs=list(design1,design2,params,...));    
    return(sum(u*post[,3]));
}
