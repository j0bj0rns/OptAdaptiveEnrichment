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
                z_etaS=qnorm(1-etaS),z_etaSC=qnorm(1-etaSC),z_ahalf=qnorm(1-alpha/2)));
}

##
# DecisionProblem. Construct a named list with a specification of the decision problem.
# There are two types: type="Int" uses integration in stage 2, whereas type="Sum" uses summation.
#
# In order to support the specification of different forms of the utility functions for different design types,
# the utility functions must be implemented as generic functions (via R's S3 object system).
# The generic function should take the arguments (deltaS,deltaSC,design1,design2,params,...).
# Method dispatch must be implemented for the two design arguments, design1 and design2.
# If type="Sum", then the dots should be replaced by the stage 1 z-values followed by the stage 2 z-values,
# as appropriate for the designs used.
# If type="Int", then the dots should be replaced by the stage 1 z-values only,
# since integration is done over the stage 2 z-values.
#
# Note that the stage 1 expected utility is always computed using the summation approach.
##

DecisionProblemSum=function(prior,utility,stage1Designs,stage2Designs,epsilon1,stepSize1,epsilon2,stepSize2) {
    return(list(prior=prior,utility=utility,
                stage1Designs=stage1Designs,stage2Designs=stage2Designs,
                epsilon1=epsilon1,stepSize1=stepSize1,
                epsilon2=epsilon2,stepSize2=stepSize2,
                type="Sum"));
}

DecisionProblemInt=function(prior,utility,stage1Designs,stage2Designs,epsilon1,stepSize1) {
    return(list(prior=prior,utility=utility,
                stage1Designs=stage1Designs,stage2Designs=stage2Designs,
                epsilon1=epsilon1,stepSize1=stepSize1,
                type="Int"));
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

## S3 class representing a single choice of no observations (n1 = 0).
Null=function() {
    return(structure(list(),class="Null"));
}

## S3 class representing a full enrichment design (lambda1 = 1 or lambda2 = 1).
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
# eu2 = Optimal stage 2 expected utility, a function mapping the stage 1 outcome to the corresponding optimal stage 2 EU.
# design2 = Optimal stage 2 design, a function mapping the stage 1 outcome to the corresponding optimal stage 2 design.
# zSmin = Min value of zS1 for the grid approximation.
# zSmax = Max value of zS1 for the grid approximation.
# stepSize1 = Step size used in grid approximation.
##
Policy=function(eu1,design1,eu2,design2,zSmin,zSmax,stepSize1) {
    return(list(eu1=eu1,design1=design1,
                eu2=eu2,design2=design2,
                zSmin=zSmin,zSmax=zSmax,stepSize1=stepSize1));
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
                   mean=c(deltaS*sqrt(trialLambda*n/sigma^2),deltaSC*sqrt((1 - trialLambda)*n/sigma^2)),
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

ComputePosterior.FullEnrichment=function(design1,prior,sigma,zS1){
    n1=design1$nS;
    
    ## Compute the sampling probabilities
    priorProbs=prior[,3];
    sampleProbs=mapply(SampleDensity1,prior[,1],MoreArgs=list(zS=zS1,n=n1,sigma=sigma));   
    
    ## Return the updated matrix representing the posterior distribution
    return(cbind(prior[,1:2,drop=FALSE],normalise(priorProbs*sampleProbs)));
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
    optimals=function(eu1,eu2,stage1Design,stage2Design,zSmin,zSmax){
        return(list(eu1=eu1,eu2=eu2,stage1Design=stage1Design,stage2Design=stage2Design,zSmin=zSmin,zSmax=zSmax))
    }
    
    opts=optimals(-Inf,NA,NA,NA,NA,NA);
   
    for (d in dp$stage1Designs) {
        if (class(d)=="Null") {
            design1=Null();
            out=Stage1EU(design1,dp,params);
            
            if(out$eu1>opts$eu1){
                opts=optimals(out$eu1,out$eu2,design1,out$design2,out$zSmin,out$zSmax);
            }
        }

        if (class(d)=="FullEnrichment") {
            for (nS in d$nS) {
                design1=FullEnrichment(nS);
                out=Stage1EU(design1,dp,params);

                if(out$eu1>opts$eu1){
                    opts=optimals(out$eu1,out$eu2,design1,out$design2,out$zSmin,out$zSmax);
                }
            }
        }
        
        if (class(d)=="PartialEnrichment") {
            for (nS in d$nS) {
                for (nSC in d$nSC) {
                    design1=PartialEnrichment(nS,nSC);
                    out=Stage1EU(design1,dp,params);

                    if(out$eu1>opts$eu1){
                        opts=optimals(out$eu1,out$eu2,design1,out$design2,out$zSmin,out$zSmax);
                    }                    
                }
            }
        }        
    }

    return(Policy(eu1=opts$eu1,design1=opts$stage1Design,eu2=opts$eu2,design2=opts$stage2Design,
                  zSmin=opts$zSmin,zSmax=opts$zSmax,stepSize1=dp$stepSize1));
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
# VecFun. Convert a 1-dimensional vector holding values at grid points into a function.
# V = the vector.
# xMin = x-position corresponding to first entry of V.
# stepSize = step size between grid points.
##
VecFun=function(V,xMin,stepSize){
    return(function(x){
        row=round((x-xMin)/stepSize)+1;
        return(V[row]);
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

        return(mapply(function(r,c){
            d=designMat[cbind(r,c)];
            if (d=="Null") { Null(); }
            else if (d=="FullEnrichment") { FullEnrichment(nSMat[r,c]); }
            else if (d=="PartialEnrichment") { PartialEnrichment(nSMat[r,c],nSCMat[r,c]); }            
        },row,col));
    });
}

##
# VecDesignFun. Construct a function returning the optimal stage 2 design, which depends on the stage 1 data.
# designVec = Vector holding the design type.
# nSVec = Optimal sample size for subgroup S, if applicable.
# nSCVec = Optimal sample size for subgroup SC, if applicable.
# xMin = x-position corresponding to first entry of the vectors.
# yMin = y-position corresponding to first column of the vectors.
# stepSize = step size between grid points (in both dimensions).
##
VecDesignFun=function(designVec,nSVec,nSCVec,xMin,stepSize){
    return(function(x){
        row=round((x-xMin)/stepSize)+1;

        return(mapply(function(r){
            d=designVec[r];
            if(d=="Null") { NullDesign(); }
            else if(d=="FullEnrichment") { FullEnrichment(nSVec[r]); }
            else if(d=="PartialEnrichment") { PartialEnrichment(nSVec[r],nSCVec[r]); }            
        },row));
    });
}

##
# Stage1EU. Compute the stage 1 expected utility for a given stage 1 design and return the optimal stage 2 policy.
##
Stage1EU=function(design1,...) {
    UseMethod("Stage1EU",design1);
}

Stage1EU.Null=function(design1,dp,params){
    out=OptStage2EU(design1,dp,params);

    return(list(eu1=out$eu2,
                eu2=out$eu2,
                design2=out$design2,
                zSmin=0,zSmax=0));
}

Stage1EU.FullEnrichment=function(design1,dp,params){
    ## Construct a grid containing all sufficiently probable values of zS1, as deltaS varies over the prior points.    
    n1=design1$nS;
        
    hw=GridHalfWidth(dp$epsilon1);
    
    priorMeansZS=dp$prior[,1]*sqrt(n1/params$sigma^2);
        
    globalMaxStepZS=ceiling((max(priorMeansZS)+hw)/dp$stepSize1);
    globalMinStepZS=floor((min(priorMeansZS)-hw)/dp$stepSize1);
    
    vecLength=globalMaxStepZS-globalMinStepZS+1;   
        
    eu2Vec=rep(NA,vecLength);
    designVec=rep(NA,vecLength);
    nSVec=rep(NA,vecLength);
    nSCVec=rep(NA,vecLength);
    
    nrowPrior=nrow(dp$prior);
    
    ## Compute expected utility given true effect sizes and then sum over prior
    eu1=0;
    for(k in 1:nrowPrior) {
        deltaS=dp$prior[k,1];      
        priorProb=dp$prior[k,3];                
        
        maxStepZS=ceiling((priorMeansZS[k]+hw)/dp$stepSize1);
        minStepZS=floor((priorMeansZS[k]-hw)/dp$stepSize1);        
        
        is=minStepZS:maxStepZS;            
        
        condEU=0;
        for(i in is) {           
            zS1=i*dp$stepSize1;
            
            ## Check if entry has already been computed. If not, compute and update vectors for stage 2 policy. 
            ii=i-globalMinStepZS+1; ## Translate indices so they start at 1
            
            if(is.na(eu2Vec[ii])) {                                          
                out=OptStage2EU(design1,dp,params,zS1);
                
                designType=class(out$design2);
                eu2Vec[ii]=out$eu2;                                        
                designVec[ii]=designType;

                if(designType=="FullEnrichment") {
                    nSVec[ii]=out$design2$nS;
                }
                
                if(designType=="PartialEnrichment") {
                    nSVec[ii]=out$design2$nS;
                    nSCVec[ii]=out$design2$nSC;
                }                                                           
            }                    
            
            condEU=condEU+eu2Vec[ii]*SampleDensity1(deltaS,zS1,n1,params$sigma)
        }

        ## Get normalising constant for the sample density given deltaS
        normConst=sum(SampleDensity1(deltaS,dp$stepSize1*is,n1,params$sigma));
        
        eu1=eu1+priorProb*condEU/normConst;
    }
    
    ## Construct functions for the stage 2 optimal EU and the optimal policy    
    zSmin=globalMinStepZS*dp$stepSize1;
    zSmax=globalMaxStepZS*dp$stepSize1;    
    
    return(list(eu1=eu1,
                eu2=VecFun(eu2Mat,zSmin,dp$stepSize1),
                design2=VecDesignFun(designVec,nSVec,nSCVec,zSmin,dp$stepSize1),              
                zSmin=zSmin,zSmax=zSmax));
}

Stage1EU.PartialEnrichment=function(design1,dp,params){
    ## Construct a grid containing all sufficiently probable pairs (zS1, zSC1),
    ## as (deltaS, deltaSC) varies over the prior points.
    
    n1=design1$nS+design1$nSC;
    lambda1=design1$nS/n1;
        
    hw=GridHalfWidth(dp$epsilon1);  

    priorMeansZS=dp$prior[,1]*sqrt(lambda1*n1/params$sigma^2);
    priorMeansZSC=dp$prior[,2]*sqrt((1-lambda1)*n1/params$sigma^2);
        
    globalMaxStepZS=ceiling((max(priorMeansZS)+hw)/dp$stepSize1);
    globalMinStepZS=floor((min(priorMeansZS)-hw)/dp$stepSize1);
    globalMaxStepZSC=ceiling((max(priorMeansZSC)+hw)/dp$stepSize1);
    globalMinStepZSC=floor((min(priorMeansZSC)-hw)/dp$stepSize1);

    matRows=globalMaxStepZS-globalMinStepZS+1;
    matCols=globalMaxStepZSC-globalMinStepZSC+1;
        
    eu2Mat=matrix(nrow=matRows,ncol=matCols);
    designMat=matrix(nrow=matRows,ncol=matCols);
    nSMat=matrix(nrow=matRows,ncol=matCols);
    nSCMat=matrix(nrow=matRows,ncol=matCols);
    
    nrowPrior=nrow(dp$prior);
    
    ## Compute expected utility given true effect sizes and then sum over prior
    eu1=0;
    for(k in 1:nrowPrior) {
        deltaS=dp$prior[k,1];
        deltaSC=dp$prior[k,2];
        priorProb=dp$prior[k,3];                
        
        maxStepZS=ceiling((priorMeansZS[k]+hw)/dp$stepSize1);
        minStepZS=floor((priorMeansZS[k]-hw)/dp$stepSize1);
        maxStepZSC=ceiling((priorMeansZSC[k]+hw)/dp$stepSize1);
        minStepZSC=floor((priorMeansZSC[k]-hw)/dp$stepSize1);

        is=minStepZS:maxStepZS;
        js=minStepZSC:maxStepZSC;       
        
        condEU=0;
        for(i in is) {
            for(j in js) {
                zS1=i*dp$stepSize1;
                zSC1=j*dp$stepSize1;                                
                
                ## Check if entry has already been computed. If not, compute and update matrices for stage 2 policy. 
                ii=i-globalMinStepZS+1; ## Translate indices so they start at 1
                jj=j-globalMinStepZSC+1;
                
                if(is.na(eu2Mat[ii,jj])) {                                          
                    out=OptStage2EU(design1,dp,params,zS1,zSC1);

                    designType=class(out$design2);
                    eu2Mat[ii,jj]=out$eu2;                                        
                    designMat[ii,jj]=designType;

                    if(designType=="FullEnrichment") {
                        nSMat[ii,jj]=out$design2$nS;
                    }
                    
                    if(designType=="PartialEnrichment") {
                        nSMat[ii,jj]=out$design2$nS;
                        nSCMat[ii,jj]=out$design2$nSC;
                    }
                }                    
                
                condEU=condEU+eu2Mat[ii,jj]*SampleDensity(deltaS,deltaSC,zS1,zSC1,n1,lambda1,params$sigma);
            }
        }

        ## Get normalising constant for the sample density given (deltaS,deltaSC)
        normConst=sum(outer(dp$stepSize1*is,dp$stepSize1*js,
            function(x,y) SampleDensity(deltaS,deltaSC,x,y,n1,lambda1,params$sigma)));        
        
        eu1=eu1+priorProb*condEU/normConst;
    }    

    ## Construct functions for the stage 2 optimal EU and the optimal policy    
    zSmin=globalMinStepZS*dp$stepSize1;
    zSmax=globalMaxStepZS*dp$stepSize1;
    zSCmin=globalMinStepZSC*dp$stepSize1;
    zSCmax=globalMaxStepZSC*dp$stepSize1;
    
    return(list(eu1=eu1,
                eu2=MatFun(eu2Mat,zSmin,zSCmin,dp$stepSize1),
                design2=MatDesignFun(designMat,nSMat,nSCMat,zSmin,zSCmin,dp$stepSize1),              
                zSmin=zSmin,zSmax=zSmax));
}
    
##
# OptStage2EU. Optimise the stage 2 expected utility.
# Either summation or integration is used to compute the stage 2 EU,
# as indicated by the type of the decision problem (type = "Sum" or type = "Int").
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
     
        if (class(d)=="FullEnrichment") {

            if (dp$type=="Sum") {
                ## Create grid for zS2
                hw=GridHalfWidth(dp$epsilon2); ## Half width of grid square
                Gz2=seq(-hw,hw,dp$stepSize2);
                
                ## Compute stage 2 expected utility for each stage 2 design         
                for (nS in d$nS) {
                    design2=FullEnrichment(nS);
                    out=Stage2EUSum(design1,design2,dp,params,post,Gz2,...);
                    
                    if(out>opts$eu2){
                        opts=optimals(out,design2);
                    }  
                }
                
            } else if (dp$type=="Int") {
                ## Create wrapper function for applying optim
                f=function(x){
                    design2=FullEnrichment(x[1]);
                    Stage2EUInt(design1,design2,dp,params,post,...)
                }
              
                ## Set start point for local optima search
                start=(min(d$nS)+max(d$nS))/2;
                
                ## Find optimal combination of n2S and n2SC and return the result
                opt=optim(start,fn=f,method="L-BFGS-B",
                    lower=min(d$nS),upper=max(d$nS),control=list(fnscale=-1));                 
                
                if(opt$value>opts$eu2){                
                    opts=optimals(opt$value,FullEnrichment(opt$par[1]));
                }
            }          
        }
        
        if (class(d)=="PartialEnrichment") {
            
            if (dp$type=="Sum") {
                ## Create grid for (zS2,zSC2)
                hw=GridHalfWidth(dp$epsilon2); ## Half width of grid square
                Gz2=as.matrix(expand.grid(seq(-hw,hw,dp$stepSize2),seq(-hw,hw,dp$stepSize2)));
                colnames(Gz2)=NULL;
            
                ## Compute stage 2 expected utility for each stage 2 design         
                for (nS in d$nS) {
                    for (nSC in d$nSC) {
                        design2=PartialEnrichment(nS,nSC);
                        out=Stage2EUSum(design1,design2,dp,params,post,Gz2,...);

                        if(out>opts$eu2){
                            opts=optimals(out,design2);
                        }  
                    }
                }
                
            } else if (dp$type=="Int") {
                ## Create wrapper function for applying optim
                f=function(x){
                    design2=PartialEnrichment(x[1],x[2]);
                    Stage2EUInt(design1,design2,dp,params,post,...)
                }
                
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
# Stage2EUSum. Compute the stage 2 expected utility given the stage 1 data using summation.
# sd = Discretised sample distribution for zS2 or (zS2, zSC2) with center 0 or (0,0).
# Each row in the submatrix consisting of the first (or first two) columns represents a possible point.
# The last column in the same row assigns the probability for the point.
##
Stage2EUSum=function(design1,design2,dp,params,post,...){
    UseMethod("Stage2EUSum",design2);
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
