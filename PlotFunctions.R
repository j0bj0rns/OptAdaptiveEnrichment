## PlotFunctions.

## Contains functions for plotting the results of the optimization.

## Generalize the contour plot functions of Thomas for stage 2 and stage 1
PlotStage2ContourSponsor=function(){
}

PlotStage2ContourPublicHealth=function(){
}

## Plot the optimal stage 2 EU, sample size and prevalence against the population prevalence,
## for both the sponsor and the public health viewpoint.
PlotStage2Optimals=function(){
}

PlotStage1ContourSponsor=function(){
}

PlotStage1ContourPublicHealth=function(){
}

## Plot the optimal stage 1 EU, sample size and prevalence against the population prevalence,
## for both the sponsor and the public health viewpoint.
PlotStage1Optimals=function(){
}

PlotPolicy=function(policy,zSC1s){

    zS1s=seq(policy$zSmin,policy$zSmax,policy$stepSize1)

    ## Plot optimal stage 2 EU as a function of zS1
    eu2max=max(sapply(zSC1s,function(zSC1) max(policy$eu2(zS1s,zSC1))));
    eu2min=min(sapply(zSC1s,function(zSC1) min(policy$eu2(zS1s,zSC1))));
    x11();
    plot(c(policy$zSmin,policy$zSmax),c(eu2min,eu2max),
         type="n",xlab=expression(z[S]^1),ylab="Stage 2 Optimal EU");
   
    for(i in 1:length(zSC1s)){
        lines(zS1s,policy$eu2(zS1s,zSC1s[[i]]),lty=i);    
    }  

    legend("bottomleft",
           legend=sapply(zSC1s, function(zSC1) as.expression(substitute(zSC1 == x, list(x = as.name(zSC1))))),
           lty=1:length(zSC1s),
           bg="white")   

    ## Plot optimal stage 2 sample size as a function of zS1
    n2max=max(sapply(zSC1s,function(zSC1) max(policy$n2(zS1s,zSC1))));
    n2min=min(sapply(zSC1s,function(zSC1) min(policy$n2(zS1s,zSC1))));
    x11();
    plot(c(policy$zSmin,policy$zSmax),c(n2min,n2max),
         type="n",xlab=expression(z[S]^1),ylab="Stage 2 Optimal Sample Size");
   
    for(i in 1:length(zSC1s)){
        lines(zS1s,policy$n2(zS1s,zSC1s[[i]]),lty=i);   
    }

    legend("bottomleft",
           legend=sapply(zSC1s, function(zSC1) as.expression(substitute(zSC1 == x, list(x = as.name(zSC1))))),
           lty=1:length(zSC1s),
           bg="white")
    
    ## Plot optimal stage 2 prevalence as a function of zS1
    lambda2max=max(sapply(zSC1s,function(zSC1) max(policy$lambda2(zS1s,zSC1))));
    lambda2min=min(sapply(zSC1s,function(zSC1) min(policy$lambda2(zS1s,zSC1))));
    x11(); 
    plot(c(policy$zSmin,policy$zSmax),c(lambda2min,lambda2max),
         type="n",xlab=expression(z[S]^1),ylab="Stage 2 Optimal Prevalence");
   
    for(i in 1:length(zSC1s)){
        lines(zS1s,policy$lambda2(zS1s,zSC1s[[i]]),lty=i);
    }

    legend("bottomleft",
           legend=sapply(zSC1s, function(zSC1) as.expression(substitute(zSC1 == x, list(x = as.name(zSC1))))),
           lty=1:length(zSC1s),
           bg="white")
}
