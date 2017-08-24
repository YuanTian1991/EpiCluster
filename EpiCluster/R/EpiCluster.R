EpiCluster <- function(X,nIter=60,K=NULL,EE=0.0001)
{
    if(class(X) == "ExpressionSet"){
        X <- exprs(X)
    }else{
        X <- X
    }
    if(length(which(is.na(X)))!=0)
    {
        print("NA values are not allowed.")
        return(FALSE)
    }
    
    bgNMF.output <- BGNMF(Data=X,nIter=nIter,K=K,EE=EE)
    
    print("Run BGNMF success!")
    Basis <- bgNMF.output$alphaMatrix / 
    (bgNMF.output$alphaMatrix+bgNMF.output$betaMatrix)

    betaRPMM <- blcTree(Basis)
    betaRPMMClasses <- blcTreeLeafClasses(betaRPMM)
    print("Run RPBMM success!")

    EpiCluster.output <- new("EpiClusterResult",
                             bgNMF=bgNMF.output,
                             Cluster=betaRPMMClasses,
                             betaRPMM=betaRPMM)
    return(EpiCluster.output)
}
