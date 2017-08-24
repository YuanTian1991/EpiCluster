GenSimData <- function(Ncpg=10000,
                       Npheno=3,
                       Nsample=10,
                       alpha_p=1,
                       beta_p=3,
                       Nsig=1000)
{
    SimMatrix <- matrix(rbeta(Ncpg*Npheno*Nsample,
                              alpha_p,beta_p),Ncpg,Npheno*Nsample)

    DoDuplicate <- 0
    sig.index <- list()
    Part <- sample(1:Ncpg,Nsig)
    DoDuplicate <- c(DoDuplicate,Part)
    sig.index[[1]] <- Part
    SimMatrix[Part,1:Nsample] <- rbeta(Nsample*Nsig,beta_p,alpha_p)
    for(i in 2:Npheno)
    {
        Part <- sample(c(1:Ncpg)[-DoDuplicate],Nsig)
        DoDuplicate <- c(DoDuplicate,Part)
        sig.index[[i]] <- Part
        SimMatrix[Part,((i-1)*Nsample+1):(i*Nsample)] <- 
            rbeta(Nsample*Nsig,beta_p,alpha_p)
    }
    pheno.v <- as.factor(rep(1:Npheno,each=Nsample))
    rownames(SimMatrix) <- paste("cg",as.character(c(1:Ncpg)),sep="")
    colnames(SimMatrix) <- paste("sample",
                as.character(c(1:(Nsample*Npheno))),sep="")
    return(list(beta=SimMatrix,SigCpG=sig.index,pheno.v=pheno.v));
}
