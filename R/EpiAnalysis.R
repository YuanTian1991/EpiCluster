EpiAnalysis <- function(EpiCluster.output,
                        PhenoTypes,
                        maxlevel=NULL,
                        threshold=10)
{
    if(class(EpiCluster.output)!="EpiClusterResult")
    {
        print(class(EpiCluster.output),
              "is not correct EpiCluster Result Format.")
        return(FALSE)
    }
## Get Basis Matrix and Excitation Matrix
    Basis <- EpiCluster.output@bgNMF$alphaMatrix / 
    (EpiCluster.output@bgNMF$alphaMatrix+EpiCluster.output@bgNMF$betaMatrix)
    Excitation <- EpiCluster.output@bgNMF$excitationMatrix
    Cluster <- EpiCluster.output@Cluster
    if(is.null(maxlevel))
    {
        clusters <- Cluster
        maxlevel <- "Inf"
    }
    else
        clusters <- substr(Cluster,1,maxlevel)

    message("\n=============== EpiCluster Result ===============")
    message(paste("bgNMF Detected",as.character(dim(Basis)[2]),"Components."))
    message(paste("bgNMF Detected",as.character(length(table(clusters))),
              "Clusters under maxlevel",as.character(maxlevel),"."))
    message("There are",sum(table(clusters)>=threshold),
        "Clusters contain more than",threshold,"Samples in it.\n")
## Detect all PhenoTypes here:
    if(mode(PhenoTypes)=="list")
    {
        message("======= More Than one Covariate in PhenoType ======")
        message("===  All covariate will be Analysis one by one ===")
        message("== All Analysis Result will be returned in output ==\n")
        message("--------------------- START ------------------------")
        Analysis <- list()
        if("" %in% names(PhenoTypes))
        {
            names(PhenoTypes) <- paste("Factor_",1:length(PhenoTypes),sep="")
        }
        PMatrix <- matrix(0,dim(Basis)[2],length(PhenoTypes))
        rownames(PMatrix) <- paste("Components",1:dim(Basis)[2],sep="-")
        colnames(PMatrix) <- names(PhenoTypes)
        for(element in 1:length(PhenoTypes))
        {
            cov.name <- names(PhenoTypes)[element]
            covariant <- PhenoTypes[,element]
            message(paste(">>>>>>",cov.name,"<<<<<<"))
            message("---------------------------------------")
            if(class(PhenoTypes[,element])=="numeric")
            {
                Analysis[[cov.name]] <- NumericAnalysis(Basis,
                                                        clusters,
                                                        covariant,
                                                        threshold,
                                                        maxlevel,
                                                        cov.name=cov.name)
                PMatrix[,element] <- Analysis[[cov.name]]$cor.spearman
            }
            else
            {
                Analysis[[cov.name]] <- CategoricalAnalysis(Basis,
                                                            clusters,
                                                            covariant,
                                                            threshold,
                                                            maxlevel,
                                                            cov.name=cov.name)
                PMatrix[,element] <- Analysis[[cov.name]]$AOV
            }
            message("---------------------------------------")
            message("\n")
        }
        return(output=list(PMatrix=PMatrix,Analysis=Analysis))
    }
    else
    {
        message("========= Only one Covariate in PhenoType =========")
        message("== All Analysis Result will be returned in output ==")
        message("--------------------- START ------------------------")
        cov.name <- "PhenoType"
        if(class(PhenoTypes)=="numeric")
        {
            output <- NumericAnalysis(Basis,
                                      clusters,
                                      PhenoTypes,
                                      threshold,
                                      maxlevel,
                                      cov.name=cov.name)
        }
        else
        {
            output <- CategoricalAnalysis(Basis,
                                          clusters,
                                          PhenoTypes,
                                          threshold,
                                          maxlevel,
                                          cov.name=cov.name)
        }
        return(output)
    }
}
