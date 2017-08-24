CategoricalAnalysis <- function(Basis,
                                clusters,
                                covariant,
                                threshold,
                                maxlevel,
                                cov.name="This")
{
    cat(paste(cov.name,"is a categorical covariate.\n (1) ANOVA Test will be
              conducted between this covariate and",dim(Basis)[2],"estimated
              components each.\n     (2) Chisquare Test will be conducted for"
              ,sum(table(clusters)>=threshold),"clusters contain more than",
              threshold,"samples.\n"))
    cat("---------------------------------------\n\n")

    cat(paste("(1) ANOVA test between each components to",cov.name,"\n"))
    component.aov <- apply(Basis,2,function(x)
                           summary(aov(x ~ covariant))[[1]][1,5])
    cat(paste("There are",as.character(length(which(component.aov<0.05))),
              "Components show significance to",cov.name,":\n"))
    ANOVA.p.value=component.aov[which(component.aov<0.05)]

    Clusters_Above_Threshold <- clusters[clusters %in%
            names(which(table(clusters) >= threshold))]
    Covariate_Corresponding_To_Cluster <- covariant[clusters %in%
            names(which(table(clusters) >= threshold))]
    cat("\n")

    cat(paste("(2) Chisquare Test for",
              length(table(Clusters_Above_Threshold)),
              "clusters on",cov.name,"\n"))
    cat("---------------------------------------\n")
    cat(paste("After Filering, there are",
              length(table(Clusters_Above_Threshold)),
              "clusters(based on maxlevel",maxlevel,") contain",
              threshold,"or more Samples, we will only do Chisquare 
              Test on these",length(table(Clusters_Above_Threshold)),
              "clusters. Corresponding to these clusters, we get",
              length(Covariate_Corresponding_To_Cluster),"samples.\n"))
    cat("---------------------------------------\n")
                                                                                                            
    cat("Chisquare Test:\n")
    print(t(table(Clusters_Above_Threshold,
                  Covariate_Corresponding_To_Cluster)))
    ChisquareTest <- chisq.test(table(Covariate_Corresponding_To_Cluster,
                                      Clusters_Above_Threshold))
    print(ChisquareTest)
                                                                                                                                                    
    return(list(cor.spearman=NULL,
                cor.pearson=NULL,
                Krustal=NULL,
                AOV=component.aov,
                ChisquareTest=ChisquareTest))
}
