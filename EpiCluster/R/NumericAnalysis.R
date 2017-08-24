NumericAnalysis <- function(Basis,
                            clusters,
                            covariant,
                            threshold,
                            maxlevel,
                            cov.name="This")
{
    cat(paste(cov.name,"is a numeric covariate.\n (1) Correlation Test will
              be conducted between this covariate and",
              dim(Basis)[2],"estimated components each.\n
              (2) Krustal Test/ANOVA will be conducted for",
              sum(table(clusters)>=threshold),
              "clusters contain more than",threshold,"samples.\n"))
    cat("---------------------------------------\n\n")

    cat(paste("(1) Correlation between each components to",cov.name,"\n"))
    cor.spearman <- apply(Basis,2,
            function(x) cor.test(x,covariant,method="spearman")$p.value)
    cor.pearson <- apply(Basis,2,
            function(x) cor.test(x,covariant,method="pearson")$p.value)
    cat(paste("There are",
              as.character(length(which(cor.spearman<0.05))),
              "Components significantly Correlate with",cov.name,":\n"))
    correlation.output <- data.frame(cor.spearman=cor.spearman,
                                     cor.pearson=cor.pearson)
    cat("\n")
    ## ==============================================================
    Clusters_Above_Threshold <- clusters[clusters %in% 
                        names(which(table(clusters) >= threshold))]
    Covariate_Corresponding_To_Cluster <- covariant[clusters %in% 
                        names(which(table(clusters) >= threshold))]

    cat(paste("(2) Krustal Test/ANOVA for",
              length(table(Clusters_Above_Threshold)),
              "clusters on",cov.name,"\n"))
    cat("---------------------------------------\n")
    cat(paste("After Filering, there are",
              length(table(Clusters_Above_Threshold)),
              "clusters(based on maxlevel",maxlevel,") contain",
              threshold,"or more Samples, 
              we will only do Krustal Test and ANOVA  on these",
              length(table(Clusters_Above_Threshold)),
              "clusters. Corresponding to these clusters, we get",
              length(Covariate_Corresponding_To_Cluster),"samples.\n"))
    cat("---------------------------------------\n")
    Krustal <- kruskal.test(Clusters_Above_Threshold,
                            Covariate_Corresponding_To_Cluster)
    cat("Krustal Test:\n")
    print(Krustal)
    cat("ANOVA Test:\n")
    AOV <- aov(Covariate_Corresponding_To_Cluster ~ Clusters_Above_Threshold)
    print(AOV)
    print(summary(AOV))

    return(list(cor.spearman=cor.spearman,
                cor.pearson=cor.pearson,
                Krustal=Krustal,
                AOV=AOV,
                ChisquareTest=NULL))
}
