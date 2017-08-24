\name{EpiDraw}
\alias{EpiDraw}
\title{
Draw basis matrix heatmap and clustering result.
}
\description{
This function generates a graphical representation of the RPBMM clustering result on the basis matrix inferred using EpiCluster. In more detail, the function generates a heatmap of the basis matrix, with the clusters indicated by colors from rainbow(n) function, where n represents the numbers of clusters. The numbers of clusters may be controlled by the paramater \emph{maxlevel}. In addition, function can generate boxplots or barplots to illustrate correlations of latent variables or clustering result with only one covariate, which is depended on types of covariate. For more analysis between components and covariates, one may turn to EpiAnalysis() function for help.
}
\usage{
EpiDraw(EpiCluster.output, PhenoTypes = NULL, maxlevel = NULL,threshold=10,
        legend = TRUE, samplecluster = TRUE, maintitle = "",
        colors = c("yellow", "navy"),clusterplot = TRUE,ComponentMatrix = TRUE)
}
\arguments{
  \item{EpiCluster.output}{
An EpiClusterResult object from EpiCluster. Must be provided. The bgNMG slot in EpiCluster.output S4 object contains paramater matrixes A and B, which can be used to approximately reconstruct Basis Matrix. The Basis Matrix contains value of components(latent variables), which is also the principle of RPMM unsupervised clustering.
}
  \item{PhenoTypes}{
A phenotypes list or vector object, input as for EpiAnalysis function. If it's a list, a p value heatmap will be ploted between each component and each covariate, each vector in the list will be treated as a covariate and will be ploted one by one in boxplots (for numeric covariates) or barplots (for category covariates).
}
  \item{maxlevel}{
Max level (height) of the binary clustering tree generated from RPBMM. An integer parameter. The default value is NULL meaning that all branchs of clusters will be shown and no assembling will be conducted on clusters.
}
  \item{threshold}{
After EpiCluster function, RPBMM may return clusters will different contained samples (even after combination with maxlevel parameter above), if there is only very little samples in one cluster, it may be pointless to do analysis(Krustal or Chisquare Test) between this cluster and various covariates. Thus user may set minimum threshold here to filter clusters contain numbers of samples less than this threshold. The default threshold is 10, which means, the function will automatically ignore all clusters(even after combination of maxlelve parameter above) contain numbers of sample less than 10.

}
  \item{legend}{
A boolean parameter to specify whether the legend of heatmap should be drawn. TRUE or FALSE. The default value is TRUE.
}
  \item{samplecluster}{
A boolean parameter to specify if a bar indicating the cluster assignment should be drawn under the heatmap. The default value is TRUE. The colorbar for cluster will be indicated by rainbow() function, just like EpiAnalysis function.
}
  \item{maintitle}{
A character parameter specifying the name of the plot.
}
  \item{colors}{
A two character parameter specifying the colors used in heatmap, the default is "yellow" and "navy". Users may specify two colors, like c("green","red"), the first one will be assigned to lowest value and the later one will be assigned to highest value.
}
  \item{clusterplot}{
A boolean parameter to specify if component plots should be drawn. The default is TRUE if PhenoType is provided.
}
  \item{ComponentMatrix}{
A boolean parameter to specify if heatmap will be drawn. The default if TRUE.
}
}
\details{
For numeric covariates, boxplot of each cluster will be draw, while for category PhenoTypes, barplot will be shown.
}
\value{
No output at all. Only plots.
}
\references{
Yuan T, Ma Z, Beck S, Teschendorff AE. (2015). A fast variational Bayes dimensional reduction and clustering algorithm for Epigenome-Wide Association Studies (EWAS). Under Review.
}
\author{
Yuan Tian, Zhanyu Ma, Andrew Teschendorff
}
\examples{
    DataSet <- GenSimData(Ncpg=1000,Nsig=100,Npheno=5)
    EpiCluster.Result <- EpiCluster(DataSet$beta,nIter=20)
    EpiDraw(EpiCluster.output=EpiCluster.Result,PhenoTypes=DataSet$pheno.v,
            maintitle="Category DataSet",maxlevel=3,threshold=10)
}
