\name{EpiAnalysis}
\alias{EpiAnalysis}
\title{
Correlation analyse the output from EpiCluster.
}
\description{
This function provide analysis that correlates the K inferred basis components, each a vector of length equal to the number of samples, to factors. These factors could represent phenotypes (e.g. age, tissue type etc) or technical effects (batch, chip etc). Correlations are calculated for numeric covariates, while ANOVA will be applied to categorical covariates. Thus users could use this function to detect if estimated componetes have effect on each covariates. Apart from that, this function also provide enrichment analysis between clusters from RPBMM and each covariate, for this analysis we will use Krustal Test for numeric covariates, and Chisquare Test for categorical covariates. This function will employ NumericAnalysis() and CategoricalAnalysis() functions automatically for different covariates.
}
\usage{
EpiAnalysis(EpiCluster.output, PhenoTypes, maxlevel = NULL,threshold = 10)
}
\arguments{
  \item{EpiCluster.output}{
      Output of EpiCluster function, an EpiClusterResult S4 object containing bgNMF result and clustering result. MUST be provided.
}
  \item{PhenoTypes}{
      A list or vecter object, which entries representing vectors specifying phenotypes of samples. These could be numeric, factors or characters are acceptable. Function will automatically detect which type the phenotype is and perform correlative analysis according to the type. If the input is a vector on certain covariate, the analysis will be conducted on that specific covariate. Users should make sure the phenotypes inputed are "numeric","character" or "factor", or function may not work properly. Also, if the inputted object is a list, we strong recommend user assign name of each covariate beforehead, otherwise, this function will assign all covariates' name as "factor_1,factor_2,factor_3...".
}
  \item{maxlevel}{
      Parameter control levels of RPBMM trees node will be analysis. After EpiCluster method, maybe too many branchs of trees will be generated, so we may need to assemble some of them by higher branch node, which is decided by users' research. The maxlevel parameter determine the maxium levels of binary trees will be specified. The default value for maxlevel is NULL, which means all branchs will be considered as a cluster, which might result to too many clusters. Thus, users should assign a proper maxlevel according to their aim.
}
  \item{threshold}{
      After EpiCluster function, RPBMM may return clusters will different contained samples (even after combination with maxlevel parameter above), if there is only very little samples in one cluster, it may be pointless to do analysis(Krustal or Chisquare Test) between this cluster and various covariates. Thus user may set minimum threshold here to filter clusters contain numbers of samples less than this threshold. The default threshold is 10, which means, the function will automatically ignore all clusters(even after combination of maxlelve parameter above) contain numbers of sample less than 10.
}
}
\details{
    EpiAnalysis function will employ NumericAnalysis function and CategoricalAnalysis function for two different types of covariates.
}
\value{
    if the "Phenotypes" parameter is assigned as a vector, a list contain following five value will be returned. For two different kind of covariates (numeric and category), two different result will be given.
    \item{cor.spearman}{Correlation Analysis Result with Spearman method for each component (latent variables). Numeric covariates only.}
    \item{cor.pearson}{Correlation Analysis Result with Pearson Correlation method for each component (latent variables). Numeric covariates only.}
    \item{Krustal}{Krustal Test result conducted between filtered clusters and numeric covariates,Numeric covariates only.}
    \item{AOV}{For numeric covariates, ANOVA will be conducted between filtered clusters and numeric covariates. But for categorical covariates, ANOVA will be conducted between estimated latent components and each categorical phenotypes.}
    \item{ChisquareTest}{Chisquare Test applied between filtered clusters and categorical covariates. Categorical covariates only.}

    if the "Phenotypes" parameter is assigned as a list. For each covariate, the analysis will be done as above and saved in a list object named Analysis, additionally, a P value matrix will be returned, which indicate the significance p value between each estimated components and each covariates. For numeric covariats, p value are calculated by Spearman Correlation, for categorical covariates, p value are calculated by ANOVA.
    \item{PMatrix}{P value of significance between each component and covariates.}
    \item{Analysis}{A list contain all analysised result for each covariate, each of them will be treated as a "numeric" vector or "category" vector, and the analysis result will be stored into this list.}

}
\references{
Yuan T, Ma Z, Beck S, Teschendorff AE. (2015). A fast variational Bayes dimensional reduction and clustering algorithm for Epigenome-Wide Association Studies (EWAS). Under Review.
}
\author{
Yuan Tian, Zhanyu Ma, Andrew Teschendorff
}
\note{
For all statistic test invloved in this function, some output will be printed, user may use sink() function to record all of them for future analysis.
}
\examples{
    DataSet <- GenSimData(Ncpg=1000,Npheno=4,Nsig=100)
    EpiCluster.Result <- EpiCluster(DataSet$beta,nIter=20)
    EpiAnalysis.Result <- EpiAnalysis(EpiCluster.Result,
                                      PhenoTypes=DataSet$pheno.v,
                                      maxlevel=3,
                                      threshold=10)}
\keyword{ Analysis }
