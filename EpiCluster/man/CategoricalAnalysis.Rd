\name{CategoricalAnalysis}
\alias{CategoricalAnalysis}
\title{
    Correlation Analysis for each categorical covariate.
}
\description{
This function provide ANOVA analysis for the K inferred basis components, each a vector of length equal to the number of samples, to categorical factors. These factors could represent categorical phenotypes (e.g. tissue type etc) or technical effects (batch, chip etc). ANOVA will be applied to categorical covariates. Thus users could use this function to detect if estimated componetes have effect on each covariates. Apart from that, this function also provide enrichment analysis between clusters from RPBMM and each covariate, for this analysis we will use Chisquare Test for categorical covariates.
}
\usage{
CategoricalAnalysis(Basis,clusters,covariant,threshold,maxlevel,
                    cov.name="This")
}
\arguments{
    \item{Basis}{
        Basis matrix from bgNMF function, which contain information of each estimated components.
        }     
    \item{clusters}{
        clusters detected by EpiCluster function, shall be gathered already by maxlevel value.
        }     
    \item{covariant}{
        A vector of character or factor value, indicating one covariate of these samples, for example tissue-type.
        }   
    \item{threshold}{
        After EpiCluster function, RPBMM may return clusters will different contained samples (even after combination with maxlevel parameter above), if there is only very little samples in one cluster, it may be  pointless to do analysis(Krustal or Chisquare Test) between this cluster and various covariates. Thus user may set minimum threshold here to filter clusters contain numbers of samples less than this threshold.   The default threshold is 10, which means, the function will automatically ignore all clusters(even after combination of maxlevel parameter above) contain numbers of sample less than 10.
        }
    \item{maxlevel}{
        Parameter control levels of RPBMM trees node will be analysis. After EpiCluster method, maybe too many branchs of trees will be generated, so we may need to assemble some of them by higher branch node, which is decided by users' research. The maxlevel parameter determine the maxium levels of binary trees will be specified. The default value for maxlevel is NULL, which means all branchs will be considered as a  cluster, which might result to too many clusters. Thus, users should assign a proper maxlevel according to their aim.
        }
    \item{cov.name}{
        name for this covariate. If not specificed, it will be assigned as "This".
        }
}
\details{
    Normally this function will not be used by user. It will be employed inside EpiAnalysis function automatically.
}
\value{
    Totally, there are five returned value, but for categorical phenotypes, only last two values will be assigned, the first three will be omitted.
    \item{cor.spearman}{Correlation Analysis Result with Spearman method for each component (latent variables). It will be assigned NULL here for categorical PhenoTypes.}
    \item{cor.pearson}{Correlation Analysis Result with Pearson Correlation method for each component (latent variables). It will be assigned NULL here for categorical PhenoTypes.}
    \item{Krustal}{Krustal Test result conducted between filtered clusters and numeric covariates, It will be assigned NULL here for categorical PhenoTypes.}
    \item{AOV}{Anove analysis applied on each component to test if they are specificly detected different phenotypes, this result could be used to indicate phenotype-specified components.}
    \item{ChisquareTest}{Chisquare Test applied between filtered clusters and categorical covariates.}
}
\references{
Yuan T, Ma Z, Beck S, Teschendorff AE. (2015). A fast variational Bayes dimensional reduction and clustering algorithm for Epigenome-Wide Association Studies (EWAS). Under Review.
}
\author{
Yuan Tian, Zhanyu Ma, Andrew Teschendorff
}
\examples{
    DataSet <- GenSimData(Ncpg=1000,Npheno=4,Nsig=100)
    EpiCluster.Result <- EpiCluster(DataSet$beta,nIter=20)
    EpiAnalysis.Result <- EpiAnalysis(EpiCluster.Result
        ,PhenoTypes=DataSet$pheno.v,maxlevel=3,threshold=10)
    bgNMGResult <- getbgNMF(EpiCluster.Result)
    BasisMatrix <- bgNMGResult$alphaMatrix/(bgNMGResult$alphaMatrix+
    bgNMGResult$betaMatrix)
    Cluster <- getCluster(EpiCluster.Result)
    clusters <- substr(Cluster,1,3)
    CategoricalAnalysis(Basis=BasisMatrix,
                        clusters=clusters,
                        covariant=DataSet$pheno.v,
                        maxlevel=3,
                        threshold=10)
}
\keyword{ Analysis}
\keyword{ Categorical}
