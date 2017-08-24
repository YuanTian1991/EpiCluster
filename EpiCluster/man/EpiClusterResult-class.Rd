\name{EpiClusterResult-class}
\Rdversion{1.1}
\docType{class}
\alias{EpiClusterResult-class}
\alias{getbgNMF,EpiClusterResult-method}
\alias{setbgNMF,EpiClusterResult-method}
\alias{getCluster,EpiClusterResult-method}
\alias{setCluster,EpiClusterResult-method}
\alias{getbetaRPMM,EpiClusterResult-method}
\alias{setbetaRPMM,EpiClusterResult-method}
\alias{getbgNMF}
\alias{setbgNMF}
\alias{getCluster}
\alias{setCluster}
\alias{getbetaRPMM}
\alias{setbetaRPMM}

\title{Class \code{"EpiClusterResult"}}
\description{
A S4 object holds result from EpiCluster function, including result of bgNMF and RPBMM. And this object is specificly used in EpiAnalysis and EpiDraw function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EpiClusterResult", ...)}.
This S4 object contains result of bgNMF (a list contain 5 element), a S3 class object of RPBMM (a blcTree object from RPMM package, need to load RPMM package) and clustering result of EpiCluster. Not that this object is what will be used in EpiDraw() and EpiAnalysis() function. So it's better not to modified its internal structure.
}
\section{Slots}{
  \describe{
    \item{\code{bgNMF}:}{Object of class \code{"list"} Result of bgNMF function, a list object which contains Parameter Matrix and Excitation Matrix.}
    \item{\code{Cluster}:}{Object of class \code{"factor"} clustering result of RPMM, this is a factor vector in which each cluster was named as rL,rR,rLL,rLR..., each cluster is a branch for the whole RPBMM binary tree. For some users may not want so much branch, can use maxlevel parameter in EpiDraw and EpiAnalysis function to assemble some branches into a bigger cluster.}
    \item{\code{betaRPMM}:}{Object of class \code{"ANY"} Special S3 object provided by RPMM package , which contain all information and clustering parameter for each cluster (and each node in clustering binary tree), for more information, users may check RPMM package for help}
  }
}
\section{Methods}{
    \describe{
        \item{getbgNMF}{\code{signature(x = "EpiClusterResult")}: ... }
        \item{getCluster}{\code{signature(x = "EpiClusterResult")}: ... }
        \item{getbetaRPMM}{\code{signature(x = "EpiClusterResult")}: ... }
        \item{setbgNMF}{\code{signature(x = "EpiClusterResult")}: ... }
        \item{setCluster}{\code{signature(x = "EpiClusterResult")}: ... }
        \item{setbetaRPMM}{\code{signature(x = "EpiClusterResult")}: ... }
        }
}

\section{Accessors}{
    In the code snippets below, \code{x} is an EpiClusterResult object.
    \describe{
        \item{}{
            \code{getbgNMF(x)}:
            Return the bgNMF list result in \code{x}.
            }
        \item{}{
            \code{getCluster(x)}:
            Return the Cluster factor vector result in \code{x}.
            }
        \item{}{
            \code{getbetaRPMM(x)}:
            Return the betaRPMM object in \code{x}.
            }
        \item{}{
            \code{setbgNMF(x,value)}:
            Reset betaNMF list result in \code{x}.
            }
        \item{}{
            \code{setCluster(x,value)}:
            Reset Cluster factor vector in \code{x}.
            }
        \item{}{
            \code{setbetaRPMM(x,value)}:
            Reset the betaRPMM object in \code{x}.
            }
    }
}

\references{Yuan T, Ma Z, Beck S, Teschendorff AE. (2015). A fast variational Bayes dimensional reduction and clustering algorithm for Epigenome-Wide Association Studies (EWAS). Under Review.
}
\author{
Yuan Tian, Zhanyu Ma, Andrew Teschendorff
}
\examples{
showClass("EpiClusterResult")

Data <- GenSimData(Ncpg=1000,Nsig=100)
EpiCluster.Result <- EpiCluster(Data$beta,nIter=20,EE=0.005)

bgNMF_result <- getbgNMF(EpiCluster.Result)
Cluster_result <- getCluster(EpiCluster.Result)
betaRPMM_result <- getbetaRPMM(EpiCluster.Result)

EpiCluster.Result <- setbgNMF(EpiCluster.Result,bgNMF_result)
EpiCluster.Result <- setCluster(EpiCluster.Result,Cluster_result)
EpiCluster.Result <- setbetaRPMM(EpiCluster.Result,betaRPMM_result)
}
\keyword{classes}
\keyword{EpiCluster}
