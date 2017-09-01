\name{BGNMF}
\alias{BGNMF}
\title{
Run Beta-Gamma Non-negative Matrix Factorization (bgNMF) on a DNA methylation dataset.
}
\description{
This functions performs a dimensional reduction of a DNA methylation data matrix using a Beta-Gamma Non-negative Matrix Factorization (bgNMF) algorithm which preserves the bounded nature of the data. The user can adjust the number of iterations, a sparsity parameter (EE) and the number of latent variables to infer. The function BGNMF is an internal function of this package. The output of this function is used by the other functions provided in this package.
}
\usage{
BGNMF(Data, nIter = 60, K = NULL, EE = 1e-04)
}
\arguments{
  \item{Data}{
DNA methylation data matrix, with number of rows labeling features (e.g. CpGs/Genes), and columns labeling samples. NA values are NOT allowed, so imputation is necessary before application of this function. Also ExpressionSet object is supported by this function. If ExpressionSet is detected, EpiCluster will automatically extract data matrix contained in this dataset.
}
  \item{nIter}{
Number of iterations used in the optimisation. Default number is 60. The smaller this number the less accurate the final estimate, but larger values of this parameter will lead to longer running times.
}
  \item{K}{
Number of latent variables to infer. If not known (usual scenario), the function will use Random Matrix Theory to estimate the number of latent variables.
}
  \item{EE}{
Sparsity parameter used to control the degree of sparsity in the NMF procedure. A smaller EE will lead to a higher sparsity level. The default parameter is 0.0001. Note that sparsity is imposed on features.
}
}
\details{
Normally this function will not be used by users, since it is called by the EpiCluster() function.
}
\value{
A list will be returned, which contain 5 slots:
\item{X_reconstruct}{X_reconstruct is the factorized approximation to the original data matrix, obtained by multiplying the estimated basis matrix alphaMatrix/(alphaMatrix+betaMatrix) and excitation matrix excitationMatrix.}
\item{alphaMatrix}{alphaMatrix is one part of the basis matrix, defining the alpha parameter of the basis matrix, which is beta-valued.}
\item{betaMatrix}{betaMatrix is the other part of the basis Matrix, defining the beta parameter of the basis Matrix. So, estimated basis matrix is alphaMatrix/(alphaMatrix+betaMatrix).}
\item{excitationMatrix}{excitationMatrix. is the excitation matrix, sparsity of which is controled by parameter EE.}
\item{ee}{ee is the sparsity parameted EE used by bgNMF.}
}
\references{
Ma, Z., Teschendorff, A. E., Leijon, A., Qiao, Y., Zhang, H., and Guo, J. (2014b). Variational bayesian matrix factorization for bounded support data. IEEE Trans. on Pattern Analysis and Machine Intelligence, 35(4), 876-889.

Yuan T, Ma Z, Beck S, Teschendorff AE. (2015). A fast variational Bayes dimensional reduction and clustering algorithm for Epigenome-Wide Association Studies (EWAS). Under Review.

Teschendorff AE, Zhuang J, Widschwendter M. (2011). Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray studies. Bioinformatics. 2011 Jun 1;27(11):1496-505.
}
\author{
Yuan Tian, Zhanyu Ma, Andrew Teschendorff
}
\note{
NA value are NOT allowed, so imputation may need to be conducted to methylation dataset beforehead.
}

\examples{
    Data <- GenSimData(Ncpg=1000,Nsig=100)
    BG.Result <- BGNMF(Data$beta,nIter=10,EE=0.005)
}
\keyword{ bgNMF }
