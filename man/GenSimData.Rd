\name{GenSimData}
\alias{GenSimData}
\title{
A function generate Artificial Simulation Dataset.
}
\description{
This is a function used to generate simulation data. In the simulation data, number of phenotypes, number of samples in each phenotype and number of significant CpGs can be assigned. The function randomly generates significant differential CpGs in each phenotype. This dataset can be used to test bgNMF function and EpiCluster.
}
\usage{
GenSimData(Ncpg = 10000, Npheno = 3, Nsample = 10, alpha_p = 1,
         beta_p = 3, Nsig = 1000)
}
\arguments{
  \item{Ncpg}{
Number of total CpGs will be generated. Number of rows for generated Beta matrix. The default number is 10000.
}
  \item{Npheno}{
Number of PhenoTypes in artificial dataset. Note that number of total samples is the number of phenotypes multiply number of samples in each phenotype. Each phenotype will share same number of sample. The default number is 3.
}
  \item{Nsample}{
Number of samples in each PhenoTypes in aritificial dataset. The default number is 10 which means each phenotype will contain 10 samples.
}
  \item{alpha_p}{
One parameter control the beta distribution of artificial simulation data, which will be used to generated beta-distributed data as rbeta(N,alpha_p,beta_p). The default number is 1.
}
  \item{beta_p}{
one parameter control the beta distribution of artificial simulation data, which will be used to generate beta-distributed data as rbeta(N,alpha_p,beta_p). The default number is 3.
}
  \item{Nsig}{
Number of significant CpGs in each PhenoTypes among all Simulation Dataset. The default number is 1000. Normally, number of significant CpGs should be around 10\% of total CpG.
}
}
\value{
A list will be returned, which contain three following information. This simulation data was merely used to test EpiCluster package.
\item{beta}{The aritificial simulation beta matrix generated from this function, which can be used to do EpiCluster clustering.}
\item{SigCpG}{A list recording significant CpGs selected in each phenotype. Users may use this to estimate optimised ee value and effect of bgNMF.}
\item{pheno.v}{PhenoTypes for each sample. This maybe used in EpiDraw function or EpiAnalysis function.}
}
\references{
Yuan T, Ma Z, Beck S, Teschendorff AE. (2015). A fast variational Bayes dimensional reduction and clustering algorithm for Epigenome-Wide Association Studies (EWAS). Under Review.
}
\author{
Yuan Tian, Zhanyu Ma, Andrew Teschendorff
}
\examples{
    Data <- GenSimData(Ncpg=20000,Npheno=5,Nsig=1200)

    Data <- GenSimData(Ncpg=10000,Npheno=3,Nsample=30,Nsig=1000)
}
\keyword{ Simulation }
