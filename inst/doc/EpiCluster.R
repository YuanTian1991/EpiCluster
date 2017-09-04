## ----eval=FALSE----------------------------------------------------------
#  install.package("RPMM")
#  install.package("Rlab")
#  install.package("isva")
#  install.package("pheatmap")

## ----eval=FALSE----------------------------------------------------------
#  R CMD INSTALL EpiCluster_1.0.tar.gz

## ----eval=FALSE----------------------------------------------------------
#  library("EpiCluster")

## ----eval=FALSE----------------------------------------------------------
#  # Here we generated a Artificial Simulation Data, which contain 10000 CpGs and 4 classes/phenotypes. Each class contains 20 samples, and 1000 class-specific DMCs will be generated.
#  SimData <- GenSimData(Ncpg=10000,Nsig=1000,Npheno=4,Nsample=20)
#  library("pheatmap")
#  pheatmap(SimData$beta)

## ---- out.width = 800,out.height = 350, fig.retina = NULL,echo=F---------
knitr::include_graphics("./Figure1.pdf")

## ----eval=FALSE----------------------------------------------------------
#  SimData <- GenSimData(Ncpg=1000,Nsig=100,Npheno=4,Nsample=30)
#  # Here we conducted EpiCluster on simulation generated above, we set number of iteration as 20, and EE parameter as 0.005, which is relatively big so more significant CpGs shall be find after EpiCluster. We did not assigned K value here so EpiCluster will automaticly detect number of latent viariables, based on Random Matrix Theory in isva pacakge.
#  EpiCluster.Result <- EpiCluster(SimData$beta,nIter=20,EE=0.005)
#  
#  [1] "Iteration times:  20"
#  [1] "Components number:  4"
#  [1] "ee parameter:  0.005"
#  [1] 1
#  [1] 2
#  [1] 3
#  [1] 4
#  [1] 5
#  [1] 6
#  [1] 7
#  [1] 8
#  [1] 9
#  [1] 10
#  [1] 11
#  [1] 12
#  [1] 13
#  [1] 14
#  [1] 15
#  [1] 16
#  [1] 17
#  [1] 18
#  [1] 19
#  [1] 20
#  [1] "Run BGNMF success!"
#  [1] "root"
#  [1] Inf
#  [1] 0.06998368
#  [1] 0.05706815
#  ...
#  [1] 3.154375e-07
#  [1] 365.676
#  [1] "rL"
#  [1] Inf
#  [1] 0.06250856
#  [1] 0.005655091
#  [1] 2.249948e-09
#  [1] 368.4922
#  [1] "rLL"
#  [1] Inf
#  [1] 0.006686889
#  [1] 0.0100536
#  ...
#  [1] 0.0002950578
#  [1] 208.3635
#  [1] "rLR"
#  [1] Inf
#  [1] 0.00893239
#  [1] 0.009745804
#  ...
#  [1] 2.059788e-05
#  [1] 216.4498
#  [1] "rR"
#  [1] Inf
#  [1] 0.06133592
#  [1] 5.436303e-05
#  [1] 0
#  [1] 380.732
#  [1] "rRL"
#  [1] Inf
#  [1] 0.002652312
#  [1] 0.003085925
#  ...
#  [1] 0.00085122
#  [1] 218.0377
#  [1] "rRR"
#  [1] Inf
#  [1] 0.004746239
#  [1] 0.007349823
#  ...
#  [1] 217.5082
#  [1] "Run RPBMM success!"
#  The simulated data looks like below:

## ---- out.width = 700, out.height = 700, fig.retina = NULL,echo=F--------
knitr::include_graphics("./Figure2.pdf")

## ----eval=FALSE----------------------------------------------------------
#  slotNames(EpiCluster.Result)
#  [1] "bgNMF" "Cluster" "betaRPMM"

## ----eval=FALSE----------------------------------------------------------
#  # Here we only input one covariate into Phenotypes parameter because it's not easy to contruct multi-covariate simulation data. We will show a more comprehensive samples later. Also, in this example, we ignored parameter maxlevel here, which means all clusters from EpiCluster will be analysed.
#  EpiAnalysis.Result <- EpiAnalysis(EpiCluster.Result,PhenoTypes=SimData$pheno.v,threshold=10)
#  
#  ========= Only one Covariate in PhenoType =========
#  == All Analysis Result will be returned in output ==
#  --------------------- START ------------------------
#  
#  =============== EpiCluster Result ===============
#  bgNMF Detected 4 Components.
#  bgNMF Detected 4 Clusters under maxlevel Inf .
#  There are 4 Clusters contain more than 10 Samples in it.
#  
#  ---------------------------------------
#  PhenoType is a categorical covariate.
#       (1) ANOVA Test will be conducted between this covariate and 4 estimated components each.
#       (2) Chisquare Test will be conducted for 4 clusters contain more than 10 samples.
#  ---------------------------------------
#  
#  (1) ANOVA test between each components to PhenoType
#  There are 4 Components show significance to PhenoType :
#  
#  (2) Chisquare Test for 4 clusters on PhenoType
#  ---------------------------------------
#  After Filering, there are 4 clusters(based on maxlevel Inf ) contain 10 or more Samples, we will only do Chisquare Test on these 4 clusters. Corresponding to    these clusters, we get 120 samples.
#  ---------------------------------------
#  Chisquare Test:
#                                    Clusters_Above_Threshold
#  Covariate_Corresponding_To_Cluster rLLL rLLR rLR rR
#                                   1   30    0   0  0
#                                   2    0    0   0 30
#                                   3    0    0  30  0
#                                   4    0   30   0  0
#  
#          Pearson's Chi-squared test
#  
#  data:  table(Covariate_Corresponding_To_Cluster, Clusters_Above_Threshold)
#  X-squared = 360, df = 9, p-value < 2.2e-16

## ----eval=FALSE----------------------------------------------------------
#  EpiDraw(EpiCluster.Result,PhenoTypes=SimData$pheno.v)

## ---- out.width = 750, out.height = 1200, fig.retina = NULL,echo=F-------
knitr::include_graphics("./Figure4.pdf")

## ----eval=FALSE----------------------------------------------------------
#  names(EpiAnalysis.Result)
#  [1] "PMatrix"  "Analysis"

## ----eval=FALSE----------------------------------------------------------
#  EpiDraw(EpiCluster.Result,PhenoTypes=PhenoTypes.lv,maxlevel=3,threshold=10)

## ---- out.width = 700, out.height = 700, fig.retina = NULL,echo=F--------
knitr::include_graphics("./Figure7.pdf")

## ---- out.width = 800, out.height = 500, fig.retina = NULL,echo=F--------
knitr::include_graphics("./Figure8.pdf")

