setClass("EpiClusterResult",
         slots=list(bgNMF="list",
                    Cluster="factor",
                    betaRPMM="ANY"))

.valid.EpiClusterResult.bgNMF <- function(x)
{
    if (length(x@bgNMF) != 5)
        return("'bgNMF' must be a 'list' of length 5.")
    if (!all(names(x@bgNMF) %in% c("ee", "X_reconstruct", "alphaMatrix", 
                                   "betaMatrix","excitationMatrix")))
        return("'names(bgNMF)' must be c(\"ee\", \"X_reconstruct\", 
               \"alphaMatrix\", \"betaMatrix\",\"excitationMatrix\")")
    NULL
}
setValidity("EpiClusterResult", .valid.EpiClusterResult.bgNMF)

### Getters.

setGeneric("getbgNMF", function(x) standardGeneric("getbgNMF"))
setMethod("getbgNMF", "EpiClusterResult", function(x) x@bgNMF)

setGeneric("getCluster", function(x) standardGeneric("getCluster"))
setMethod("getCluster", "EpiClusterResult", function(x) x@Cluster)

setGeneric("getbetaRPMM", function(x) standardGeneric("getbetaRPMM"))
setMethod("getbetaRPMM", "EpiClusterResult", function(x) x@betaRPMM)

### Setters.
setGeneric("setbgNMF", function(x,...) standardGeneric("setbgNMF"))
setMethod("setbgNMF", signature(x="EpiClusterResult"),
          function(x,value)
          {
              if(length(value)==5 & names(value) == c("ee",
                                                      "X_reconstruct",
                                                      "alphaMatrix",
                                                      "betaMatrix",
                                                      "excitationMatrix")){
                  x@bgNMF <- value
              }else{
                  stop("'bgNMF' value must be a list contain 5 elements.")
              }
              x
          }
)

setGeneric("setCluster", function(x,...) standardGeneric("setCluster"))
setMethod("setCluster", "EpiClusterResult",
          function(x,value)
          {
              if(class(value)=="factor"){
                  x@Cluster <- value
              }else{ 
                  stop("'Cluster' value must be a factor vector.")
              }
              x
                  
          }
)


setGeneric("setbetaRPMM", function(x,...) standardGeneric("setbetaRPMM"))
setMethod("setbetaRPMM", "EpiClusterResult",
          function(x,value)
          {
              if(class(value)=="blcTree"){
                  x@betaRPMM <- value
              }else{ 
                  stop("'betaRPMM' value must be a blcTree object.")
              }
              x
          }
)
