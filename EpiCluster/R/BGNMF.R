### BGNMF function

BGNMF <- function(Data,nIter=60,K=NULL,EE=0.0001)
{
    if(class(Data) == "ExpressionSet"){
        X <- t(exprs(Data))
    }else{
        X <- t(Data)
    }

	N <- nIter
	if(is.null(K))
    {
        tmp.o <- Data-rowMeans(Data)
        rmt.o <- EstDimRMT(tmp.o)
		K <- rmt.o$dim+1
    }
    print(paste("Iteration times: ",as.character(N)))
    print(paste("Components number: ",as.character(K)))
    print(paste("ee parameter: ",as.character(EE)))

	Prow <- nrow(X)
	Tcol <- ncol(X)

	X[which(X < 1e-04)] <- 1e-04
	X[which(X > 0.9999)] <- 0.9999


	M <- matrix(1, nrow = Prow, ncol = Tcol)


	ee <- EE
	Zeta0 <- 1
	Rho0 <- ee

	Alpha0 <- ee
	Beta0 <- ee
	Mu0 <- 1
	Nu0 <- 1

	aa <- rgamma(Prow * K, Mu0, rate = Alpha0)
	bb <- rgamma(Prow * K, Nu0, rate = Beta0)
	zz <- rgamma(K * Tcol, Rho0, rate = Zeta0)

	A <- matrix(aa, Prow, K)
	B <- matrix(bb, Prow, K)
	z <- matrix(zz, K, Tcol)

	A[which(A < ee^2)] <- ee^2
	B[which(B < ee^2)] <- ee^2
	z[which(z < ee^2)] <- ee^2

	Atilde <- A
	Btilde <- B
	ztilde <- z

    message("Start Iteration:")
    message("--------------------------")

	ELogPOut <- matrix(0, nrow = N, ncol = 1)
	scale = 1

	for (i in 1:N) {                                                                         
	    Alpha <- Alpha0 - (scale * M * log(X)) %*% t(z)
	    Mu <- (scale * M * (psigamma((A + Btilde) %*% ztilde) -
                            psigamma(A %*% ztilde))) %*% t(ztilde) * A + Mu0
	    
	    Beta <- Beta0 - (scale * M * log(1 - X)) %*% t(z)
	    Nu <- (scale * M * (psigamma((Atilde + B) %*% ztilde) -
                            psigamma(B %*% ztilde))) %*% t(ztilde) * B + Nu0
	    
	    Zeta <- Zeta0 - scale *(t(A) %*% (M * log(X)) + t(B) %*% (M * log(1 - X)))
	    Rho <- Rho0 + (scale * t(Atilde)) %*% (M * (psigamma((Atilde + Btilde) %*% 
                z) - psigamma(Atilde %*% z))) * z + (scale * t(Btilde)) %*% 
	        (M * (psigamma((Atilde + Btilde) %*% z) - psigamma(Btilde %*% z)))* z
	    
	    A <- Mu/Alpha
	    B <- Nu/Beta
	    z <- Rho/Zeta
	    
	    Atilde <- A
	    Btilde <- B
	    ztilde <- z
	    
	    tempmat <- (A %*% z - 1) * log(X) * M + (B %*% z - 1) *
        log(1 - X) * M - lbeta(A %*% z, B %*% z) * M
	    ELogPOut[i, 1] <- sum(apply(tempmat, 2, sum))
	    
	    print(i)
	    
	    X_hat <- (A %*% z)/((A + B) %*% z)
	} 
    bgNMF.output <- list(ee=EE,
                         X_reconstruct=X_hat,
                         alphaMatrix=A,
                         betaMatrix=B,
                         excitationMatrix=z)
	return(bgNMF.output)
}
