EpiDraw <- function(EpiCluster.output,
                    PhenoTypes=NULL,
                    maxlevel=NULL,
                    threshold=10,
                    legend=TRUE,
                    samplecluster=TRUE,
                    maintitle="",
                    colors=c("yellow","navy"),
                    clusterplot=TRUE,
                    ComponentMatrix=TRUE)
{
    if(class(EpiCluster.output)!="EpiClusterResult")
    {
        print(class(EpiCluster.output),
              "is not correct EpiCluster Result Format.")
        return(FALSE)
    }

    Basis <- EpiCluster.output@bgNMF$alphaMatrix / 
    (EpiCluster.output@bgNMF$alphaMatrix+EpiCluster.output@bgNMF$betaMatrix)
    Exitation <- EpiCluster.output@bgNMF$exitationMatrix
    Cluster <- EpiCluster.output@Cluster
    SampleColfunc <- colorRampPalette(c("deepskyblue","firebrick1"))


    if(is.null(maxlevel))
    {
        Q <- as.character(Cluster)
        maxlevel <- "Inf"
    }else
        Q <- substr(Cluster,1,maxlevel)

    Q2 <- Q[Q %in% names(which(table(Q) >= threshold))]
    if(length(Q2)==0)
    {
        cat("Sorry There is not cluster(based on maxlevel",maxlevel,
            ") contain ",threshold,"or more Samples. You may want to
            set lower the threshold or increase the maxlevel.")
    }
    S2 <- c(1:dim(Basis)[1])[Q %in% names(which(table(Q) >= threshold))]



    if(ComponentMatrix==TRUE)
    {
        par(mar=c(6,3,3,4))
        hr <- hclust(as.dist(1-cor(Basis[S2,], 
                    method="pearson")),method="complete")
        colfunc <- colorRampPalette(colors)
        image(Basis[S2,][order(Q2),hr$order],col=colfunc(50),
              xlab="",ylab="",axes=FALSE)
        box()
        mtext(side=2,paste(as.character(ncol(Basis[S2,])),
            "Estimated Components"),line=1,cex=1,font=2)
        axis(2,at=seq(0,1,length.out=ncol(Basis[S2,])),labels= hr$order,
             mgp=c(1,0.1,0),cex.axis=0.8,las=2,tick=FALSE,font=2)
        mtext(side=3,maintitle,line=1,cex=1.5,font=2)
        for(i in cumsum(table(Q2)))
            abline(v=i/length(Q2),lwd=2)

        if(is.null(PhenoTypes)!=TRUE &
           mode(PhenoTypes)!="list" &
           class(PhenoTypes)!="numeric")
        {
            mai <- par("mai")
            fin <- par("fin")
            x.legend.fig <- c(mai[2]/fin[1],1-(mai[4]/fin[1]))
            y.legend.fig <- c(mai[1]/fin[2], 1 - (mai[3]/fin[2]))
            x.legend.plt <- x.legend.fig
            y.legend.plt <- c(y.legend.fig[2],y.legend.fig[2] + 
                              (0.02 * (y.legend.fig[2] - y.legend.fig[1])))
            par(new = TRUE, pty = "m", plt = c(x.legend.plt, y.legend.plt))

            Pheno <- PhenoTypes[S2]
            index <- c(1:length(table(Pheno)))
            names(index) <- names(table(Pheno))
            q <- index[Pheno]
            q <- q[order(Q2)]
            image(matrix(q,ncol=1),col = SampleColfunc(length(table(Q2))), 
                  xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        }
    

    if(legend==TRUE)
    {
        zlim=range(Basis[S2,])
        mai <- par("mai")
        fin <- par("fin")
        x.legend.fig <- c(1 - (mai[4]/fin[1]), 1)
        y.legend.fig <- c(mai[1]/fin[2], 1 - (mai[3]/fin[2]))
        x.legend.plt <- c(x.legend.fig[1] + (0.09 * (x.legend.fig[2] -
            x.legend.fig[1])), x.legend.fig[2] - 
            (0.6 * (x.legend.fig[2] - x.legend.fig[1])))
        y.legend.plt <- y.legend.fig
        cut.pts <- seq(zlim[1], zlim[2], length = length(colfunc(50)) + 1)
        z <- (cut.pts[1:length(colfunc(50))] + 
              cut.pts[2:(length(colfunc(50)) + 1)])/2
        z <- c(zlim[1],z,zlim[2])
        y.legend.plt[1] <- 0.5
        par(new = TRUE, pty = "m", plt = c(x.legend.plt, y.legend.plt))
        image(x = 1, y = z, z = matrix(z, nrow = 1, 
            ncol = length(colfunc(52))),col = colfunc(52), xlab = "", 
            ylab = "", xaxt = "n", yaxt = "n")
        axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis = 0.6, tcl = -0.1)
    }

    if(samplecluster == TRUE)
    {
        ncluster <- length(table(Q2))
        mai <- par("mai")
        fin <- par("fin")
        x.legend.fig <- c(mai[2]/fin[1],1-(mai[4]/fin[1]))
        y.legend.fig <- c(mai[1]/fin[2], 1 - (mai[3]/fin[2]))
        x.legend.plt <- x.legend.fig
        y.legend.plt <- c(y.legend.fig[1] - (0.02 * (y.legend.fig[2] -
                        y.legend.fig[1])), y.legend.fig[1])
        par(new = TRUE, pty = "m", plt = c(x.legend.plt, y.legend.plt))
        q <- rep(1:ncluster,table(Q2))
        image(matrix(q,ncol=1),col = rainbow(ncluster), 
              xlab = "", ylab = "", xaxt = "n", yaxt = "n") 
    }

    {
        if(is.null(PhenoTypes)!=TRUE & 
           mode(PhenoTypes)!="list" & 
           class(PhenoTypes)!="numeric")
        {
            par(new = TRUE, mar=c(1,3,1,4))
            plot(1:10,type="n",axes=FALSE)
            npheno <- length(table(PhenoTypes[S2]))
            legend("bottomleft",legend = names(table(PhenoTypes[S2])),
                   fill=SampleColfunc(npheno),ncol=npheno/2,cex = 0.75)
        }
        par(new = TRUE, mar=c(1,3,1,4))
        plot(1:10,type="n",axes=FALSE)
        legend("bottomright",legend = names(table(Q2)),
               fill=rainbow(length(table(Q2))),ncol=ncluster/2,cex = 0.75)
    }

    }

    if(clusterplot == TRUE & mode(PhenoTypes) == "list")
    {
        if("" %in% names(PhenoTypes))
        {
            names(PhenoTypes) <- paste("Factor_",1:length(PhenoTypes),sep="")
        }
        for(element in 1:length(PhenoTypes))
        {
            covariate <- PhenoTypes[S2,element]
            cov.name <- names(PhenoTypes)[element]
            if(class(covariate) == "numeric")
            {
                par(mar=c(5,5,5,5))
                ncluster <- length(table(Q2))
                boxplot(covariate ~ Q2,col=rainbow(ncluster),cex.axis=1,
                        ann=FALSE,boxlwd=1,ylab=cov.name)
                mtext(side=2,cov.name,line=1.7,cex=0.7,font=2)
            }
            else
            {
                par(mar=c(5,5,5,5))
                pheno.cluster <- table(covariate,Q2)
                barplot(pheno.cluster,col = 
                        SampleColfunc(length(table(covariate))),
                        cex.axis=1,cex.names=1,ylab=cov.name)
                legend("topright",legend = rownames(pheno.cluster),
                    fill=SampleColfunc(length(table(covariate))),cex=0.75)
            }
        }
    }

    if(clusterplot == TRUE & mode(PhenoTypes) != "list")
    {
        if(is.null(PhenoTypes)==TRUE)
        {
            print("Please provide PhenoType information.")
        }
        else if(class(PhenoTypes)=="numeric")
        {
            par(mar=c(5,5,5,5))
            ncluster <- length(table(Q2))
            boxplot(PhenoTypes[S2] ~ Q2,col=rainbow(ncluster),
                    cex.axis=1,ann=FALSE,boxlwd=1)
            mtext(side=2,"PhenoType",line=1.7,cex=0.7,font=2)
        }
        else
        {
            par(mar=c(5,5,5,5))
            pheno.cluster <- table(PhenoTypes[S2],Q2)
            barplot(pheno.cluster,
                    col=SampleColfunc(length(table(PhenoTypes[S2]))),
                    ylab="PhenoType",cex.axis=1,cex.names=1)
            legend("topright",legend = rownames(pheno.cluster),
                fill=SampleColfunc(length(table(PhenoTypes[S2]))),cex =0.75)
        }
    }
}
