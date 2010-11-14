getTdtSnps <- function(x){
	chrom <- unique(x[, "CHR"])
	dirs <- getDirs(CIDR=FALSE, MIN.SAMPLES=3)
	aa <- list()
	bb <- list()
	gg <- list()

	for(i in seq(along=chrom)){
		CHR <- chrom[i]
		xx <- x[x[, "CHR"] == CHR, ]
		id <- xx[, "SNP"]
		a <- list()
		b <- list()
		gt <- list()
		for(p in seq(along=dirs)){
			message("Loading:", file.path(dirs[p], paste("crlmmSetList_", CHR, ".rda", sep="")))
			load(file.path(dirs[p], paste("crlmmSetList_", CHR, ".rda", sep="")))
			callSetPlus <- as(crlmmSetList, "SnpCallSetPlus")
			if(!"familyId" %in% varLabels(callSetPlus)){
				phenoData(callSetPlus) <- addSampleSheet(callSetPlus)
			}
			callSetPlus <- callSetPlus[, -grep("CIDR", callSetPlus$Sample.Name)]
			a[[p]] <- A(callSetPlus)[match(id, featureNames(callSetPlus)), ]
			b[[p]] <- B(callSetPlus)[match(id, featureNames(callSetPlus)), ]
			gt[[p]] <- calls(callSetPlus)[match(id, featureNames(callSetPlus)), ]
		}
		aa[[i]] <- do.call("cbind", a)
		bb[[i]] <- do.call("cbind", b)
		gg[[i]] <- do.call("cbind", gt)
		rm(a, b, gt)
	}
	results <- list(aa=aa, bb=bb, gg=gg)
	save(results, file="results.rda")
}

tdtMatrix <- function(){
	tmp <- matrix(NA, nrow=3^3, ncol=7)
	colnames(tmp) <- c("father", "mother", "offspring", "transmittedB", "untransmittedB", "transmittedA", "untransmittedA")
	counter <- 1
	for(i in c("AA", "AB", "BB")){##genotype of father
		for(j in c("AA", "AB", "BB")){ ## genotype of mother
			for(k in c("AA", "AB", "BB")){ ##genotype of offspring
				tmp[counter, 1:3] <- c(i, j, k)
				counter <- counter+1
			}
		}
	}
	tmp[, "transmittedB"] <-   c(0, 0, 0, 0, 1, 1, 0, 1, 2, 0, 1,1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2)
	tmp[, "untransmittedB"] <- c(0, 0, 0, 1, 0, 0, 2, 1, 0, 1, 0, 0, 2, 1, 0, 3, 2, 1, 2, 1, 0, 3, 2, 1, 4, 3, 2)
	tmp[, c(6, 7)] <-   matrix(c(2, 2,
				     1, 3,
				     0, 4,
				     2, 1,
				     1, 1,
				     0, 3,
				     2, 0,
				     1, 1,
				     0, 2,
				     2, 1,
				     1, 2,
				     0, 3,
				     2, 0,
				     1, 1,
				     0, 2,
				     1, 0,
				     1, 0,
				     0, 1,
				     2, 0,
				     1, 1,
				     0, 2,
				     1, 0,
				     1, 0,
				     0, 1,
				     0, 0,
				     0, 0,
				     0, 0), byrow=TRUE, ncol=2)
	return(tmp)
}


computeTdt <- function(crlmmCalls){
	##match(colnames(crlmmCalls), substr(samplesheet[, "Sample.Name"], 1, 8))
	col.index <- split(1:ncol(crlmmCalls), substr(colnames(crlmmCalls), 1, 5))
	tdtM <- tdtMatrix()
	TU <- vector("list", length(col.index))
	for(m in seq(along=col.index)){
		if(length(col.index[[m]]) < 3) next()  ##next family
		gts <- crlmmCalls[, col.index[[m]]]
		colnames(gts) <- who3(colnames(gts))
		gts <- gts[, match(c("father", "mother", "offspring"), colnames(gts))]
		##possibilities
		tmp <- matrix(NA, nrow(crlmmCalls), 2)
		colnames(tmp) <- c("T", "U")
		counter <- 1
		for(i in 1:3){##genotype of father
			for(j in 1:3){ ## genotype of mother
				for(k in 1:3){ ##genotype of offspring
					index <- which(gts[, 1] == i & gts[, 2] == j & gts[, 3] == k)
					tmp[index, "T"] <- as.integer(tdtM[counter, "transmittedB"]) + as.integer(tdtM[counter, "transmittedA"])
					tmp[index, "U"] <- as.integer(tdtM[counter, "untransmittedB"]) + as.integer(tdtM[counter, "untransmittedA"])
					counter <- counter+1
				}
			}
		}
		TU[[m]] <- tmp
	}
	matFun <- function(lis, FUN){
		lis <- lis[!sapply(lis, is.null)]
		if(!is.list(lis) || !all(sapply(lis, is.matrix)))
			stop("'lis' must be a list containing 2-dimensional arrays")
		dims <- sapply(lis, dim)
		n <- dims[1, 1]
		p <- dims[2, 1]
		if(!all(n == dims[1, ]) || !all(p == dims[2, ]))
			stop("the matrices must have the same dimensions")
		mat <- matrix(unlist(lis), n * p, length(lis))
		matrix(apply(mat, 1, FUN), n, p)
	}
	err <- tryCatch(tmp <- matFun(TU, "sum"), error=function(e) NULL)
	if(!is.null(err)){
		TU <- tmp
	} ##otherwise, keep TU  in list form
	return(TU)
}

my.qq.func=function(p,f,pval=TRUE,gc=FALSE,tn,hc=0.99,hm=100,...){
  if(missing(f)) f=rep(0,length(p))
  wh=which(is.na(p)|is.na(f))
  if(length(wh)>0){
    p=p[-wh]
    f=f[-wh]
  }
  wh=order(p)
  p=p[wh]
  f=f[wh]
  if(gc){
    z=qchisq(p,1,lower=F)
    gcp=median(z)/qchisq(0.5,1)
    cat("The genomic control parameter is ",gcp,"\n")
    z=z/gcp
    pchisq(z,1,lower=F)
  }
  n=length(p)
  x1=1:n
  x2=n+1-x1
  x=x1/(n+1)
  up=qbeta(0.975,x1,x2)
  lo=qbeta(0.025,x1,x2)
  if(pval==TRUE){
    uu=up
    up=-log10(lo)
    lo=-log10(uu)
    x=-log10(x)
    z=-log10(p)
  }
  else{
    up=qchisq(up,1)
    lo=qchisq(lo,1)
    x=qchisq(x,1)
    z=qchisq(p,1,lower=F)
  }
  z=rev(z)
  f=rev(f)
  x=rev(x)
  up=rev(up)
  lo=rev(lo)
  tt=NULL
  if(!missing(tn)){
    mx=10^tn
    tt=list(tn=mx,k=rev(x)[mx])
  }
  if(hm>1){
    nc=round(n*hc)
    wh=seq(1,nc,hm)
    wh=c(wh,(nc+1):n)
    z=z[wh]
    f=f[wh]
    x=x[wh]
    up=up[wh]
    lo=lo[wh]
  }
  return(list(z=z,f=f,x=x,up=up,lo=lo,tt=tt))
}

my.qq.plot=function(zz,cut,rmx=0,mt="",mt.cex=1,mt.line=NA,mgp=c(3,1,0),tn.cex=1,plab=T,xr,yr,...){
  if(!missing(cut)) zz$z[zz$z>cut]=cut
  if(missing(xr)) xr=c(0,1.02*max(zz$x))
  if(missing(yr)) yr=c(0,1.02*max(c(zz$up,zz$z,rmx)))
  if(plab){
    xl=expression(paste("expected  ",-log[10]," (p-value)",sep=""))
    yl=expression(paste("observed  ",-log[10]," (p-value)",sep=""))
  }
  else{
    xl=""
    yl=""
  }
  par(las=1)
  print(xr)
  plot(range(zz$x),range(c(0,zz$up)),type="n",xlim=xr,ylim=yr,xlab=xl,ylab=yl,...)
  axis(2,0:ceiling(yr)[2],...)
  polygon(c(zz$x,rev(zz$x)),c(zz$lo,rev(zz$up)),col="lightgrey",border=F)
  lines(c(0,max(zz$x)),c(0,max(zz$x)))
  cls=c("blue","red")
  cls=cls[zz$f+1]
  points(zz$x,zz$z,pch=20,col=cls,cex=0.5)
  par(mgp=mgp)
  if(length(zz$tt)==2){
    options(scipen=7)
    axis(3,zz$tt$k,zz$tt$tn,cex.axis=tn.cex)
    options(scipen=0)
  }
  title(mt,cex.main=mt.cex,line=mt.line)
  par(mgp=c(3,1,0))
}


