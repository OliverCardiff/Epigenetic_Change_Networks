PlotFitQQ <- function(ot, otBig, cexp=1, xlb="", ylb="")
{
	plot(ot[[1]], ot[[4]], col=rgb(0,0,0,0), xlab=xlb, ylab=ylb)
	
	hrng <- ot[[9]] - ot[[7]]
	lrng <- ot[[7]] - ot[[8]]
	
	alph <- ot[[7]]/max(ot[[7]])
	
	arng <- ot[[7]] - ot[[12]]
	
	alw <- pmax(arng - lrng, 0)/lrng
	
	ahi <- abs(pmin(hrng + arng, 0))/hrng
	
	ah <- max(ahi)
	al <- max(alw)
	
	ahi <- ahi/max(ahi)
	alw <- alw/max(alw)
	
	for(i in 1:nrow(ot[[4]]))
	{
		for(j in 1:nrow(ot[[4]]))
		{
			b <- 1 + (2 * (i-1))
			b2 <- 1 + (2 * (j-1))
			inc1 <- 2
			inc2 <- 2

			if(i == nrow(ot[[4]]))
			{
				inc1 <- 1
			}
			if(j == nrow(ot[[4]]))
			{
				inc2 <- 1
			}
			
			al <- (ahi[i,j] + alw[i,j]) * 0.8
			
			cc <- rgb(alw[i,j], ahi[i,j], 0.7, 0.2 + al)
			
			xs <- c(otBig[[1]][b,b2], otBig[[1]][b+inc1,b2], otBig[[1]][b+inc1,b2], otBig[[1]][b,b2])
			ys <- c(otBig[[4]][b,b2], otBig[[4]][b,b2], otBig[[4]][b,b2+inc2], otBig[[4]][b,b2+inc2])
			polygon(xs, ys, border=rgb(0,0,0,0), col=cc)
		}
	}
}

PlotStack <- function(ot, cexp=1, xl="")
{
	div <- ncol(ot[[7]])
	xs <- c(0.5/div, (seq(1/div, 9/div, 1/div) + 0.5/div))
	par(mfrow=c(div+1,1), mar=c(0,6,0.3,0), cex=cexp)
	
	cmin <- min(c(ot[[12]], ot[[7]], ot[[8]], ot[[9]]))
	cmax <- max(c(ot[[12]], ot[[7]], ot[[8]], ot[[9]]))

	ctxt <- cmin + ((cmax - cmin) * 0.9)
	
	for(i in 1:ncol(ot[[9]]))
	{
		FC <- signif(mean(ot[[4]][,i]), 3)
		plot(xs, ot[[7]][i,], type="l", ylim=c(cmin,cmax), xaxt="n", yaxt="n", ylab=paste("Q", xs[i], "\nGrid Counts",sep=""), bty="n")
		axis(2, cex=0.7*cexp)
		grid(col="black")
		polygon(c(xs, rev(xs)), c(ot[[8]][i,], rev(ot[[9]][i,])), border=rgb(0,0,0,0), col=rgb(0.1,0.7,0.9))
		points(xs, ot[[7]][i,], type="l")
		points(xs, ot[[12]][i,], type="l", col="red", lwd=2)
		text(median(xs), ctxt, paste("IQ Mean: ", FC, sep=""), cex=0.8)
	}
	par(mar=c(5,6,0,0), cex=0.8*cexp)
	rs <- signif(rowMeans(ot[[1]]), 3)
	plot(xs, ot[[7]][i,], col=rgb(0,0,0,0),type="l", xaxt="n", yaxt="n", ylab="", ylim=c(cmin,cmax), bty="n", xlab=xl)
	axis(1, at=xs, labels=as.character(rs), cex=0.9*cexp)
}
Plot2DQQ <- function(ot, ymx=NA, ymin=NA, xmx=NA, xmin=NA, inc=0.4, xlb="", ylb="")
{
	for(i in 1:length(ot))
	{
		ot[[i]][!is.finite(ot[[i]])] <- 0
	}
	
	if(is.na(ymx))
	{
		ymx = max(ot[[4]])
	}
	if(is.na(ymin))
	{
		ymin = min(ot[[4]])
	}
	if(is.na(xmx))
	{
		xmx = max(ot[[1]])
	}
	if(is.na(xmin))
	{
		xmin = min(ot[[1]])
	}
	
	cmin <- min(c(ot[[12]], ot[[7]], ot[[8]], ot[[9]]))
	cmax <- max(c(ot[[12]], ot[[7]], ot[[8]], ot[[9]]))

	dcnt <- pmax((ot[[7]] - cmin)/(cmax-cmin),0)
	dcl <- pmax((ot[[8]] - cmin)/(cmax-cmin),0)
	dch <- pmax((ot[[9]] - cmin)/(cmax-cmin),0)

	tcnt <- pmax((ot[[12]] - cmin)/(cmax-cmin),0)
	
	#print(tcnt)
	#print(dch)
	#print(dcl)
	#print(dcnt)


	plot(as.vector(ot[[1]]), as.vector(ot[[4]]), pch=".", col=rgb(0,0,0,0), xlim=c(xmin, xmx), ylim=c(ymin, ymx), ylab=ylb, xlab=xlb)
	symbols(c(-10000,-10000,as.vector(ot[[1]])), c(-10000,-10000,as.vector(ot[[4]])), lwd=1, bg=rgb(0,0,0,0.3), fg=rgb(1,1,1,0.0), circles=c(0,1,as.vector(dch)), inches=inc, add=TRUE)
	symbols(c(-10000,-10000,as.vector(ot[[1]])), c(-10000,-10000,as.vector(ot[[4]])), lwd=1, circles=c(0,1,as.vector(dcnt)), inches=inc, add=TRUE)
	#symbols(c(-10000,-10000,as.vector(ot[[1]])), c(-10000,-10000,as.vector(ot[[4]])), lwd=1, bg=rgb(1,1,1,1), fg=rgb(1,1,1,0.2), circles=c(0,1,as.vector(dcl)), inches=inc, add=TRUE)
	segments(x0=ot[[2]], x1=ot[[3]], y0=ot[[4]], y1=ot[[4]], lwd=2)
	segments(x0=ot[[1]], x1=ot[[1]], y0=ot[[5]], y1=ot[[6]], lwd=2)
	segments(x0=ot[[1]], x1=ot[[10]], y0=ot[[4]], y1=ot[[11]], col="red")

	symbols(c(-10000,-10000,as.vector(ot[[1]])), c(-10000,-10000,as.vector(ot[[4]])), lwd=1, bg=rgb(1,0,0,0.5),fg="red", circles=c(0,1,as.vector(tcnt)), inches=inc, add=TRUE)

	abline(v=ot[[13]], lty=2, col=rgb(0,0,0,0.2))
	abline(h=ot[[14]], lty=2, col=rgb(0,0,0,0.2))

}

Plot2DQQGrid <- function(ot, ymx=NA, ymin=NA, xmx=NA, xmin=NA, inc=0.4, ylb="Expression Quantiles", xlb="Methylation Quantiles")
{
	if(is.na(ymx))
	{
		ymx = max(ot[[4]])
	}
	if(is.na(ymin))
	{
		ymin = min(ot[[4]])
	}
	if(is.na(xmx))
	{
		xmx = max(ot[[1]])
	}
	if(is.na(xmin))
	{
		xmin = min(ot[[1]])
	}
	cmin <- min(c(ot[[12]], ot[[7]], ot[[8]], ot[[9]]))
	cmax <- max(c(ot[[12]], ot[[7]], ot[[8]], ot[[9]]))

	dcnt <- pmax((ot[[7]] - cmin)/(cmax-cmin),0)
	dcl <- pmax((ot[[8]] - cmin)/(cmax-cmin),0)
	dch <- pmax((ot[[9]] - cmin)/(cmax-cmin),0)

	tcnt <- pmax((ot[[12]] - cmin)/(cmax-cmin),0)
	
	#print(tcnt)
	#print(dch)
	#print(dcl)
	#print(dcnt)
	
	ys <- rep(1, nrow(ot[[4]]))
	
	for(i in 2:nrow(ot[[4]]))
	{
		ys <- c(ys, rep(i, nrow(ot[[4]])))
	}
	
	xs <- rep(c(1:nrow(ot[[4]])), nrow(ot[[4]]))


	plot(xs, ys, pch=".", col=rgb(0,0,0,0), ylab=ylb, xlab=xlb)
	symbols(c(-10000,-10000,xs), c(-10000,-10000,ys), lwd=1, bg=rgb(0,0,0,0.3), fg=rgb(1,1,1,0.0), circles=c(0,1,as.vector(dch)), inches=inc, add=TRUE)
	symbols(c(-10000,-10000,xs), c(-10000,-10000,ys), lwd=1, circles=c(0,1,as.vector(dcnt)), inches=inc, add=TRUE)
	#symbols(c(-10000,-10000,xs), c(-10000,-10000,ys), lwd=1, bg=rgb(1,1,1,1), fg=rgb(1,1,1,0.2), circles=c(0,1,as.vector(dcl)), inches=inc, add=TRUE)
	#segments(x0=ot[[2]], x1=ot[[3]], y0=ot[[4]], y1=ot[[4]], lwd=2)
	#segments(x0=ot[[1]], x1=ot[[1]], y0=ot[[5]], y1=ot[[6]], lwd=2)
	#segments(x0=ot[[1]], x1=ot[[10]], y0=ot[[4]], y1=ot[[11]], col="red")

	symbols(c(-10000,-10000,xs), c(-10000,-10000,ys), lwd=1, bg=rgb(1,0,0,0.5),fg="red", circles=c(0,1,as.vector(tcnt)), inches=inc, add=TRUE)

	#abline(v=ot[[13]], lty=2, col=rgb(0,0,0,0.2))
	#abline(h=ot[[14]], lty=2, col=rgb(0,0,0,0.2))

}

Run2DQQ <- function(d1, d2, div, boots, trim=0)
{
	outs <- list()

	if(trim == 1)
	{
		d2 <- d2[d1 > min(d1)  & d1 < max(d1)]
		d1 <- d1[d1 > min(d1)  & d1 < max(d1)]
	}
	if(trim == 2)
	{
		d1 <- d1[d2 > min(d2)  & d2 < max(d2)]
		d2 <- d2[d2 > min(d2)  & d2 < max(d2)]
	}
	sndat <- rep(0, div*div*boots)
	CNTmat <- array(sndat, c(div, div, boots))
	xMNmat <- array(sndat, c(div, div, boots))
	yMNmat <- array(sndat, c(div, div, boots))
	
	fmatX <- matrix(0,div,div)
	fmatY <- matrix(0,div,div)
	fmatC <- matrix(0,div,div)

	qnts <- seq(1/div, (div-1)/div, 1/div)

	qn1 <- quantile(d1, qnts)
	qn2 <- quantile(d2, qnts)

	mids1 <- qn1
	mids2 <- qn2

	for(i in 1:(length(qnts)+1))
	{
		l1 <- 0
		h1 <- max(d1)
		l2 <- 0
		h2 <- max(d2)

		if(i > 1)
		{
			l1 <- qn1[i-1]
			l2 <- qn2[i-1]
		}
		if(i <= length(qnts))
		{
			h1 <- qn1[i]
			h2 <- qn2[i]
		}

		mids1[i] <- l1 + ((h1-l1)/2)
		mids2[i] <- l2 + ((h2-l2)/2)
	}

	for(j in 1:div)
	{
		jLow <- min(d1)
		jHigh <- max(d1)
		if(j > 1)
		{
			jLow <- qn1[j - 1]
		}
		if(j < div)
		{
			jHigh <- qn1[j]
		}
		for(k in 1:div)
		{
			kLow <- min(d2)
			kHigh <- max(d2)
			if(k > 1)
			{
				kLow <- qn2[k - 1]
			}
			if(k < div)
			{
				kHigh <- qn2[k]
			}
			
			mnx <- mids1[j]
			mny <- mids2[k]
			
			for(i in 1:boots)
			{
				s1 <- sample(d1)
				s2 <- sample(d2)

				slcts <- s1 >= jLow & s1 <= jHigh & s2 >= kLow & s2 <= kHigh
				tx <- s1[slcts]
				ty <- s2[slcts]				
			
				if(length(tx) > 0)
				{
					mnx <- mean(tx)
				}
				if(length(ty) > 0)
				{
					mny <- mean(ty)
				}

				xMNmat[j,k,i] <- mnx
				yMNmat[j,k,i] <- mny
				CNTmat[j,k,i] <- length(tx)
			}
		
			mnx <- mids1[j]
			mny <- mids2[k]
			
			slcts <- d1 >= jLow & d1 <= jHigh & d2 >= kLow & d2 <= kHigh
			
			tx <- d1[slcts]
			ty <- d2[slcts]
			
			if(length(tx) > 0)
			{
				mnx <- mean(tx)
			}
			if(length(ty) > 0)
			{
				mny <- mean(ty)
			}

			fmatX[j,k] <- mnx
			fmatY[j,k] <- mny

			fmatC[j,k] <- length(tx)
		}
	}

	for(i in 1:9)
	{
		outs[[i]] <- matrix(0,div,div)
	}

	for(j in 1:div)
	{
		for(k in 1:div)
		{
			vofi <- xMNmat[j,k,]

			outs[[1]][j,k] <- mean(vofi)
			outs[[2]][j,k] <- quantile(vofi, .05)
			outs[[3]][j,k] <- quantile(vofi, .95)

			vofi <- yMNmat[j,k,]

			outs[[4]][j,k] <- mean(vofi)
			outs[[5]][j,k] <- quantile(vofi, .05)
			outs[[6]][j,k] <- quantile(vofi, .95)

			vofi <- CNTmat[j,k,]

			outs[[7]][j,k] <- mean(vofi)
			outs[[8]][j,k] <- quantile(vofi, .05)
			outs[[9]][j,k] <- quantile(vofi, .95)
		}
	}

	outs[[10]] <- fmatX
	outs[[11]] <- fmatY
	outs[[12]] <- fmatC
	outs[[13]] <- mids1
	outs[[14]] <- mids2
	
	outs
}