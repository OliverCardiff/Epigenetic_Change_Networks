library(gtools)
library(RColorBrewer)

trpl <- read.table("gexp_genome/transp_uniq.txt", sep="\t", header=FALSE)
mt <- floor(trpl[,18:29])

nm2 <- c("V->M.1", "M->M.1", "V->M.2", "M->M.2", "V->V.1", "M->V.1", "V->V.2", "M->V.2", "V->V.3", "M->V.3", "M->M.3", "V->M.3")

colnames(mt) <- nm2

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

arthMean <- apply(mt, 1, FUN=function(x){gm_mean(x)})
rats <- apply(mt, 2, FUN=function(x){x/arthMean})

facts <- rep(0, ncol(rats))

for(i in 1:ncol(rats))
{
	tst <- rats[,i]
	tst <- tst[ tst > 0 ]
	facts[i] <- median(tst)

	mt[,i] <- mt[,i] / facts[i]

}

mt <- t(apply(mt, 1, collaPS))
mt <- t(apply(mt, 1, FUN=function(x){x/max(x)}))
mt <- t(apply(mt, 1, FUN=function(x){x - gm_mean(x)}))

res$pvalue[is.na(res$pvalue)] <- 1
resO$pvalue[is.na(resO$pvalue)] <- 1
resD$pvalue[is.na(resD$pvalue)] <- 1
resM$pvalue[is.na(resM$pvalue)] <- 1
resV$pvalue[is.na(resV$pvalue)] <- 1

mtx <- mt[res$pvalue < 0.05 | resO$pvalue < 0.05 | resD$pvalue < 0.05 | resV$pvalue < 0.05 | resM$pvalue < 0.05,]

mtx <- mt[sample(1:nrow(mt), nrow(mt)/3),]

mtx[is.na(mtx)] <- 0

#colnames(mtx) <- nm2

source("clusters/heatmap.3.R")
cols <- colorRampPalette(brewer.pal(9, "PuBuGn"))(256)

pdf("RNA_heat.pdf", width=9, height=9)
par(cex=0.5)
heatmap.3(mtx, trace="none", col=cols, main="RNA Gene-Normalised Read Count", cex.main=0.9)
dev.off()

collaPS <- function(x)
{
	ref <- seq(1,length(x),1)
	y <- x[rev(order(x))]
	ref <- ref[rev(order(x))]
	y[1] <- y[2]
	x <- y[rev(order(ref))]
	x
}