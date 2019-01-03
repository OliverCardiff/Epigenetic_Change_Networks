library(gtools)
library(RColorBrewer)

nms <- c("1F", "1M", "3F", "3M", "6F", "6M", "8F", "8M", "10F", "10M", "12F", "12M")
nm2 <- c("V->V.1", "M->V.1", "V->V.2", "M->V.2", "V->V.3", "M->V.3", "V->M.1", "M->M.1", "V->M.2", "M->M.2", "V->M.3", "M->M.3")



mrna <- list()

hl <- read.table("Counts_novel_and_mirbase/high_list.txt", header=FALSE)

for(i in 1:length(nms))
{
	mrna[[i]] <- read.table(paste("Counts_novel_and_mirbase/", nms[i], "_finalcounts2b.txt", sep=""), sep="\t", header=FALSE)
}

c_mat <- merge(hl, mrna[[1]], by="V1")

accu <- c_mat$V2

for(i in 2:length(mrna))
{
	c_mat <- merge(hl, mrna[[i]], by="V1", all.x = TRUE )

	c_mat$V2[is.na(c_mat$V2)] <- 0

	accu <- cbind(accu, c_mat$V2)
}

colnames(accu) <- nm2

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

arthMean <- apply(accu, 1, FUN=function(x){gm_mean(x)})
rats <- apply(accu, 2, FUN=function(x){x/arthMean})

facts <- rep(0, ncol(rats))

for(i in 1:ncol(rats))
{
	tst <- rats[,i]
	tst <- tst[ tst > 0 ]
	facts[i] <- median(tst)

	accu[,i] <- accu[,i] / facts[i]

}

accu <- t(apply(accu, 1, collaPS))
accu <- t(apply(accu, 1, FUN=function(x){x/max(x)}))
accu <- t(apply(accu, 1, FUN=function(x){x - mean(x)}))

mi_res$pvalue[is.na(mi_res$pvalue)] <- 1
mi_resO$pvalue[is.na(mi_resO$pvalue)] <- 1
mi_resD$pvalue[is.na(mi_resD$pvalue)] <- 1
mi_resV$pvalue[is.na(mi_resV$pvalue)] <- 1
mi_resM$pvalue[is.na(mi_resM$pvalue)] <- 1

accuX <- accu[mi_res$pvalue < 0.05 | mi_resO$pvalue < 0.05 | mi_resD$pvalue < 0.05 | mi_resV$pvalue < 0.05 | mi_resM$pvalue < 0.05,]

accuX[is.na(accuX)] <- 0

source("clusters/heatmap.3.R")
cols <- colorRampPalette(brewer.pal(9, "PiYG"))(256)

pdf("miRNA_heat_fix.pdf", width=9, height=9)
par(cex=0.5)
heatmap.3(accuX, trace="none", col=cols, main="miRNA Expression", cex.main=0.9)
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
