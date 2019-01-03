library(gtools)
library(RColorBrewer)
mth_g <- read.table("meth_nams.bed", sep="\t", header=FALSE)

mth_g <- mth_g[mth_g$V4 == "gene",]
mg <- (mth_g[,7:11] * 150)/(mth_g[,3] - mth_g[,2])

mg <- apply(mg, 2, FUN=function(x){(x/sum(x)) * 1000000})

nms <- c("V->M", "M->M", "V->V", "M->V", "Input")

colnames(mg) <- nms

mg <- t(apply(mg, 1, FUN=function(x){x/max(x)}))
mg  <- t(apply(mg, 1, FUN=function(x){x - gm_mean(x)}))

me_res$pvalue[is.na(me_res$pvalue)] <- 1
me_resO$pvalue[is.na(me_resO$pvalue)] <- 1
me_resD$pvalue[is.na(me_resD$pvalue)] <- 1

mgX <- mg[me_res$pvalue > 0.05 & me_resO$pvalue > 0.05 & me_resD$pvalue > 0.05,]

mgX <- mg[sample(1:nrow(mg), nrow(mg)/2),]

mgX[is.na(mgX)] <- 0

source("clusters/heatmap.3.R")
cols <- colorRampPalette(brewer.pal(9, "PuBuGn"))(256)

pdf("Meth_heat_fix.pdf", width=9, height=9)
par(cex=0.5)
heatmap.3(mgX, trace="none", col=cols, main="MEDIP Gene-Normalised Read Count", cex.main=0.9)
dev.off()

