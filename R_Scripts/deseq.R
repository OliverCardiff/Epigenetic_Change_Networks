library(DESeq2)
/
trpl <- read.table("gexp_genome/transp_uniq.txt", sep="\t", header=FALSE)
mt <- floor(trpl[,18:29])

row.names(mt) <- trpl$V1

origin <- c("V", "M", "V", "M", "V", "M", "V", "V", "M", "V", "M", "V")

destination <- c("M", "M", "M", "M", "V", "V","V", "V", "V", "V", "M", "M")

transp <- c("VM", "MM", "VM", "MM", "VV", "MV", "VV", "MV", "VV", "MV", "VM", "MM")

change <- c("C", "S", "C", "S", "S", "C", "S", "C", "S", "C", "C", "S")


VM <- c(1,3,11)
MM <- c(2, 4, 12)
MV <- c(6, 8, 10)
VV <- c(5, 7, 9)

coldata <- cbind(origin=(origin), destination=(destination), transp=transp, change=change)

cVdata <- coldata[c(VM, VV),]
cMdata <- coldata[c(MV,MM),]

mtV <- mt[,c(VM, VV)]
mtM <- mt[,c(MV,MM)]

ddsM <- DESeqDataSetFromMatrix(countData = mtM,
                              colData = cMdata,
                              design= ~ transp)

ddsM <- DESeq(ddsM)
resultsNames(ddsM)

ddsV <- DESeqDataSetFromMatrix(countData = mtV,
                              colData = cVdata,
                              design= ~ transp)

ddsV <- DESeq(ddsV)
resultsNames(ddsV)

dds <- DESeqDataSetFromMatrix(countData = mt,
                              colData = coldata,
                              design= ~ origin + destination + change)

dds <- DESeq(dds)
resultsNames(dds)


res <- lfcShrink(dds, coef="change_S_vs_C", type="normal")
resD <- lfcShrink(dds, coef="destination_V_vs_M", type="normal")
resO <- lfcShrink(dds, coef="origin_V_vs_M", type="normal")

resV <- lfcShrink(ddsV, coef="transp_VV_vs_VM", type="normal")
resM <- lfcShrink(ddsM, coef="transp_MV_vs_MM", type="normal")

pcnt <- c(countSig(res),countSig(resO),countSig(resD), countSig(resV),countSig(resM))

pbar2 <- data.frame(countSig(res),countSig(resO),countSig(resD), countSig(resV),countSig(resM))

labs <- c("Change vs Static", "Origin (M) vs Origin (V)", "Destination (M) vs Destination (V)", "(V -> M) vs (V -> V)", "(M -> V) vs (M -> M)")

nms <- c("(1)", "(2)", "(3)", "(4)", "(5)")

colrs <- c(rgb(139/255, 101/255, 8/255, 1), rgb(139/255, 101/255, 8/255, 0.8), rgb(139/255, 101/255, 8/255, 0.6))

res <- res[res$baseMean > 0,]
write.table(res, file="SvC.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=FALSE)
resD <- resD[resD$baseMean > 0,]
write.table(resD, file="../Networks/DvD.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=FALSE)
res <- res[res$baseMean > 0,]
write.table(res, file="SvC.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=FALSE)

pdf("RNA_p_bar_scale.pdf", width=3, height=3)
par(cex=0.7)
barplot(as.matrix(pbar2), col=colrs, ylab="Gene Count", xlab="Tests (1-5)", main="p < 0.05, Genes per Differential Test")
dev.off()

pdf("Test_results.pdf", width=6, height=9)
par(mfrow=c(3,2), cex=0.7)
plotMA(res, ylim=c(-2,2), main="(1). Change vs Static")
plotMA(resO, ylim=c(-2,2), main="(2). Origin (M) vs Origin (V)")
plotMA(resD, ylim=c(-2,2), main="(3). Destination (M) vs Destination (V)")
plotMA(resV, ylim=c(-2,2), main="(4). (V -> M) vs (V -> V)")
plotMA(resM, ylim=c(-2,2), main="(5). (M -> V) vs (M -> M)")
barplot(pcnt, names=nms_lab, col="darkred", ylab="Gene Counts", xlab="Tests (1-5)", main="p < 0.05, Genes per Differential Test")
dev.off()

mth_p <- read.table("gexp_genome/meth_nm_promo.bed", sep="\t", header=FALSE)
mth_g <- read.table("meth_nams.bed", sep="\t", header=FALSE)
mth_g <- mth_g[mth_g$V4 == "gene",]
mth_g <- mth_g[order(mth_g$V12),]
mg <- (mth_g[,7:11] * 150)/(mth_g[,3] - mth_g[,2])

mth_p1K <- mth_p[seq(2, 53902, 2),]
mth_p1C <- mth_p[seq(1, 53901, 2),]

mth_p1K <- mth_p1K[order(mth_p1K$V12),]
mth_p1C <- mth_p1C[order(mth_p1C$V12),]

#10f 10m 6f 6m VM MM VV MV

fc_MV_mK <- SafeFC(p1K[,4],p1K[,2])
fc_VM_mK <- SafeFC(p1K[,1],p1K[,3])

fc_O_mK <- SafeFC(p1K[,3] + p1K[,1], p1K[,2] + p1K[,2])
fc_D_mK <- SafeFC(p1K[,1] + p1K[,2], p1K[,3] + p1K[,4])
fc_S_mK <- SafeFC(p1K[,1] + p1K[,4], p1K[,3] + p1K[,2])

fc_MV_mC <- SafeFC(p1C[,4],p1C[,2])
fc_VM_mC <- SafeFC(p1C[,1],p1C[,3])

fc_O_mC <- SafeFC(p1C[,3] + p1C[,1], p1C[,2] + p1C[,2])
fc_D_mC <- SafeFC(p1C[,1] + p1C[,2], p1C[,3] + p1C[,4])
fc_S_mC <- SafeFC(p1C[,1] + p1C[,4], p1C[,3] + p1C[,2])

fc_MV_m <- SafeFC(mg[,4],mg[,2])
fc_VM_m <- SafeFC(mg[,1],mg[,3])

fc_O_m <- SafeFC(mg[,3] + mg[,1], mg[,2] + mg[,4])
fc_D_m <- SafeFC(mg[,1] + mg[,2], mg[,3] + mg[,4])
fc_S_m <- SafeFC(mg[,1] + mg[,4], mg[,3] + mg[,2])




countSig <- function(x)
{
	tst <- x$pvalue
	tst[is.na(tst)] <- 1
	tsA <- tst[tst <= 0.05 & tst > 0.01]
	tsB <- tst[tst <= 0.01 & tst > 0.001]
	tsC <- tst[tst <= 0.001]
	c(length(tsA), length(tsB), length(tsC))
}

SafeFC <- function(a,b)
{
	c <- pmin(a/b, 1000)
	c[is.na(c)] <- 0
	c <- log(c)
	c <- pmax(c, -max(c))
}


