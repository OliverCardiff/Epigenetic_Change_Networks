library(DESeq2)

nms <- c("1F", "1M", "3F", "3M", "6F", "6M", "8F", "8M", "10F", "10M", "12F", "12M")

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

colnames(accu) <- nms

destination <- c("V", "V", "V", "V", "V", "V", "M", "M", "M", "M", "M", "M")

origin <- c("V", "M", "V", "M", "V", "M","V", "M", "V", "M", "V", "M")

transp <- c("VV", "MV", "VV", "MV", "VV", "MV", "VM", "MM", "VM", "MM", "VM", "MM")

change <- c("S", "C", "S", "C", "S", "C", "C", "S", "C", "S", "C", "S")

MV <- c(2,4,6)
MM <- c(8, 10, 12)
VM <- c(7, 9, 11)
VV <- c(1, 3, 5)

coldata <- cbind(origin=(origin), destination=(destination), transp=transp, change=change)

cVdata <- coldata[c(VM, VV),]
cMdata <- coldata[c(MV,MM),]

colnames(accu) <- NULL

accuV <- accu[,c(VM, VV)]
accuM <- accu[,c(MV,MM)]

mi_ddsM <- DESeqDataSetFromMatrix(countData = accuM,
                              colData = cMdata,
                              design= ~ transp)

mi_ddsM <- DESeq(mi_ddsM)

mi_ddsV <- DESeqDataSetFromMatrix(countData = accuV,
                              colData = cVdata,
                              design= ~ transp)

mi_ddsV <- DESeq(mi_ddsV)

mi_dds <- DESeqDataSetFromMatrix(countData = accu,
                              colData = coldata,
                              design= ~ origin + destination + change)

mi_dds <- DESeq(mi_dds)


mi_res <- lfcShrink(mi_dds, coef="change_S_vs_C", type="normal")
mi_resD <- lfcShrink(mi_dds, coef="destination_V_vs_M", type="normal")
mi_resO <- lfcShrink(mi_dds, coef="origin_V_vs_M", type="normal")

mi_resV <- lfcShrink(mi_ddsV, coef="transp_VV_vs_VM", type="normal")
mi_resM <- lfcShrink(mi_ddsM, coef="transp_MV_vs_MM", type="normal")

mi_pcnt <- c(countSig(mi_res),countSig(mi_resO),countSig(mi_resD), countSig(mi_resV),countSig(mi_resM))
labs <- c("Change vs Static", "Origin (M) vs Origin (V)", "Destination (M) vs Destination (V)", "(V -> M) vs (V -> V)", "(M -> V) vs (M -> M)")

mi_pbar <- data.frame(countSig(mi_res),countSig(mi_resO),countSig(mi_resD), countSig(mi_resV),countSig(mi_resM))

colrs <- c(rgb(0, 0.703, 0, 1), rgb(0, 0.703, 0, 0.8), rgb(0, 0.703, 0, 0.6))

nms_lab <- c("(1)", "(2)", "(3)", "(4)", "(5)")

pdf("miRNA_p_bar_scale.pdf", width=3, height=3)
par(cex=0.7)
barplot(as.matrix(mi_pbar), names=nms_lab, col=colrs, ylab="miRNA Count", xlab="Tests (1-5)", main="p < 0.05, miRNA Read Count Diff.")
dev.off()

pdf("mi_test_results.pdf", width=6, height=9)
par(mfrow=c(3,2), cex=0.7)
plotMA(mi_res, ylim=c(-2,2), main="(1). Change vs Static")
plotMA(mi_resO, ylim=c(-2,2), main="(2). Origin (M) vs Origin (V)")
plotMA(mi_resD, ylim=c(-2,2), main="(3). Destination (M) vs Destination (V)")
plotMA(mi_resV, ylim=c(-2,2), main="(4). (V -> M) vs (V -> V)")
plotMA(mi_resM, ylim=c(-2,2), main="(5). (M -> V) vs (M -> M)")
barplot(mi_pcnt, names=nms_lab, col="darkred", ylab="miRNA Counts", xlab="Tests (1-5)", main="Counts of (p < 0.05) miRNA per Test")
dev.off()



