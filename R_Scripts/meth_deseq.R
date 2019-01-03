mth_p <- read.table("gexp_genome/meth_nm_promo.bed", sep="\t", header=FALSE)
mth_g <- read.table("meth_nams.bed", sep="\t", header=FALSE)

mth_g <- mth_g[mth_g$V4 == "gene",]
mth_g <- mth_g[order(mth_g$V12),]
mg <- (mth_g[,7:10] * 150)/(mth_g[,3] - mth_g[,2])

mth_p1K <- mth_p[seq(2, 53902, 2),]
mth_p1C <- mth_p[seq(1, 53901, 2),]
mth_p1K <- mth_p1K[order(mth_p1K$V12),]
mth_p1C <- mth_p1C[order(mth_p1C$V12),]
p1K <- mth_p1K[,7:10]
p1C <- mth_p1C[,7:10]
row.names(p1K) <- mth_p1K$V12
row.names(p1C) <- mth_p1C$V12
row.names(mg) <- mth_g$V12

#thrd <- Qdisp(mg[,1], mg[,2], 10)

p1K <- floor(apply(p1K, 2, FUN=function(x){(x/sum(x)) * 1000000}))
p1C <- floor(apply(p1C, 2, FUN=function(x){(x/sum(x)) * 1000000}))
mg <- floor(apply(mg, 2, FUN=function(x){(x/sum(x)) * 1000000}))

origin <- c("V", "V", "M", "M", "V", "M", "X", "X", "X", "X")
destination <- c("V", "M", "V", "M", "X", "X", "V", "M", "X", "X")
change <- c("S", "C", "C", "S", "X", "X", "X", "X", "S", "C")

mg <- fourway(mg)
p1C <- fourway(p1C)
p1K <- fourway(p1K)

coldata <- cbind(origin=(origin), destination=(destination), change=change)

me_dds <- DESeqDataSetFromMatrix(countData = mg,
                              colData = coldata,
                              design= ~ origin + destination + change)

me_dds <- DESeq(me_dds)

pC_dds <- DESeqDataSetFromMatrix(countData = p1C,
                              colData = coldata,
                              design= ~ origin + destination + change)

pC_dds <- DESeq(pC_dds)

p1K_dds <- DESeqDataSetFromMatrix(countData = p1K,
                              colData = coldata,
                              design= ~ origin + destination + change)

p1K_dds <- DESeq(p1K_dds)

me_res <- lfcShrink(me_dds, coef="change_S_vs_C", type="normal")
me_resD <- lfcShrink(me_dds, coef="destination_V_vs_M", type="normal")
me_resO <- lfcShrink(me_dds, coef="origin_V_vs_M", type="normal")

p1K_res <- lfcShrink(p1K_dds, coef="change_S_vs_C", type="normal")
p1K_resD <- lfcShrink(p1K_dds, coef="destination_V_vs_M", type="normal")
p1K_resO <- lfcShrink(p1K_dds, coef="origin_V_vs_M", type="normal")

p1C_res <- lfcShrink(pC_dds, coef="change_S_vs_C", type="normal")
p1C_resD <- lfcShrink(pC_dds, coef="destination_V_vs_M", type="normal")
p1C_resO <- lfcShrink(pC_dds, coef="origin_V_vs_M", type="normal")

me_pcnt <- c(countSig(me_res),countSig(me_resO),countSig(me_resD))
pC_pcnt <- c(countSig(p1C_res),countSig(p1C_resO),countSig(p1C_resD))
p1K_pcnt <- c(countSig(p1K_res),countSig(p1K_resO),countSig(p1K_resD))

all_pcnt <- cbind(me_pcnt, pC_pcnt, p1K_pcnt)

me_pbar <- data.frame(countSig(me_res),countSig(me_resO),countSig(me_resD))
colrs <- c(rgb(67/255, 110/255, 238/255, 1), rgb(67/255, 110/255, 238/255, 0.8), rgb(67/255, 110/255, 238/255, 0.6))

nms_lab <- c("(1)", "(2)", "(3)")

pdf("meth_p_bar_scale.pdf", width=3, height=3)
par(cex=0.7)
barplot(as.matrix(me_pbar), names=nms_lab, col=colrs, ylab="Gene Count", xlab="Tests (1-3)", main="p < 0.05, Gene body Methylation Diff.")
dev.off()

bcol <- c("purple","darkblue","darkred")

pdf("me_test_results.pdf", width=6, height=6)
par(mfrow=c(2,2), cex=0.7)
plotMA(me_res, ylim=c(-5,5), main="(1). Change vs Static")
plotMA(me_resO, ylim=c(-5,5), main="(2). Origin (M) vs Origin (V)")
plotMA(me_resD, ylim=c(-5,5), main="(3). Destination (M) vs Destination (V)")
barplot(t(all_pcnt), ylim=c(0, 6500), beside=TRUE, names=nms_lab, col=bcol, ylab="miRNA Counts", xlab="Tests (1-5)", main="Counts of (p < 0.05) Genes per Test")
legend("topleft", col=bcol, pch=c(15,15,15), legend=c("Gene Body", "Promotor 100bp", "Promotor 1Kbp"), bty="n")
dev.off()

write.table(me_res, file="deseqOut/meth_SvC.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(me_resO, file="deseqOut/meth_OvO.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(me_resD, file="deseqOut/meth_DvD.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(p1C_res, file="deseqOut/p1C_SvC.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(p1C_resO, file="deseqOut/p1C_OvO.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(p1C_resD, file="deseqOut/p1C_DvD.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(p1K_res, file="deseqOut/p1K_SvC.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(p1K_resO, file="deseqOut/p1K_OvO.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(p1K_resD, file="deseqOut/p1K_DvD.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(mi_res, file="deseqOut/mi_SvC.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(mi_resO, file="deseqOut/mi_OvO.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(mi_resD, file="deseqOut/mi_DvD.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(mi_resM, file="deseqOut/mi_MvV.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(mi_resV, file="deseqOut/mi_VvM.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(res, file="deseqOut/g_SvC.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(resO, file="deseqOut/g_OvO.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(resD, file="deseqOut/g_DvD.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(resM, file="deseqOut/g_MvV.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(resV, file="deseqOut/g_VvM.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

fourway <- function(x, d=10)
{
	o1 <- QD2 (x[,1],x[,2])
	o2 <- QD2 (x[,3],x[,4])
	
	d1 <- QD2 (x[,1],x[,3])
	d2 <- QD2 (x[,2],x[,4])
	
	s1 <- QD2 (x[,1],x[,4])
	s2 <- QD2 (x[,2],x[,3])

	floor(cbind(x,o1,o2,d1,d2,s1,s2))
}

QD2 <- function(x,y)
{
	z <- cbind(x,y)
	zM <- apply(z, 1, mean)
	zSD <- apply(z, 1, sd)

	w <- pmax(0,rnorm(length(zM), zM, zSD))
}

Qdisp <- function(x, y, div=10)
{
	z <- cbind(x,y)
	zM <- apply(z, 1, mean)
	zSD <- apply(z, 1, sd)

	qnts <- seq(1/div, (div-1)/div, 1/div)
	qn1 <- quantile(zM, qnts)

	oldSlct <- rep(FALSE, length(zM))

	fin <- vector()

	for(i in 1:div)
	{
		Low <- min(zM)
		High <- max(zM)
		if(i > 1)
		{
			Low <- qn1[i - 1]
		}
		if(i < div)
		{
			High <- qn1[i]
		}

		slcts <- zM >= Low & zM <= High

		s2 <- xor(oldSlct, slcts)

		slcts <- s2 & slcts

		qM <- zM[slcts]
		qS <- zSD[slcts]

		df <- data.frame(qM, qS)
		#fit <- lm(qM ~ qS, data=df)
		#ps <- predict(fit, df)
		#ps[is.na(ps)] <- qS[is.na(ps)]
		preds <- rnorm(length(qM), qM, qS)
		preds[is.na(preds)] <- qM[is.na(preds)]
		fin <- c(fin, preds)
		oldSlct <- slcts
	}

	pmax(0,fin)
}



