mth_p <- read.table("gexp_genome/meth_nm_promo.bed", sep="\t", header=FALSE)
mth_g <- read.table("meth_nams.bed", sep="\t", header=FALSE)
trpl <- read.table("gexp_genome/transp_uniq.txt", sep="\t", header=FALSE)
tser <- read.table("gexp_genome/tser_uniq.txt", sep="\t", header=FALSE)

mth_g <- mth_g[mth_g$V4 == "gene",]
mth_g <- mth_g[order(mth_g$V12),]
mg <- (mth_g[,7:11] * 150)/(mth_g[,3] - mth_g[,2])

mth_p1K <- mth_p[seq(2, 53902, 2),]
mth_p1C <- mth_p[seq(1, 53901, 2),]

mth_p1K <- mth_p1K[order(mth_p1K$V12),]
mth_p1C <- mth_p1C[order(mth_p1C$V12),]


mts <- tser[,44:83]

g_csm <- colSums(mt)
g_csmts <- colSums(mts)

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

arthMean <- apply(mts, 1, FUN=function(x){gm_mean(x)})
rats <- apply(mts, 2, FUN=function(x){x/arthMean})

facts <- rep(0, ncol(rats))

for(i in 1:ncol(rats))
{
	tst <- rats[,i]
	tst <- tst[ tst > 0 ]
	facts[i] <- median(tst)

	mts[,i] <- mts[,i] / facts[i]
}

p1K <- mth_p1K[,7:11]
p1C <- mth_p1C[,7:11]

p1K <- apply(p1K, 2, FUN=function(x){(x/sum(x)) * 1000000})
p1C <- apply(p1C, 2, FUN=function(x){(x/sum(x)) * 1000000})
mg <- apply(mg, 2, FUN=function(x){(x/sum(x)) * 1000000})

row.names(mg) <- mth_g$V12
row.names(mt) <- trpl$V1
row.names(mts) <- tser$V1
row.names(p1K) <- mth_p1K$V12
row.names(p1C) <- mth_p1C$V12

names(g_csmts) <- c(1:40)
names(g_csm2ts) <- c(1:40)
names(g_csm) <- c(1:14)
names(g_csm2) <- c(1:14)
g_csm2 <- colSums(mt)
g_csm2ts <- colSums(mts)


write.table(mg, file="clusters/meth_gene_body_norm.txt", sep="\t", quote=FALSE, col.names=FALSE)
write.table(mt, file="clusters/gexp_trpl_norm.txt", sep="\t", quote=FALSE, col.names=FALSE)
write.table(mts, file="clusters/gexp_tser_norm.txt", sep="\t", quote=FALSE, col.names=FALSE)
write.table(p1K, file="clusters/meth_1Kb_promo_norm.txt", sep="\t", quote=FALSE, col.names=FALSE)
write.table(p1C, file="clusters/meth_34bp_promo_norm.txt", sep="\t", quote=FALSE, col.names=FALSE)

VM <- c(1,3,11)
MM <- c(2, 4, 12)
MV <- c(6, 8, 10)
VV <- c(5, 7, 9)

p1K_c <- p1K
p1C_c <- p1C
mg_c <- mg

for(i in 1:4)
{
	p1K_c[,i] <- pmax(p1K[,i] - p1K[,5], 0)
	p1C_c[,i] <- pmax(p1C[,i] - p1C[,5], 0)
	mg_c[,i] <- pmax(mg[,i] - mg[,5], 0)
}
mts_tst <- mts

for(i in 1:nrow(mts))
{
	mts_tst[i,] <- sample(mts[i,])
}

mts_t <- apply(mts_tst, 1, FUN=function(x){cor(1:40, x, method="kendall")})
mts_t[is.na(mts_t)] <- 0

mts_sp <- apply(mts, 1, FUN=function(x){cor(1:40, x, method="spearman")})
mts_sp[is.na(mts_sp)] <- 0

mts_pe <- apply(mts, 1, FUN=function(x){cor(1:40, x, method="pearson")})
mts_pe[is.na(mts_pe)] <- 0

mts_kd <- apply(mts, 1, FUN=function(x){cor(1:40, x, method="kendall")})
mts_kd[is.na(mts_kd)] <- 0

plot(density(mts_kd))
points(density(mts_t), type="l")
abline(v=0)


VM_m <- rowMeans(mt[,VM])
MM_m <- rowMeans(mt[,MM])
VV_m <- rowMeans(mt[,VV])
MV_m <- rowMeans(mt[,MV])

VM_sd <- apply(mt[,VM], 1, sd)
MM_sd <- apply(mt[,MM], 1, sd)
VV_sd <- apply(mt[,VV], 1, sd)
MV_sd <- apply(mt[,MV], 1, sd)

p_MV <- pnorm(VM_m, mean=MM_m, sd=MM_sd, lower.tail=FALSE)
p_MM <- pnorm(VM_m, mean=MV_m, sd=MV_sd, lower.tail=FALSE)
p_VV <- pnorm(VV_m, mean=MM_m, sd=MM_sd, lower.tail=FALSE)
p_VM <- pnorm(VM_m, mean=VV_m, sd=VV_sd, lower.tail=FALSE)

samps3 <- c("V->V vs M->M", "V->M vs V->V", "M->V vs M->M", "M->V vs V->M")
colr <- c(rgb(0.2,0.5,0.8,0.5), rgb(0.6,0.7,0.2,0.5), rgb(0.2,0.4,0.2,0.5), rgb(0.9,0.2,0.3,0.5), rgb(0.9,0.2,0.9,0.5), rgb(0.9,0.9,0.1,0.5))

pdf("PV_test_diffs.pdf", width=6,height=5)
par(cex=0.7)
plot(density(p_VV, cut=0), col=colr[1], lwd=2, xlim=c(0,1.3), main="Z-Score Density of Differential Tests",
xlab="left -> right, normal distribution z-scores")
points(density(p_VM, cut=0), type="l", col=colr[2], lwd=2)
points(density(p_MV, cut=0), type="l", col=colr[3], lwd=2)
points(density(p_MM, cut=0), type="l", col=colr[4], lwd=2)
abline(v=0.05, lwd=2, lty=2)
abline(v=0.95, lwd=2, lty=2)
legend("topright", col=c(colr[1:4], "black"), legend=c(samps3, "z < 0.05 || z > 0.95"), lwd=c(2,2,2,2,2), lty=c(1,1,1,1,2))
dev.off()

pdf("FC_test_diffs.pdf", width=6,height=5)
par(cex=0.7)
plot(density(fc_VV, cut=0), col=colr[1], lwd=2, ylim=c(0,1.6), main="Log FC Density Differential Tests",
xlab="Log Fold Change")
points(density(fc_VM, cut=0), type="l", col=colr[2], lwd=2)
points(density(fc_MV, cut=0), type="l", col=colr[3], lwd=2)
points(density(fc_MM, cut=0), type="l", col=colr[4], lwd=2)
abline(v=-6.7, lwd=2, lty=2)
abline(v=6.7, lwd=2, lty=2)
legend("topright", bg="white", col=c(colr[1:4], "black"), legend=c(samps3, "FC < 1000 || FC > 1000"), lwd=c(2,2,2,2,2), lty=c(1,1,1,1,2))
dev.off()

vfix <- p_VV

for(i in 1:length(p_VV))
{
	if(p_VV[i] > 0.5)
	{
		vfix[i] <- 1 - p_VV[i]
	}
}

pdf("fc_vs_zs.pdf", width=6,height=5)
plot(fc_VV, -log10(vfix), col=rgb(0.1,0.1,0.1,0.2), main="z-score vs log FC, V->V vs M->M",
xlab="log FC", ylab="Negative Log10 P-Value", ylim=c(0,10))
dev.off()

fc_MV <- SafeFC(MV_m,MM_m)
fc_VV <- SafeFC(VV_m,MM_m)
fc_VM <- SafeFC(VM_m,VV_m)
fc_MM <- SafeFC(VM_m,MV_m)

Ac_MV <- (MV_m - MM_m)
Ac_VV <- (VV_m - MM_m)
Ac_VM <- (VM_m - VV_m)
Ac_MM <- (VM_m - MV_m)

#10f 10m 6f 6m VM MM VV MV

fc_MV_mK <- SafeFC(p1K[,4],p1K[,2])
fc_VV_mK <- SafeFC(p1K[,3],p1K[,2])
fc_VM_mK <- SafeFC(p1K[,1],p1K[,3])
fc_MM_mK <- SafeFC(p1K[,4],p1K[,1])

fc_MV_mC <- SafeFC(p1C[,4],p1C[,2])
fc_VV_mC <- SafeFC(p1C[,3],p1C[,2])
fc_VM_mC <- SafeFC(p1C[,1],p1C[,3])
fc_MM_mC <- SafeFC(p1C[,4],p1C[,1])

fc_MV_m <- SafeFC(mg[,4],mg[,2]) #rev
fc_VV_m <- SafeFC(mg[,3],mg[,2]) #actually MM
fc_VM_m <- SafeFC(mg[,1],mg[,3]) #rev
fc_MM_m <- SafeFC(mg[,4],mg[,1]) #acutally VV

fc_MV_m <- SafeFC(mg[,4],mg[,2])
fc_VV_m <- SafeFC(mg[,3],mg[,2])
fc_VM_m <- SafeFC(mg[,1],mg[,3])
fc_MM_m <- SafeFC(mg[,4],mg[,1])

rsm <- rowMeans(mg[,1:4])
p1sm <- rowMeans(p1K[,1:4])
p2sm <- rowMeans(p1C[,1:4])

exp <- log(rowMeans(cbind(VM_m, MM_m, VV_m, MV_m)))

plot(fc_MV_m, mts_sp, pch=".", ylim=c(-1,1), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_MM_m, mts_sp, pch=".", ylim=c(-1,1), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_VV_m, mts_sp, pch=".", ylim=c(-1,1), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_VM_m, mts_sp, pch=".", ylim=c(-1,1), xlim=c(-8,8))
abline(h=0, v=0, col="red")

plot(fc_MV_m, fc_MV, pch=".", ylim=c(-4,4), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_MM_m, fc_MM, pch=".", ylim=c(-4,4), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_VM_m, fc_VM, pch=".", ylim=c(-4,4), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_VV_m, fc_VV, pch=".", ylim=c(-4,4), xlim=c(-8,8))
abline(h=0, v=0, col="red")

plot(fc_MV_m, Ac_MV, pch=".", ylim=c(-4000,4000), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_MM_m, Ac_MM, pch=".", ylim=c(-4000,4000), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_VM_m, Ac_VM, pch=".", ylim=c(-4000,4000), xlim=c(-8,8))
abline(h=0, v=0, col="red")
plot(fc_VV_m, Ac_VV, pch=".", ylim=c(-4000,4000), xlim=c(-8,8))
abline(h=0, v=0, col="red")

d1 <- sample(exp)
d2 <- sample(rsm)

plot(exp, rsm, ylim=c(0, 4000))
plot(d1,d2, ylim=c(0, 4000))

ot20 <- Run2DQQ(rsm, exp, 20, 40)
ot40 <- Run2DQQ(rsm, exp, 40, 40)

ot1k20 <- Run2DQQ(p1sm, exp, 20, 40)
ot1k40 <- Run2DQQ(p1sm, exp, 40, 40)

otMM20 <- Run2DQQ(fc_MM_m, fc_MM, 20, 40)
otMM40 <- Run2DQQ(fc_MM_m, fc_MM, 40, 40)

KotMM20 <- Run2DQQ(fc_MM_mK, fc_MM, 20, 40, trim=1)
KotMM40 <- Run2DQQ(fc_MM_mK, fc_MM, 40, 40, trim=1)

ot <- Run2DQQ(rsm, exp, 10, 40)
otT <- Run2DQQ(rsm, mts_sp, 10, 40)
ot1k <- Run2DQQ(p1sm, exp, 10, 40)
ot1c <- Run2DQQ(p2sm, exp, 10, 40)

otMV <- Run2DQQ(fc_MV_m, fc_MV, 10, 40)
otMM <- Run2DQQ(fc_MM_m, fc_MM, 10, 40)
otVV <- Run2DQQ(fc_VV_m, fc_VV, 10, 40)
otVM <- Run2DQQ(fc_VM_m, fc_VM, 10, 40)

otMVt <- Run2DQQ(fc_MV_m, mts_sp, 10, 40)
otMMt <- Run2DQQ(fc_MM_m, mts_sp, 10, 40)
otVVt <- Run2DQQ(fc_VV_m, mts_sp, 10, 40)
otVMt <- Run2DQQ(fc_VM_m, mts_sp, 10, 40)

KotMV <- Run2DQQ(fc_MV_mK, fc_MV, 10, 40, trim=1)
KotMM <- Run2DQQ(fc_MM_mK, fc_MM, 10, 40, trim=1)
KotVV <- Run2DQQ(fc_VV_mK, fc_VV, 10, 40, trim=1)
KotVM <- Run2DQQ(fc_VM_mK, fc_VM, 10, 40, trim=1)

CotMV <- Run2DQQ(fc_MV_mC, fc_MV, 10, 40, trim=1)
CotMM <- Run2DQQ(fc_MM_mC, fc_MM, 10, 40, trim=1)
CotVV <- Run2DQQ(fc_VV_mC, fc_VV, 10, 40, trim=1)
CotVM <- Run2DQQ(fc_VM_mC, fc_VM, 10, 40, trim=1)

pdf("quants/ot.pdf", width=4, height=8)
PlotStack(ot, 0.67, "Log Mean TPM Quantiles")
dev.off()

pdf("quants/ot1k.pdf", width=4, height=8)
PlotStack(ot1k, 0.67, "Log Mean TPM Quantiles")
dev.off()

pdf("quants/ot1c.pdf", width=4, height=8)
PlotStack(ot1c, 0.67, "Log Mean TPM Quantiles")
dev.off()

pdf("quants/otVMt.pdf", width=4, height=8)
PlotStack(otVMt, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/otMMt.pdf", width=4, height=8)
PlotStack(otMMt, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/otMVt.pdf", width=4, height=8)
PlotStack(otMVt, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/otVVt.pdf", width=4, height=8)
PlotStack(otVVt, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()

pdf("quants/otVM.pdf", width=4, height=8)
PlotStack(otVM, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/otMM.pdf", width=4, height=8)
PlotStack(otMM, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/otMV.pdf", width=4, height=8)
PlotStack(otMV, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/otVV.pdf", width=4, height=8)
PlotStack(otVV, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()

pdf("quants/KotVM.pdf", width=4, height=8)
PlotStack(KotVM, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/KotMM.pdf", width=4, height=8)
PlotStack(KotMM, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/KotMV.pdf", width=4, height=8)
PlotStack(KotMV, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/KotVV.pdf", width=4, height=8)
PlotStack(KotVV, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()

pdf("quants/CotVM.pdf", width=4, height=8)
PlotStack(CotVM, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/CotMM.pdf", width=4, height=8)
PlotStack(CotMM, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/CotMV.pdf", width=4, height=8)
PlotStack(CotMV, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()
pdf("quants/CotVV.pdf", width=4, height=8)
PlotStack(CotVV, 0.67, "Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()


PlotStack(otMMt, 0.7)
PlotStack(otMVt, 0.7)
PlotStack(otVVt, 0.7)
PlotStack(otVMt, 0.7)

source("m_func.R")

pdf("Spatial_qmap.pdf", width=8, height=8)
par(mfrow=c(2,2), cex=0.7)
Plot2DQQ(ot, inc=0.15, ylb="Gene Body Methylation RPM", xlb="All-Sample Mean RNASeq TPM Quantiles")
Plot2DQQ(otMM, inc=0.15, ylb="Gene Body Methylation Log FC", xlb="Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
Plot2DQQ(ot1k, inc=0.15, ylb="1kb TSS Upstream Methylation RPM", xlb="All-Sample Mean RNASeq TPM Quantiles")
Plot2DQQ(KotMM, inc=0.15, ylb="1kb TSS Upstream Methylation Log FC", xlb="Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()

pdf("Grid_qmap.pdf", width=8, height=8)
par(mfrow=c(2,2), cex=0.7)
Plot2DQQGrid(ot, inc=0.2, ylb="Gene Body Meth. RPM Quantiles", xlb="All-Sample Mean RNASeq TPM Quantiles")
Plot2DQQGrid(otMM, inc=0.2, ylb="Gene Body Meth. Log FC Quantiles", xlb="Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
Plot2DQQGrid(ot1k, inc=0.2, ylb="1kb TSS Upstream Meth. RPM Quantiles", xlb="All-Sample Mean RNASeq TPM Quantiles")
Plot2DQQGrid(KotMM, inc=0.2, ylb="1kb TSS Upstream Meth. RPM Quantiles", xlb="Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()

pdf("tile_qmap.pdf", width=8, height=8)
par(mfrow=c(2,2), cex=0.7)
PlotFitQQ(ot20, ot40, ylb="Gene Body Methylation RPM", xlb="All-Sample Mean RNASeq TPM Quantiles")
PlotFitQQ(otMM20, otMM40, ylb="Gene Body Methylation Log FC", xlb="Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
PlotFitQQ(ot1k20, ot1k40, ylb="1kb TSS Upstream Methylation RPM", xlb="All-Sample Mean RNASeq TPM Quantiles")
PlotFitQQ(KotMM20, KotMM40, ylb="1kb TSS Upstream Methylation Log FC", xlb="Gene Expression Quantiles: TPM A / (TPM_A + TPM_B)")
dev.off()

sfc <- SafeFC(p1K[,2], p1K[,4])
summary(sfc)

SafeFC <- function(a,b)
{
	c <- pmin(a/b, 1000)
	c[is.na(c)] <- 0
	c <- log(c)
	c <- pmax(c, -max(c))
}

ChangeFC <- function(a,b)
{
	c <- abs(a-b)/pmax(a,b)
	c[is.na(c)] <- 0
	c
}



