trpl <- read.table("nam_norms/transp_nm_cnv.txt", sep="\t", header=FALSE)
tser <- read.table("nam_norms/tseries_nm_cnv.txt", sep="\t", header=FALSE)
meth_mat <- read.table("nam_norms/meth_nm_cnv.txt", sep="\t", header=FALSE)

bd_mt <- read.table("bed_mat.txt", sep="\t", header=FALSE)

trpl <- trpl[meth_mat[,4] == "exon",]

VM <- c(1,3,11)
MM <- c(2, 4, 12)
MV <- c(6, 8, 10)
VV <- c(5, 7, 9)

VM <- VM + 3
VV <- VV + 3
MV <- MV + 3
MM <- MM + 3

vvSD <- apply(trpl[,VV], 1, sd)
mvSD <- apply(trpl[,MV], 1, sd)
vmSD <- apply(trpl[,VM], 1, sd)

allSD <- apply(trpl[,4:15], 1, sd)

mmSD <- apply(trpl[,MM], 1, sd)

mthSD <- apply(meth_mat[,4:7], 1, sd)

bd_mt <- bd_mt[bd_mt[,4] == "exon",]
bd2 <- (bd_mt[,7:10] * 150)/(bd_mt[,3] - bd_mt[,2])

mthSUM <- apply(bd2, 1, sum)

diffVM <- rowSums(trpl[,VM]) - rowSums(trpl[,VV])
d2 <- bd2[,1] - bd2[,4]

plot(log(abs(d2)), diffVM, pch=".", ylim=c(-5,5))

plot(mthSUM, allSD, pch=".", xlim=c(0,10000), ylim=c(0,0.7), col=rgb(0,0,0,0.4))
plot(mthSD, mvSD, pch=".", xlim=c(0,1.6), ylim=c(0,0.7), col=rgb(0,0,0,0.4))
plot(mthSD, vmSD, pch=".", xlim=c(0,1.6), ylim=c(0,0.7), col=rgb(0,0,0,0.4))
plot(mthSD, vvSD, pch=".", xlim=c(0,1.6), ylim=c(0,0.7), col=rgb(0,0,0,0.4))


