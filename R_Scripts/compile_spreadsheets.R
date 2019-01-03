nms <- c("MM", "MV", "VM", "VV")
nm2 <- c("mm", "mv", "vm", "vv")

sms <- c("bp4", "cc4", "mf4")

bps <- list()
ccs <- list()
mfs <- list()

bps_p <- list()
ccs_p <- list()
mfs_p <- list()

for(i in 1:(length(nms)))
{
	bps[[i]] <- read.table(paste("davids_", nms[i], "/", nm2[i], "_neg_bp4.txt", sep=""), header=TRUE, sep="\t")
	ccs[[i]] <- read.table(paste("davids_", nms[i], "/", nm2[i], "_neg_cc4.txt", sep=""), header=TRUE, sep="\t")
	mfs[[i]] <- read.table(paste("davids_", nms[i], "/", nm2[i], "_neg_mf4.txt", sep=""), header=TRUE, sep="\t")

	bps_p[[i]] <- read.table(paste("davids_", nms[i], "/", nm2[i], "_pos_bp4.txt", sep=""), header=TRUE, sep="\t")
	ccs_p[[i]] <- read.table(paste("davids_", nms[i], "/", nm2[i], "_pos_cc4.txt", sep=""), header=TRUE, sep="\t")
	mfs_p[[i]] <- read.table(paste("davids_", nms[i], "/", nm2[i], "_pos_mf4.txt", sep=""), header=TRUE, sep="\t")
}

bps <- cllps(bps, length(nms))
ccs <- cllps(ccs, length(nms))
mfs <- cllps(mfs, length(nms))

bps_p <- cllps(bps_p, length(nms))
ccs_p <- cllps(ccs_p, length(nms))
mfs_p <- cllps(mfs_p, length(nms))

for(i in 1:length(nms))
{
	bps[[i]] <- merge(bps[[i]], bps_p[[i]], all=TRUE, by="Term", suffixes=c(".Neg", ".Pos"))
	ccs[[i]] <- merge(ccs[[i]], ccs_p[[i]], all=TRUE, by="Term", suffixes=c(".Neg", ".Pos"))
	mfs[[i]] <- merge(mfs[[i]], mfs_p[[i]], all=TRUE, by="Term", suffixes=c(".Neg", ".Pos"))
}

for(i in 2:length(nms))
{
	bps[[1]] <- merge(bps[[1]], bps[[i]], all=TRUE, by="Term", suffixes=c("", nms[i]))
	ccs[[1]] <- merge(ccs[[1]], ccs[[i]], all=TRUE, by="Term", suffixes=c("", nms[i]))
	mfs[[1]] <- merge(mfs[[1]], mfs[[i]], all=TRUE, by="Term", suffixes=c("", nms[i]))
}

bio_proc <- bps[[1]]
cell_com <- ccs[[1]]
mole_func <- mfs[[1]]

bio_proc[is.na(bio_proc)] <- 0
cell_com[is.na(cell_com)] <- 0
mole_func[is.na(mole_func)] <- 0

write.table(bio_proc, file="INTER_Bio_proc.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(cell_com, file="INTER_cell_comp.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(mole_func, file="INTER_molec_func.tsv", sep="\t", row.names=FALSE, quote=FALSE)

cllps <- function(x, k)
{
	for(i in 1:k)
	{
		x[[i]] <- data.frame(Term=x[[i]]$Term, Count=x[[i]]$Count, Benjamini=x[[i]]$Benjamini)
	}
	x
}