mth <- read.table("bed_mat.txt", sep="\t", header=FALSE)
trpl <- read.table("transp.bed", sep="\t", header=FALSE)
tser <- read.table("tseries.bed", sep="\t", header=FALSE)

nams <- mth$V4
VM <- c(7,9,17)
MM <- c(8, 10, 18)
MV <- c(12, 14, 16)
VV <- c(11, 13, 15)

mth <- correctit(mth)
trpl <- correctit(trpl)
tser <- correctit(tser)

trpl <- exlvls(trpl, nams)
tser <- exlvls(tser, nams)

write.table(trpl, "ex_lvls_trpl.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(tser, "ex_lvls_tser.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

exlvls <- function(x, nams)
{
	geCov <- vector()
	for(i in 1:nrow(x))
	{
		if(nams[i] == "exon")
		{
			x[i,] <- x[i,]/geCov
		}	

		if(nams[i] == "gene")
		{
			geCov <- x[i,]
		}
	}
	x[is.na(x)] <- 1
	x
}

correctit <- function(x)
{
	lns <- x[,3] - x[,2]

	scs <- colSums(x[,7:ncol(x)])

	for(i in 7:ncol(x))
	{
		x[,i] <- (x[,i] * 1000000)/scs[i-6]
		x[,i] <- (x[,i] * 150) / lns
	}
	x[,7:ncol(x)]
}

