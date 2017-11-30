#!/usr/bin/env Rscript

####
#### Correctly scales MCMC parameters in stochInf output table given multiplier
####

args.vec <- commandArgs(TRUE)
if(length(args.vec) < 2) { stop("Usage: scale_MCMC_parameters.R TABLEFILE(character) NUMPARAMS(numeric) [MULT](numeric)") }
tablefile <- args.vec[1]
numparams <- as.numeric(args.vec[2])
if(length(args.vec) >= 3) {
	multiplier <- as.numeric(args.vec[3])
} else {
	multiplier <- 1
}

table <- read.table(tablefile, sep=" ", header=TRUE, as.is=TRUE, check.names=FALSE)
table <- table[, -ncol(table)]
table[, 2:(numparams+1)] <- table[, 2:(numparams+1)] * multiplier
table[, (numparams+2):(2*numparams+1)] <- table[, (numparams+2):(2*numparams+1)] + log(multiplier)
write.table(table, file=tablefile, sep=" ", row.names=FALSE, col.names=TRUE, quote=FALSE)

