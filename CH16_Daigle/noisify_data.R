#!/usr/bin/env Rscript

####
#### Adds noise to discretized data trajectory
####

args.vec <- commandArgs(TRUE)
if(length(args.vec) < 3) { stop("Usage: noisify_data.R DATAFILE(character) OUTFILE(character) SIGMA(numeric) [SEED](numeric)") }
datafile <- args.vec[1]
outfile <- args.vec[2]
sigma <- as.numeric(args.vec[3])
if(length(args.vec) >= 4) {
	seed <- as.numeric(args.vec[4])
} else {
	seed <- 101
}

data <- read.table(datafile, header=TRUE, sep=" ", row.names=NULL, as.is=TRUE)
set.seed(seed)
noise <- rnorm(nrow(data)*(ncol(data)-2), mean=0, sd=sigma)
data[, -c(1:2)] <- data[, -c(1:2)] + noise
write.table(data, file=outfile, sep=" ", row.names=FALSE, col.names=TRUE, quote=FALSE)

