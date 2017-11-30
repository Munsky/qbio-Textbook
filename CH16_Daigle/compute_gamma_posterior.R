#!/usr/bin/env Rscript

####
#### Print out and plot Bayesian gamma posterior distributions given MLEs and priors
####

args.vec <- commandArgs(TRUE)
if(length(args.vec) < 5) { stop("Usage: compute_gamma_posterior.R MLEFILE(character) PRIOR1(numeric) PRIOR2(numeric) OUTFILE(character) FIGURENAME(character) [TRUEPARS1](numeric)...[TRUEPARSN](numeric)") }
mlefile <- args.vec[1]
prior1 <- as.numeric(args.vec[2])
prior2 <- as.numeric(args.vec[3])
outfile <- args.vec[4]
figure <- args.vec[5]
if(length(args.vec) >= 6) {
	truepars <- as.numeric(args.vec[6:length(args.vec)])
}

### Read in MLEs; create Bayesian results matrix
mle <- as.matrix(read.table(mlefile, header=TRUE, sep="\t", as.is=TRUE))
bayes <- matrix(nrow=2, ncol=ncol(mle))

### Create plot and populate results
if(grepl(".pdf", figure, fixed=TRUE)) {
	pdf(figure, width=9, height=5)
} else {
	jpeg(figure, width=9, height=5, units="in", res=150)
}
par(mfcol=c(1,ncol(mle)), mar=c(5.1,3.6,2.1,1.5), mgp=c(2.5,1,0))
xmax <- c(10,1.5,2,1.5)
for(i in seq(ncol(mle))) {
	shape <- mle["R", i] + prior1
	rate <- mle["H", i] + prior1
	bayes[1, i] <- shape/rate
	curve(dgamma(x, shape, rate), from=0, to=xmax[i], n=101, log="", ylim=c(0,8.5), col="gray", lty=1, lwd=3, main=colnames(mle)[i], xlab="", ylab="", cex.axis=1.25, cex.main=1.5)
	if(i==1) { title(ylab="Density", cex.lab=1.35) }

	shape <- mle["R", i] + prior2
	rate <- mle["H", i] + prior2
	bayes[2, i] <- shape/rate
	curve(dgamma(x, shape, rate), from=0, to=xmax[i], n=101, log="", add=TRUE, col="black", lty=1, lwd=3)
	if(exists("truepars")) { abline(v=truepars[i], col="black", lty=2) }
}
mtext("Parameter Value", side=1, outer=TRUE, line=-2, cex=.9)
dev.off()

### Output parameter values
rownames(bayes) <- c("Mean1", "Mean2")
colnames(bayes) <- colnames(mle)
write.table(round(bayes,6), file=outfile, quote=FALSE, sep="\t")

