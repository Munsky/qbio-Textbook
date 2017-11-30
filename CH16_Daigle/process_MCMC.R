#!/usr/bin/env Rscript

####
#### Runs CODA analysis of MCMC output and writes/plots results to files
####

args.vec <- commandArgs(TRUE)
if(length(args.vec) < 6) { stop("Usage: process_MCMC.R CODAOUT(character) CODAIND(character) NPARS(numeric) OUTFILE(character) FIGURENAME(character) PLOTEXTRA(boolean) [TRUEPARS1](numeric)...[TRUEPARSN](numeric)") }
codaout <- args.vec[1]
codaind <- args.vec[2]
npars <- as.numeric(args.vec[3])
outfile <- args.vec[4]
figure <- args.vec[5]
plotextra <- as.logical(args.vec[6])
if(length(args.vec) >= 7) {
	truepars <- as.numeric(args.vec[7:length(args.vec)])
}

### Read in CODA files
require(coda)
mcmc <- read.coda(codaout, codaind, quiet=TRUE)

### Create plots
if(grepl(".pdf", figure, fixed=TRUE)) {
	pdf(figure, width=9, height=5)
} else {
	jpeg(figure, width=9, height=5, units="in", res=150)
}
par(mfcol=c(1,npars), mar=c(5.1,3.6,2.1,1.5), mgp=c(2.5,1,0))
xmax <- c(10,1.5,2,1.5)
for(i in seq(npars)) {
	plot(density(mcmc[,i], n=512, from=0, to=xmax[i]), log="", ylim=c(0,8.5), col="black", lty=1, lwd=3, main=colnames(mcmc)[i], xlab="", ylab="", cex.axis=1.25, cex.main=1.5)
	if(i==1) { title(ylab="Density", cex.lab=1.35) }
	if(exists("truepars")) { abline(v=truepars[i], col="black", lty=2) }
}
mtext("Parameter Value", side=1, outer=TRUE, line=-2, cex=.9)
dev.off()

if(plotextra) {
	jpeg("mcmcplot.jpg", width=7, height=7, units="in", res=150)
	plot(mcmc[, 1:npars])
	dev.off()

	jpeg("cumuplot.jpg", width=7, height=7, units="in", res=150)
	cumuplot(mcmc[, (npars+1):(2*npars)])
	dev.off()

	jpeg("acfplot.jpg", width=7, height=7, units="in", res=150)
	acfplot(mcmc[, (npars+1):(2*npars)])
	dev.off()

	jpeg("autoplot.jpg", width=7, height=7, units="in", res=150)
	autocorr.plot(mcmc[, (npars+1):(2*npars)])
	dev.off()
}

### Output parameter values
lower <- apply(mcmc[, 1:npars], 2, function(x) quantile(x, probs=.16))
upper <- apply(mcmc[, 1:npars], 2, function(x) quantile(x, probs=.84))
average <- summary(mcmc)$statistics[1:npars, "Mean"]
write.table(round(rbind("Lower"=lower,"Upper"=upper,"Mean"=average), 6), file=outfile, quote=FALSE, sep="\t")

