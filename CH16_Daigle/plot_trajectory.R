#!/usr/bin/env Rscript

####
#### Plots data from single trajectory for specified species 
####

args.vec <- commandArgs(TRUE)
if(length(args.vec) < 3) { stop("Usage: plot_trajectory.R DATAFILE(character) FIGURENAME(character) TYPE(character) [SPECIES1__SPECIES2__...__SPECIESN](character)") }
datafile <- args.vec[1]
data <- read.table(datafile, header=TRUE, sep=" ", row.names=NULL, as.is=TRUE)
figure <- args.vec[2]
type <- args.vec[3]
if(length(args.vec) >= 4) {
	species <- unlist(strsplit(args.vec[4], "__"))
} else {
	species <- colnames(data)[-c(1:2)]
}

xmax <- max(data[, "Time"])
#ymax <- max(data[, colnames(data) %in% species])
ymax <- 25
cex <- 2
cex.axis <- 2
pch <- 16
lwd <- 4
lty <- "solid"
#colors <- rainbow(length(species))
colors <- rep(c("gray","black"), length.out=length(species))

if(grepl(".pdf", figure, fixed=TRUE)) {
	pdf(figure, width=9, height=5)
} else {
	jpeg(figure, width=9, height=5, units="in", res=150)
}
par(xaxs="i", yaxs="i", mar=c(4, 5.25, 1, 1.5))

plot(NA, NA, xlab="", ylab="", xlim=c(-.01*xmax,1.01*xmax), ylim=c(-.01*ymax,1.01*ymax), axes=FALSE)
for(i in seq_along(species)) {
	points(data[, "Time"], data[, which(colnames(data)==species[i])], type=type, lty=lty, lwd=lwd, pch=pch, cex=cex, col=colors[i])
}
axis(1, cex.axis=cex.axis)
mtext(expression("Time"), 1, line=2.75, cex=cex)
axis(2, cex.axis=cex.axis, las=1)
mtext("Number of molecules", 2, line=3.5, cex=cex)
box()
if(type=="p") {
	legend("bottomright", legend=species, y.intersp=1.25, pch=pch, pt.cex=cex, cex=1.5, col=colors)
} else {
	legend("bottomright", legend=species, y.intersp=1.25, lty=lty, lwd=lwd, cex=1.5, col=colors)
}

dev.off()

