#!/usr/bin/env Rscript

####
#### Computes MLEs for kinetic parameters given model and complete data
####

args.vec <- commandArgs(TRUE)
if(length(args.vec) < 3) { stop("Usage: MLE_complete_data.R SBMLFILE(character) DATAFILE(character) OUTFILE(character)") }
sbmlfile <- args.vec[1]
datafile <- args.vec[2]
outfile <- args.vec[3]

### Read in model
require(SBMLR)
model <- readSBML(sbmlfile)

### Store model attributes
params <- model$globalParameters
for(i in seq_along(params)) {
	assign(names(params[i]), 1.0)
}
model.sum <- summary(model)
n <- model.sum$nStates
m <- model.sum$nReactions
V <- model.sum$incid

### Test for degenerate reaction stoichiometries
if(!identical(V, unique(V, MARGIN=2))) {
    stop("Degenerate reaction stoichiometries")
}

### Read in data
data <- as.matrix(read.table(datafile, header=TRUE, sep=" ", row.names=NULL, as.is=TRUE))

### Identify which reactions occur
changes <- diff(data)
reactions <- apply(changes[,rownames(V)], 1, function(w) apply(V, 2, function(v) identical(w,v)))

### Make sure one reaction occurred at each time point
if(!any(reactions[,ncol(reactions)])) {
    if(sum(reactions) != (ncol(reactions)-1)) stop("Each time point must coincide with one reaction firing")
} else {
    if(sum(reactions) != ncol(reactions)) stop("Each time point must coincide with one reaction firing")
}

### Sum number of reaction firings
R <- rowSums(reactions)

### Sum time-weighted hazard functions
H <- numeric(m)
for(i in 1:nrow(changes)) {
    H <- H + sapply(model$reactions, function(rxn) rxn$law(data[i,])) * changes[i,"Time"]
}

### Compute MLEs
MLE <- R/H

### Output results
results <- rbind(R, H, MLE)
colnames(results) <- names(params)
write.table(round(results,6), file=outfile, sep="\t", quote=FALSE)

