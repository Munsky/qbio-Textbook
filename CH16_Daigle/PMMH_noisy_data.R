#!/usr/bin/env Rscript

####
#### Estimates parameters from noisy data using PMMH method
####

args.vec <- commandArgs(TRUE)
if(length(args.vec) < 7) { stop("Usage: PMMH_noisy_data.R SBMLFILE(character) SIGMA(numeric) PRIOR(numeric) DATAFILE(character) NUMITER(numeric) NUMTHIN(numeric) TABLEFILE(character) [SEED](numeric)") }
sbmlfile <- args.vec[1]
sigma <- as.numeric(args.vec[2])
prior <- as.numeric(args.vec[3])
datafile <- args.vec[4]
numiter <- as.numeric(args.vec[5])
numthin <- as.numeric(args.vec[6])
tablefile <- args.vec[7]
if(length(args.vec) >= 8) {
	seed <- as.numeric(args.vec[8])
} else {
	seed <- 101
}

### Load required packages
require(SBMLR)
require(smfsb)

### Read in and summarize model
sbml <- readSBML(sbmlfile)
sbml.sum <- summary(sbml)

### Store model attributes
reactions <- sbml$reactions
th <- sbml.sum$globalVec
p <- length(th)
n <- sbml.sum$nStates
species <- sbml.sum$sIDs
Vt <- t(sbml.sum$incid)

### Create smfsb model
h <- function(x, t, th=th) {
	with(as.list(c(x, th)), {
		return(unname(sapply(reactions, function(rxn) eval(parse(text=tail(rxn$law, 2)[1])))))
	})
}
model <- list("Pre" = matrix(0,nrow(Vt),ncol(Vt)), "Post" = unname(Vt), "h" = h)

### Read in and reformat data
data <- as.matrix(read.table(datafile, header=TRUE, sep=" ", row.names=NULL, as.is=TRUE))
time <- data[, "Time"]
data <- data[, -c(1:2)]
rownames(data) <- time

### Define data likelihood function
dataLik <- function(x, t, y, log=TRUE, ...) {
	ll <- sum(dnorm(y, x, sigma, log=TRUE))
	return(ifelse(log, ll, exp(ll)))
}

### Define sampler for prior on the initial state
simx0 <- function(N, t0, ...) {
	mat <- matrix(rpois(N*n, prior), nrow=N, ncol=n)
	colnames(mat) <- species
	return(mat)
}

### Create particle filter-based marginal log-likelihood function
set.seed(seed)
mLLik <- pfMLLik(100, simx0, 0, StepGillespie(model), dataLik, data)

### Run PMMH
tune <- .05
ll <- -1e99
thmat <- matrix(0, nrow=numiter, ncol=p)
colnames(thmat) <- names(th)
for(i in seq(numiter)) {
	message(paste(i,""), appendLF=FALSE)
	for(j in seq(numthin)) {
		thprop <- th * exp(rnorm(p,0,tune))
		llprop <- mLLik(thprop)
		if(log(runif(1)) < llprop-ll) {
			th <- thprop
			ll <- llprop
		}
	}
	thmat[i, ] <- th
}
message("Done!")

### Reformat and output results
results <- cbind(seq(from=numthin-1, by=numthin, length.out=numiter), thmat, log(thmat))
colnames(results) <- c("Iter", colnames(thmat), paste0("log(", colnames(thmat), ")"))
write.table(results, file=tablefile, sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)

