#----------------------------------------------------------------------------------------------------------------
# Case Study 2: Hospital admissions due to respiratory related health problems in Greater Glasgow & Clyde
#----------------------------------------------------------------------------------------------------------------
# R libraries
library(shapefiles)
library(sp)
library(spdep)
library(CARBayes)
library(CARBayesST)
library(Rcpp)
library(truncdist)
library(coda)
library(boot)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(colorRamps)
library(LearnBayes)
library(MCMCpack)
library(gtools)
library(grid)
library(MASS)

#----------------------------------------------------------------------------------------------------------------
# Read in the Case Study 2 data
observed <- read.csv(file="respiratory admissions.csv")
expected <- read.csv(file="respiratory admissions expected counts.csv")
rownames(observed) <- rownames(expected) <- observed[, 1]

# The Intermediate Zones for Greater Glasgow & Clyde
Glasgow.IZs <- read.csv(file = "GlasgowIZ.csv")
Glasgow.IZs <- as.character(Glasgow.IZs[, 1])

# Number of hospital admissions in Greater Glasgow & Clyde for each year
Ymat <- as.matrix(observed[Glasgow.IZs, -1])

# Expected number of hospital admissions in Greater Glasgow & Clyde for each year
Emat <- as.matrix(expected[Glasgow.IZs, -1])

# Number of intermediate zones
K <- nrow(Ymat)
# Number of time points (years)
N <- ncol(Ymat)

# Read in the shapefiles for Scotland
shp <- read.shp(shp.name="Scotland study area.shp")
dbf <- read.dbf(dbf.name="Scotland study area.dbf")

# Create the combined spatial data object
spdat <- combine.data.shapefile(Ymat/Emat, shp, dbf)

# Create the neighbourhood matrix W
W.nb <- poly2nb(spdat, row.names = rownames(Ymat))
W.list <- nb2listw(W.nb, style = "B", zero.policy = TRUE)
W <- nb2mat(W.nb, style = "B", zero.policy = TRUE)

#----------------------------------------------------------------------------------------------------------------
# CARBayesST reads in the data as vectors where the first K elements correspond to the
# the observations in the first time point for all K intermediate zones, and the
# next K elements correspond to the second time point etc.

Y <- as.numeric(Ymat)
E <- as.numeric(Emat)

# Model formula
formula <- Y ~ offset(log(E))

# MCMC setup
burnin <- 200000  
n.sample <- 300000
thin <- 10
Nchains <- 4

# Read in the functions for the mixture of trends model (Poisson version)
sourceCpp("Poisson.LTM.General.cpp")
source('Poisson.LTM.General.R')

# The model function requires the following arguments
# formula: the model formula
# W: the neighbourhood matrix
# burnin: the burn-in period of the MCMC
# n.sample: the number of iterations of the MCMC
# thin: the level of thinning of the MCMC
# Nchains: The number of MCMC chains used in the Metropolis coupling
# trends: A character vector containing the selected temporal trends, which are
### Constant - Constant trend
### LI - linear increasing trend
### LD - linear decreasing trend
### LD - linear decreasing trend
### CP - changepoint trend (linear increase to 'changepoint' followed by linear decrease)
### CT - changepoint trend (linear decrease to 'changepoint' followed by linear increase)
### MD - monotonically decreasing trend (requires number of 'knots' to be selected)
### MI - monotonically increasing trend (requires number of 'knots' to be selected)
# changepoint: required if CP or CT trends selected
# knots: the number of knots for the monotonic trends (MD/MI)

# Here the constant, linear increasing, and linear decreasing trends are selected

trends <- c("Constant", "LI", "LD")
changepoint <- NULL
knots <- NULL

set.seed(1)

model <- poisson.MixMod(formula, data=NULL, W=W, burnin=burnin, n.sample=n.sample, thin=thin, trends=trends,
                        changepoint=changepoint, knots=knots, prior.mean.beta=NULL, prior.var.beta=NULL, 
                        prior.mean.gamma=NULL, prior.var.gamma=NULL, prior.lambda=NULL, prior.tau2=NULL,
                        Nchains=Nchains, verbose=TRUE)

#----------------------------------------------------------------------------------------------------------------
# Some output from the model
# Model summary
model$summary.results

# To check trace plots
plot(model$samples$beta)
plot(model$samples$lambda)
plot(model$samples$tau2)
plot(model$samples$rho)

# For checking the gamma parameters corresponding to the trends
plot(model$samples$gamma[2, ], type = "l") # gamma.LD
plot(model$samples$gamma[3, ], type = "l") # gamma.LI
# model$samples$gamma[1, ] correspond to the constant trend and so is just zero. Due to how the gamma parameters
# are brought together. This output will be improved onced incorporated into the CARBayesST package

# Checking trend allocations
# The trend allocations are in the order the are given in the summary table,
# i.e. 1: constant, 2: linear decreasing, 3: linear increasing
model$trends
table(model$trends[, 2])

# Probabilities associated with each trend for each IZ
model$trend.probs

#----------------------------------------------------------------------------------------------------------------
# Plot of posterior probabilities that each IZ is assigned to each trend (Figure 5)
trend.p <- model$trend.probs

# Reordering trends to Decreasing, Constant, Increasing
swap1 <- trend.p[, 1]
swap2 <- trend.p[, 2]
trend.p[, 1] <- swap2
trend.p[, 2] <- swap1

trends.mod <- model$trends[, 2]
swap1 <- which(trends.mod == 1)
swap2 <- which(trends.mod == 2)
trends.mod[swap1] <- 2
trends.mod[swap2] <- 1

plot.probs <- cbind(trend.p, trends.mod)
# Ordering the probabilites by trend
plot.probs <- plot.probs[order(plot.probs[, 4]), ]
# Number of areas per trend
per.trend <- table(plot.probs[, 4])

plot(1:K, rep(1, K), ylim = c(0, 1), type = "n", axes = FALSE, ylab = "", xlab = "")

for(i in 1:K)
{
  pr <- plot.probs[i, -4]
  lines(c(i, i), c(0, pr[1]), col ="blue", lwd = 3)
  lines(c(i, i), c(pr[1], sum(pr[1:2])), col ="black", lwd = 3)
  lines(c(i, i), c(sum(pr[1:2]), 1), col ="red", lwd = 3)
}

lines(c(per.trend[1] + 0.5, per.trend[1] + 0.5), c(0, 1), lwd = 3, col = "white")
lines(c(sum(per.trend[1:2]) + 0.5, sum(per.trend[1:2]) + 0.5), c(0, 1), lwd = 3, col = "white")

labs <- c(per.trend[1]/2, per.trend[1]/2 + sum(per.trend[1:2])/2, sum(per.trend[1:2])/2 + sum(per.trend)/2)

axis(1, at = labs, line = -1, tick = FALSE, labels = c("Decreasing", "Constant", "Increasing"), cex.axis = 1.5)
axis(2, at = seq(from = 0, to = 1, by = 0.1), line = -1.8, cex.axis = 1.1)
title(ylab = "Probability", line=1.1, cex.lab=1.5)
title(xlab = "Temporal trend classification", line=2.4, cex.lab=1.75)

#----------------------------------------------------------------------------------------------------------------
# Map of the area classification including measure of uncertainty from the probabilites (Figure 5)
trends.all <- rep(NA, K)
for (k in 1:K) {
  trends.all[k] <- trend.p[k, trends.mod[k]] + trends.mod[k]
}

# convert to ranges 0.33-0.5, 0.5-0.75 and 0.75-1
new.trends.all <- rep(NA, K)

for(i in 1:K)
{
  # decreasing trend
  if(trends.all[i] >= 1.33 & trends.all[i] <= 1.5)
  {
    new.trends.all[i] <- 1
  }
  if(trends.all[i] >= 1.5 & trends.all[i] <= 1.75)
  {
    new.trends.all[i] <- 2
  }
  if(trends.all[i] >= 1.75 & trends.all[i] <= 2)
  {
    new.trends.all[i] <- 3
  }
  # constant trend
  if(trends.all[i] >= 2.33 & trends.all[i] <= 2.5)
  {
    new.trends.all[i] <- 4
  }
  if(trends.all[i] >= 2.5 & trends.all[i] <= 2.75)
  {
    new.trends.all[i] <- 5
  }
  if(trends.all[i] >= 2.75 & trends.all[i] <= 3)
  {
    new.trends.all[i] <- 6
  }
  # increasing trend
  if(trends.all[i] >= 3.33 & trends.all[i] <= 3.5)
  {
    new.trends.all[i] <- 7
  }
  if(trends.all[i] >= 3.5 & trends.all[i] <= 3.75)
  {
    new.trends.all[i] <- 8
  }
  if(trends.all[i] >= 3.75 & trends.all[i] <= 4)
  {
    new.trends.all[i] <- 9
  }
}

labels.place <- 1:9 + 0.5
labels.range <- rep(c("0.33 - 0.5", "0.5 - 0.75", "0.75 - 1"), 3)

spdat@data$trends.all <- new.trends.all

plot.lay <- c(1, 1)
l3 <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(223000, 647500), scale = 5000)

main.plot <- spplot(spdat, c("trends.all"), col = "transparent", layout = plot.lay, sp.layout = list(l3),
                    colorkey = list(labels=list(at = labels.place, labels = labels.range, cex = 1)),
                    at = seq(from = 1, to = 10, length = 10),
                    col.regions = c(hsv(0.6, seq(0.15, 1, length.out = 33), 1), gray(seq(0.95, 0.15, length.out = 33), 1), hsv(0.05, seq(0.15, 1, length.out = 33), 1)))

print(main.plot, position = c(0, 0, 1, 1))
grid.text("Decreasing", x=unit(0.95, "npc"), y=unit(0.185, "npc"), rot=-90)
grid.text("Constant", x=unit(0.95, "npc"), y=unit(0.5, "npc"), rot=-90)
grid.text("Increasing", x=unit(0.95, "npc"), y=unit(0.81, "npc"), rot=-90)

#----------------------------------------------------------------------------------------------------------------
# Plots of the linear trends (Figure 4)
Const.trend.median <- rep(exp(median(model$samples$beta)), N)
Const.trend.lower <- rep(exp(quantile(model$samples$beta, prob = 0.025)), N)
Const.trend.upper <- rep(exp(quantile(model$samples$beta, prob = 0.975)), N)

Dec.trend.median <- exp(median(model$samples$beta) + (median(model$samples$gamma[2, ]) * 1:N))
Dec.trend.lower <- exp(quantile(model$samples$beta, prob = 0.025) + (quantile(model$samples$gamma[2, ], prob = 0.025) * 1:N))
Dec.trend.upper <- exp(quantile(model$samples$beta, prob = 0.975) + (quantile(model$samples$gamma[2, ], prob = 0.975) * 1:N))

Inc.trend.median <- exp(median(model$samples$beta) + (median(model$samples$gamma[3, ]) * 1:N))
Inc.trend.lower <- exp(quantile(model$samples$beta, prob = 0.025) + (quantile(model$samples$gamma[3, ], prob = 0.025) * 1:N))
Inc.trend.upper <- exp(quantile(model$samples$beta, prob = 0.975) + (quantile(model$samples$gamma[3, ], prob = 0.975) * 1:N))

trend.range <- range(c(Dec.trend.lower, Inc.trend.lower, Dec.trend.upper, Inc.trend.upper))
years <- colnames(Ymat)
years <- sub("Y", "", years)

plot(1:N, Inc.trend.median, ylim = trend.range, xlab = "Year", axes = FALSE, col = "red", type = 'n',
     ylab = "Risk of respiratory health problems")
axis(1, at = 1:N, labels = years)
axis(2)

lines(1:N, Const.trend.median, col = "black", lwd = 2)
lines(1:N, Const.trend.lower, col = "black", lty = 2, lwd = 2)
lines(1:N, Const.trend.upper, col = "black", lty = 2, lwd = 2)

lines(1:N, Inc.trend.median, col = "red", lwd = 2)
lines(1:N, Inc.trend.lower, col = "red", lty = 2, lwd = 2)
lines(1:N, Inc.trend.upper, col = "red", lty = 2, lwd = 2)

lines(1:N, Dec.trend.median, col = "blue", lwd = 2)
lines(1:N, Dec.trend.lower, col = "blue", lty = 2, lwd = 2)
lines(1:N, Dec.trend.upper, col = "blue", lty = 2, lwd = 2)

box()

legend(1, 0.65, title = "Temporal trend", c("Constant", "Decreasing", "Increasing"), lty = c(1, 1, 1),
       lwd = c(2, 2, 2), col = c("black", "blue", "red"), bty = "n")

#----------------------------------------------------------------------------------------------------------------