## build_emulator.R
# Author: Nathan Owen
# Last updated: 14/01/2019

setwd("C:/Users/neo204/OneDrive - University of Exeter/NEVO/Rushton Model/Paper/Tidy Code/")

## Clear workspace and load DiceKriging and mixexp packages
rm(list = ls())
require(DiceKriging)
require(mixexp)

## Source loo_km.R and transform_emulator.R functions
source("loo_km.R")
source("transform_emulator.R")

## Load inputs (experimental design) and outputs (model runs)
inputs = read.table("clhs_design.txt",header=TRUE)
outputs = read.table("model_runs.txt")

## Define useful quantities
n = nrow(inputs)      # Size of training design
m = ncol(outputs)     # Length of time series output

## Use principal component analysis on the logarithm of the output (centered by mean) to create output variables
PCA = prcomp(log(outputs), center=TRUE)

## Choose k principal components to explain proportion p = 0.99 of variance and store in outputPCA
summary(PCA)$importance
variance_explained = summary(PCA)$importance[3,]
p = 0.99
k = min(which(variance_explained > p))
outputsPCA = matrix(PCA$x[,1:k], ncol=k)

## Create data frame for fitting models to (design + PC outputs)
dat = cbind(inputs, outputsPCA)
colnames(dat)[-(1:4)] = paste("PC",1:k,sep="")

## Fit a mixture model to each PC output independently. Store in a list structure.
## model = 1 assumes a linear mixture model (first order terms - constant)
mixmod = list()
for (i in 1:k) {
  mixmod[[i]] = MixModel(frame = dat, 
                         response = paste("PC",i,sep=""), 
                         mixcomps = paste("x",1:4,sep=""), 
                         model = 1)
}

## Fit Gaussian process emulators with covariance 'ctype' to each PC output independently
# Use the mixture model as a mean function and assume coefficients known
# Estimate other parameters as normal and fix nugget = 1e-7
ctype = "gauss"
#ctype = "matern5_2"
#ctype = "matern3_2"
#ctype = "exp"

gpmod = list()
for (i in 1:k) {
  gpmod[[i]] = km(formula = formula(mixmod[[i]]),
                  design = inputs,
                  response = dat[, ncol(inputs)+i],
                  covtype = ctype,
                  nugget = 1e-7,
                  coef.trend = coefficients(mixmod[[i]]))
}
gpmod

# Look at log-likelihood for each Gaussian process
gpmod[[1]]@logLik
gpmod[[2]]@logLik

## Leave-one-out validation
LOOgpmoda = loo_km(gpmod[[1]], plotyes=FALSE)
LOOgpmodb = loo_km(gpmod[[2]], plotyes=FALSE)
pLOO = list(LOOgpmoda,LOOgpmodb)

## Plot: leave-one-out validation emulator mean and uncertainty against principal component output
par(mfrow = c(1,2), mar = c(4,4,1,1), ps = 12)

plot(1, type = "n", xlab = "", ylab = "", xlim = c(min(LOOgpmoda$lower95), max(LOOgpmoda$upper95)), ylim = c(min(LOOgpmoda$lower95), max(LOOgpmoda$upper95)))
mtext(expression(tilde(y)[1]), side = 1, line = 2.5)
mtext(expression(hat(eta)[1](bold(x))), side = 2, line = 2.5)
segments(outputsPCA[,1], LOOgpmoda$lower95, outputsPCA[,1], LOOgpmoda$upper95, col = "red", lwd = 2)
points(outputsPCA[,1], LOOgpmoda$mean, pch = 19, cex = 1)
abline(0, 1)

plot(1, type = "n", xlab = "", ylab = "", xlim = c(min(LOOgpmodb$lower95), max(LOOgpmodb$upper95)), ylim = c(min(LOOgpmodb$lower95), max(LOOgpmodb$upper95)))
mtext(expression(tilde(y)[2]), side = 1, line = 2.5)
mtext(expression(hat(eta)[2](bold(x))), side = 2, line = 2.5)
segments(outputsPCA[,2], LOOgpmodb$lower95, outputsPCA[,2], LOOgpmodb$upper95, col = "red", lwd = 2)
points(outputsPCA[,2], LOOgpmodb$mean, pch = 19, cex = 1)
abline(0, 1)

# Coverage
100*(1 - mean((outputsPCA[,1] < LOOgpmoda$lower95) + (outputsPCA[,1] > LOOgpmoda$upper95)))
100*(1 - mean((outputsPCA[,2] < LOOgpmodb$lower95) + (outputsPCA[,2] > LOOgpmodb$upper95)))

## Transform emulator predictions to scale of original data
pLOOT = transform_emulator(pLOO, PCA)

## Calculate error metrics and LOO points
# Root mean square error (RMSE)
RMSE_LOO = sqrt(rowMeans((exp(pLOOT$mean) - outputs)^2))
min(RMSE_LOO)
max(RMSE_LOO)
mean(RMSE_LOO)

# Maximum absolute error
MAE_LOO = apply(abs(exp(pLOOT$mean) - outputs), 1, max)
min(MAE_LOO)
max(MAE_LOO)
mean(MAE_LOO)

# Plot: (model time series - emulator time series) with minimum and maximum RMSE metric for LOO CV points
par(mfrow = c(2,1), mar = c(4,4,1,1), ps = 12)

i = which.min(RMSE_LOO)
min(RMSE_LOO)
inputs[i,]
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(1, m), ylim = c(min(outputs - exp(pLOOT$upper95)), max(outputs-exp(pLOOT$lower95))))
axis(side = 1, at = c(0, 1000, 2000, 3000), labels = TRUE)
mtext("time [day]", side = 1, line = 2.5)
axis(side = 2, at = c(-2, -1.5, -1, -0.5, 0, 0.5), labels = TRUE)
mtext(expression(eta(bold(x)) - hat(eta)(bold(x)) ~ "[" ~ m^3/s ~ "]"), side = 2, line = 2.5)
polygon(c(1:m, rev(1:m)), c(outputs[i,] - exp(pLOOT$lower95[i,]), rev(outputs[i,] - exp(pLOOT$upper95[i,]))), density = NA, col = "grey80")
lines(1:m, outputs[i,] - exp(pLOOT$mean[i,]), lwd=1)
legend("bottomleft", legend = c(expression(RMSE == 1.127 %*% 10^-2), expression((list(0.778,0.028,0.109,0.085)))), bty = "n", y.intersp = 1.2)
box(lwd = 1)

i = which.max(RMSE_LOO)
inputs[i,]
max(RMSE_LOO)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(1, m), ylim = c(min(outputs - exp(pLOOT$upper95)), max(outputs - exp(pLOOT$lower95))))
axis(side = 1, at = c(0, 1000, 2000, 3000), labels = TRUE)
mtext("time [day]", side = 1, line = 2.5)
axis(side = 2, at = c(-2, -1.5, -1, -0.5, 0, 0.5), labels = TRUE)
mtext(expression(eta(bold(x)) - hat(eta)(bold(x)) ~ "[" ~ m^3/s ~ "]"), side = 2, line = 2.5)
polygon(c(1:m, rev(1:m)), c(outputs[i,] - exp(pLOOT$lower95[i,]), rev(outputs[i,] - exp(pLOOT$upper95[i,]))), density = NA, col = "grey80")
lines(1:m, outputs[i,] - exp(pLOOT$mean[i,]), lwd = 1)
legend("bottomleft", legend=c(expression(RMSE == 4.387 %*% 10^-1), expression((list(0, 0.001, 0.017, 0.982)))), bty = "n", y.intersp = 1.2)
box(lwd = 1)

