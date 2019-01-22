## create_clhs_design.R
# Author: Nathan Owen
# Last updated: 14/01/2019

## Clear workspace and load 'clhs' package
rm(list = ls())
require(clhs)

# First generate a N X 4 data.frame of 4 mixture variables
# These will act as candidate points for the clhs library
N = 100000
X = matrix(NA, N, 4)
for (i in 1:N) {
  running_sum = 0
  for (j in 1:3) {
    # Generate a random number between 0 and (1 - running_sum)
    X[i, j] = runif(1, min = 0, max = 1 - running_sum)
    # Update running_sum
    running_sum = sum(X[i,], na.rm = TRUE)
  }
  # Set final variable as (1 - running_sum)
  X[i,4] = 1 - running_sum
  # Permute order of the variables
  X[i,] = X[i,][sample(1:4)]
}
# Convert to data.frame with column names (x1,x2,x3,x4)
X = data.frame(X)
colnames(X) = paste("x",1:4,sep="")

## Use 'clhs' package to sample a design of size n from candidate points in X
n = 40
rows = clhs(X, size = n)
clhs_design = X[rows,]

## Plot design
pairs(clhs_design, upper.panel = NULL, pch = 19)

## Save design to .txt file
write.table(clhs_design, file ="clhs.txt", quote = FALSE, row.names = FALSE)
