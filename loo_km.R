loo_km <- function(kmmodel, plotyes = TRUE) {
  ## loo_km.R
  # Author: Nathan Owen
  # Last updated: 14/01/2019
  
  # Extract info from km object
  n = kmmodel@n
  trend.formula = kmmodel@trend.formula
  X = kmmodel@X
  y = kmmodel@y
  name = kmmodel@covariance@name
  trend.coef = kmmodel@trend.coef
  range.val = kmmodel@covariance@range.val
  sd2 = kmmodel@covariance@sd2
  Nugget = kmmodel@covariance@nugget
  d = kmmodel@p
  
  # Pre-allocate vectors to store results
  means = numeric(n)
  sds = numeric(n)
  l95 = numeric(n)
  u95 = numeric(n)
  
  # Perform leave-one-out cross validation for each design point
  for (i in 1:n) {
    # Refit emulator with design point i withheld (parameters fixed)
    tempmod = km(formula = trend.formula,
                 design = data.frame(matrix(X[-i,], nrow = n - 1, dimnames = list(1:(n - 1), paste("x",1:d,sep="")))),
                 response=y[-i],
                 covtype = name,
                 coef.trend = trend.coef,
                 coef.cov = range.val,
                 coef.var = sd2,
                 nugget = Nugget)
    # Predict from this emulator at design point i
    ptempmod = predict(tempmod,
                       newdata = data.frame(matrix(X[i,], nrow = 1, dimnames = list(1, paste("x", 1:d, sep="")))),
                       type="UK")
    # Store results in i-th element of vectors
    means[i] = ptempmod$mean
    sds[i] = ptempmod$sd
    l95[i] = ptempmod$lower95
    u95[i] = ptempmod$upper95
  }
  
  # Calculate coverage
  coverage = 100 * (1 - mean((y < l95) + (y > u95)))
  
  # If requested, create plot
  if (isTRUE(plotyes)) {
    plot(1, type = "n", xlab = "Simulator", ylab = "Emulator", xlim = c(min(l95), max(u95)), ylim = c(min(l95), max(u95)), main = bquote(paste("coverage: ", .(round(coverage,2)), "%")))
    segments(y, l95, y, u95, col="red")
    points(y, means, pch=19)
    abline(0, 1)
  }
  
  # Print coverage and return results to user
  print(paste("coverage: ", coverage, "%", sep=""))
  return(list(mean = means, sd = sds, lower95 = l95, upper95 = u95))
  
}