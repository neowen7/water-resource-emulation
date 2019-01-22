transform_emulator <- function(pmod, PCA) {
  ## transform_emulator.R
  # Author: Nathan Owen
  # Last updated: 14/01/2019
  
  ## Function to take the predictions from k independent emulators of k principal components
  ## and convert back to original data space
  
  # pmod must be a list, with the predictions from the ith emulator stored in pmod[[i]]
  if (!is.list(pmod)) stop("pmod must be a list")
  
  if (length(names(pmod[[1]])) == 7 || length(names(pmod[[1]])) == 4) {
    # Define dimensionality
    k = length(pmod)            # number of principal components used
    N = length(pmod[[1]]$mean)  # number of prediction points
    M = length(PCA$center)      # length of time series output
    n = ncol(PCA$rotation)      # size of training design
    
    # First store the k prediction means in a matrix of size Nxk
    meanz = matrix(NA,N,k)
    for (i in 1:k) {
      meanz[,i] = pmod[[i]]$mean
    }
    
    # Transform the means back to data space using first k principal components and add center
    meany = scale(meanz%*%t(PCA$rotation[,1:k]),center=-PCA$center,scale=FALSE)
    
    # Now do uncertainty
    # Emulator uncertainty
    varGP = matrix(NA,N,M)
    for (i in 1:N) {
      if (k == 1) {
        varGP[i,] = diag(PCA$rotation[,k]%*%as.matrix(pmod[[k]]$sd[i]^2)%*%t(PCA$rotation[,k]))
      } else if (k > 1) {
        varz = numeric(k)
        for (j in 1:k) {
          varz[j] = pmod[[j]]$sd[i]^2
        }
        varGP[i,] = diag(PCA$rotation[,1:k]%*%diag(varz)%*%t(PCA$rotation[,1:k]))
      }
    }
    # Principal component reconstruction uncertainty
    varPC = diag(PCA$rotation[,-(1:k)]%*%diag(PCA$sdev[-(1:k)]^2)%*%t(PCA$rotation[,-(1:k)]))

    # Construct 95% intervals
    qtt = qt(0.975,n-4)
    lgp = meany-qtt*sqrt(varGP)
    ugp = meany+qtt*sqrt(varGP)
    
    # Construct 2 sd intervals
    lpc = meany-qtt*sqrt(varPC)
    upc = meany+qtt*sqrt(varPC)
    
    # Construct 2 sd intervals
    l = meany-qtt*sqrt(varGP+varPC)
    u = meany+qtt*sqrt(varGP+varPC)
    
    return(list(mean=meany,lower95=l,upper95=u,lowerGP95=lgp,upperGP95=ugp,lowerPC95=lpc,upperPC95=upc))
  } else {
    k = length(pmod)            # number of principal components used
    N = length(pmod[[1]])       # number of prediction points
    M = length(PCA$center)      # length of time series output
    
    # First store the k prediction means in a matrix of size Nxk
    meanz = matrix(NA,N,k)
    for (i in 1:k) {
      meanz[,i] = pmod[[i]]
    }
    
    # Transform the means back to data space using first k principal components and add center
    meany = scale(meanz%*%t(PCA$rotation[,1:k]),center=-PCA$center,scale=FALSE)
    
    return(meany)
  }
}