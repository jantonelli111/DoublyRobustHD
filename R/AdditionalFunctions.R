#################################################################
######################## MCMC functions #########################
#################################################################

SplineOutcomeMCMC = function(y, tMat, x, whichCat, df, type="continuous",
                             nBurn=1000, nScans=5000,
                             thin=4, nChains=2, a=0.001, b=0.001,
                             c=2, d=dim(x)[2], e=0.5, f=0.5) {

  ## creating design matrices for categorical variables

  if (is.vector(tMat) == TRUE) {
    dt = 1
  } else {
    dt = dim(tMat)[2]
  }

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    cols = list()
    for (j in 1 : pCont) {
      cols[[j]] = ((j-1)*df + 1):(j*df)
    }
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]

    lengthCat = c()
    cols = list()
    for (j in 1 : pCont) {
      cols[[j]] = ((j-1)*df + 1):(j*df)
    }

    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }

    for (j2 in 1 : length(whichCat)) {
      cols[[pCont + j2]] = colsCat[[j2]] + pCont*df
    }
  }

  for (j in 1 : p) {
    cols[[j]] = cols[[j]] + dt + 1
  }

  ## parameters and prior specification

  n = length(y)

  ## MCMC details
  dfY = df
  muY = rep(0, dfY)

  Group1VarY = 1

  PriorSigmaY = Group1VarY*diag(dfY)

  sigmaYPost = matrix(NA, nChains, nScans)
  gammaYPost = array(NA, dim=c(nChains, nScans, p))
  betaYPost = array(NA, dim=c(nChains, nScans, max(cols[[p]])))
  tauYPost = matrix(NA, nChains, nScans)
  sigmaBetaYPost = matrix(NA, nChains, nScans)

  sigmaYPost[,1] = 2
  gammaYPost[,1,] = rbinom(nChains*p, 1, p=0.1)
  for (nc in 1 : nChains) {
    betaYPost[nc,1,] = c(rnorm(1 + dt), rnorm(max(cols[[p]]) - dt - 1, sd=0.2)*
                           rep(gammaYPost[nc,1,], times=sapply(cols, length)))
  }
  tauYPost[,1] = 0.1
  sigmaBetaYPost[,1] = Group1VarY

  designY = cbind(rep(1,n), tMat)

  for (j in 1 : pCont) {
    tempY = scale(splines::ns(xCont[,j], dfY))
    designY = cbind(designY, tempY)
  }

  if (length(whichCat) > 0) {
    designY = cbind(designY, xCat2)
  }

  for (ni in 2 : nScans) {
    for (nc in 1 : nChains) {
      if (nc == 1 & ni %% 100 == 0) print(ni)

      ################# outcome model ##################

      ## Update sigma squared
      if (type == "continuous") {
        aStar = a + n/2 + gammaYPost[nc,ni-1,]*sapply(cols, length)/2
        bStar = b + sum((y - (designY %*% betaYPost[nc,ni-1,]))^2)/2 +
          sum(betaYPost[nc,ni-1,-c(1:(1+dt))]^2)/(2*Group1VarY)
        sigmaYPost[nc,ni] = 1/rgamma(1,aStar,bStar)
        Zy = y
      } else if (type == "binary") {
        sigmaYPost[nc,ni] = 1
        Zy = rep(NA, n)

        meanZy = designY %*% betaYPost[nc,ni-1,]

        Zy[y==1] = truncnorm::rtruncnorm(sum(y==1), a=0, mean = meanZy[y==1], sd=1)
        Zy[y==0] = truncnorm::rtruncnorm(sum(y==0), b=0, mean = meanZy[y==0], sd=1)
      }

      ## Update sigmaBeta
      sigmaBetaYPost[nc,ni] = 1/rgamma(1, e + sum(gammaYPost[nc,ni-1,]*sapply(cols, length))/2,
                                       f + sum(betaYPost[nc,ni-1,-c(1:(1+dt))]^2)/(2*sigmaYPost[nc,ni]))
      Group1VarY = sigmaBetaYPost[nc,ni]

      ## Update tau
      tauYPost[nc,ni] = rbeta(1, c + sum(gammaYPost[nc,ni-1,] == 1),
                              d + sum(gammaYPost[nc,ni-1,] == 0))

      ## Update regression coefficients and variable inclusion parameters
      tempBeta = betaYPost[nc,ni-1,]
      for (j in 1 : p) {
        tempCols = cols[[j]]

        PriorSigmaY = Group1VarY*sigmaYPost[nc,ni]*diag(length(tempCols))
        muY = rep(0, length(tempCols))

        yStar = Zy - designY[,-tempCols] %*% tempBeta[-tempCols]

        ## probability of being in group zero
        p0 = log(1 - tauYPost[nc,ni])

        ## probability of being in top group
        muVar = solve(t(designY[,tempCols]) %*% designY[,tempCols]  / sigmaYPost[nc,ni] +
                        solve(PriorSigmaY))
        muBeta = muVar %*% (t(designY[,tempCols]) %*% yStar/sigmaYPost[nc,ni] +
                              solve(PriorSigmaY) %*% muY)
        p1 = log(tauYPost[nc,ni]) + mvtnorm::dmvnorm(rep(0, length(tempCols)),
                                            mean=muY, sigma=PriorSigmaY, log=TRUE) -
          mvtnorm::dmvnorm(rep(0, length(tempCols)), mean=muBeta, sigma=muVar, log=TRUE)

        maxlog = max(p0,p1)

        p0new = exp(-maxlog + p0)
        p1new = exp(-maxlog + p1)

        gammaYPost[nc,ni,j] = sample(0:1, size=1, p=c(p0new,p1new))

        tempBeta[tempCols] = rep(0, length(tempCols))
        if (gammaYPost[nc,ni,j] == 1) tempBeta[tempCols] = mvtnorm::rmvnorm(1, muBeta, sigma=muVar)
      }
      betaYPost[nc,ni,] = tempBeta

      ## Update intercept and treatment effect
      yStar = Zy - designY[,-c(1:(1+dt))] %*% betaYPost[nc,ni,-c(1:(1+dt))]
      tempDesign = designY[,1:(1+dt)]
      betaYPost[nc,ni,c(1:(1+dt))] = mvtnorm::rmvnorm(1, solve(t(tempDesign) %*% tempDesign)
                                             %*% t(tempDesign) %*% yStar,
                                             sigma=sigmaYPost[nc,ni]*solve(t(tempDesign) %*% tempDesign))

    }
  }

  keep = seq(nBurn + 1, nScans, by=thin)

  l = list(beta = betaYPost[,keep,],
           gamma = gammaYPost[,keep,],
           sigma = sigmaYPost[,keep],
           sigmaBeta = sigmaBetaYPost[,keep])

  return(l)

}






SplineTreatmentMCMC = function(t, x, df,  whichCat, type="continuous",
                               nBurn=1000, nScans=5000,
                               thin=4, nChains=2, a=0.001, b=0.001,
                               c=2, d=dim(x)[2], e=0.5, f=0.5) {

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    cols = list()
    for (j in 1 : pCont) {
      cols[[j]] = ((j-1)*df + 1):(j*df)
    }
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]

    lengthCat = c()
    cols = list()
    for (j in 1 : pCont) {
      cols[[j]] = ((j-1)*df + 1):(j*df)
    }

    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }

    for (j2 in 1 : length(whichCat)) {
      cols[[pCont + j2]] = colsCat[[j2]] + pCont*df
    }
  }

  for (j in 1 : p) {
    cols[[j]] = cols[[j]] + 1
  }

  ## parameters and prior specification

  n = length(y)

  ## MCMC details
  dfT = df
  muT = rep(0, dfT)

  Group1VarT = 1

  sigmaTPost = matrix(NA, nChains, nScans)
  gammaTPost = array(NA, dim=c(nChains, nScans, p))
  betaTPost = array(NA, dim=c(nChains, nScans, max(cols[[p]])))
  tauTPost = matrix(NA, nChains, nScans)
  sigmaBetaTPost = matrix(NA, nChains, nScans)

  sigmaTPost[,1] = 2
  gammaTPost[,1,] = rbinom(nChains*p, 1, p=0.1)
  for (nc in 1 : nChains) {
    betaTPost[nc,1,] = c(rnorm(1), rnorm(max(cols[[p]]) - 1, sd=0.2)*
                           rep(gammaTPost[nc,1,], times=sapply(cols, length)))
  }
  tauTPost[,1] = 0.1

  designT = cbind(rep(1,n))

  for (j in 1 : pCont) {
    tempT = scale(splines::ns(xCont[,j], dfT))
    designT = cbind(designT, tempT)
  }

  if (length(whichCat) > 0) {
    designT = cbind(designT, xCat2)
  }

  for (ni in 2 : nScans) {
    for (nc in 1 : nChains) {
      if (nc == 1 & ni %% 100 == 0) print(ni)

      ## Update sigma squared
      if (type == "continuous") {
        aStar = a + n/2 + gammaTPost[nc,ni-1,]*sapply(cols, length)/2
        bStar = b + sum((t - (designT %*% betaTPost[nc,ni-1,]))^2)/2 +
          sum(betaTPost[nc,ni-1,-c(1)]^2)/(2*Group1VarT)
        sigmaTPost[nc,ni] = 1/rgamma(1,aStar,bStar)
        Zt = t
      } else if (type == "binary") {
        sigmaTPost[nc,ni] = 1
        Zt = rep(NA, n)

        meanZt = designT %*% betaTPost[nc,ni-1,]

        Zt[t==1] = truncnorm::rtruncnorm(sum(t==1), a=0, mean = meanZt[t==1], sd=1)
        Zt[t==0] = truncnorm::rtruncnorm(sum(t==0), b=0, mean = meanZt[t==0], sd=1)
      }

      ## Update sigmaBeta
      sigmaBetaTPost[nc,ni] = 1/rgamma(1, e + sum(gammaTPost[nc,ni-1,]*sapply(cols, length))/2,
                                       f + sum(betaTPost[nc,ni-1,-c(1)]^2)/(2*sigmaTPost[nc,ni]))
      Group1VarT = sigmaBetaTPost[nc,ni]

      ## Update tau
      tauTPost[nc,ni] = rbeta(1, c + sum(gammaTPost[nc,ni-1,] == 1),
                              d + sum(gammaTPost[nc,ni-1,] == 0))

      ## Update regression coefficients and variable inclusion parameters
      tempBeta = betaTPost[nc,ni-1,]
      for (j in 1 : p) {
        tempCols = cols[[j]]

        PriorSigmaT = Group1VarT*sigmaTPost[nc,ni]*diag(length(tempCols))
        muT = rep(0, length(tempCols))

        tStar = Zt - designT[,-tempCols] %*% tempBeta[-tempCols]

        ## probability of being in group zero
        p0 = log(1 - tauTPost[nc,ni])

        ## probability of being in top group
        muVar = solve(t(designT[,tempCols]) %*% designT[,tempCols]  / sigmaTPost[nc,ni] +
                        solve(PriorSigmaT))
        muBeta = muVar %*% (t(designT[,tempCols]) %*% tStar/sigmaTPost[nc,ni] +
                              solve(PriorSigmaT) %*% muT)
        p1 = log(tauTPost[nc,ni]) + mvtnorm::dmvnorm(rep(0, length(tempCols)), mean=muT,
                                            sigma=PriorSigmaT, log=TRUE) -
          mvtnorm::dmvnorm(rep(0, length(tempCols)), mean=muBeta, sigma=muVar, log=TRUE)

        maxlog = max(p0,p1)

        p0new = exp(-maxlog + p0)
        p1new = exp(-maxlog + p1)

        gammaTPost[nc,ni,j] = sample(0:1, size=1, p=c(p0new,p1new))

        tempBeta[tempCols] = 0
        if (gammaTPost[nc,ni,j] == 1) tempBeta[tempCols] = mvtnorm::rmvnorm(1, muBeta, sigma=muVar)
      }
      betaTPost[nc,ni,] = tempBeta

      ## Update intercept
      tStar = Zt - designT[,-1] %*% betaTPost[nc,ni,-1]
      betaTPost[nc,ni,1] = rnorm(1, mean(tStar), sd=sqrt(sigmaTPost[nc,ni]/n))
    }
  }

  keep = seq(nBurn + 1, nScans, by=thin)

  l = list(beta = betaTPost[,keep,],
           gamma = gammaTPost[,keep,],
           sigma = sigmaTPost[,keep],
           sigmaBeta = sigmaBetaTPost[,keep])

  return(l)

}








#################################################################
######################## MCMC functions #########################
#################################################################

GPOutcomeMCMC = function(y, tMat, x, band = 1, type="continuous",
                         whichCat,
                         nBurn=1000, nScans=5000,
                         thin=4, nChains=2, a=0.001, b=0.001,
                         c=2, d=dim(x)[2], e=0.5, f=0.5) {


  if (is.vector(tMat) == TRUE) {
    dt = 1
  } else {
    dt = dim(tMat)[2]
  }

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    pCat = 0
    nCatCols = 0
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]
    pCat = p - pCont

    lengthCat = c()
    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }

    for (j in 1 : pCat) {
      colsCat[[j]] = colsCat[[j]] + dt + 1
    }

    nCatCols = dim(xCat2)[2]
  }

  ## parameters and prior specification

  n = length(y)

  ## add noise to observations with input values too similar
  for (j in 1 : pCont) {
    if (length(unique(xCont[,j])) < n) xCont[,j] = xCont[,j] + rnorm(n, sd=.00001)
  }

  fYPost = array(NA, dim=c(nChains, nScans, n, pCont))
  sigmaYPost = array(NA, dim=c(nChains, nScans))
  gammaYPost = array(NA, dim=c(nChains, nScans, p))
  betaYPost = array(NA, dim=c(nChains, nScans, 1+dt+nCatCols))
  tauYPost = array(NA, dim=c(nChains, nScans, pCont))
  wYPost = array(NA, dim=c(nChains, nScans))
  sigmaBetaYPost = array(NA, dim=c(nChains, nScans))

  Group1VarY = 1
  sigmaBetaYPost[,1] = Group1VarY

  ## starting values

  fYPost[,1,,] = 0
  sigmaYPost[,1] = 1
  gammaYPost[,1,] = 0

  for (nc in 1 : nChains) {
    betaYPost[nc,1,] = c(mean(y), rep(0, dt), rep(0, nCatCols))
  }

  wYPost[,1] = 0.1

  if (length(whichCat) > 0) {
    designY = cbind(rep(1, n), tMat, xCat2)
  } else {
    designY = cbind(rep(1, n), tMat)
  }

  tauYPost[,1,] = 1



  ## We need to do SVD outside of MCMC first
  kernMat = array(NA, dim=c(pCont,n,n))
  vecMat = array(NA, dim=c(pCont,n,n))
  invMat = array(NA, dim=c(pCont,n,n))
  eigMat = array(NA, dim=c(pCont,n))
  nEigMat = rep(NA, pCont)
  logDet = rep(NA, pCont)


  for (j in 1 : pCont) {
    kernMat[j,,] = diag(n)
    for (n1 in 1 : n) {
      for (n2 in 1 : n1) {
        #kernMat[j,n1,n2] = exp(-((x[n1,j] - x[n2,j])^2)/band)
        kernMat[j,n1,n2] = exp(-abs(xCont[n1,j] - xCont[n2,j])/band)
        kernMat[j,n2,n1] = kernMat[j,n1,n2]
      }
    }

    eig = eigen(kernMat[j,,])
    eigMat[j,] = eig$values
    vecMat[j,,] = eig$vectors
    logDet[j] = determinant(kernMat[j,,])$modulus
    invMat[j,,] = solve(kernMat[j,,])

  }

  for (j in 1 : pCont) {
    fYPost[1,1,,j] = rep(0,n)
    fYPost[2,1,,j] = rep(0,n)
  }


  for (ni in 2 : nScans) {
    for (nc in 1 : nChains) {
      if (nc == 1 & ni %% 100 == 0) print(ni)

      ################# outcome model ##################

      ## update sigma
      if (type == "continuous") {
        SSM = 0
        for (j in 1 : pCont) {
          if (gammaYPost[nc,ni-1,j] == 1) {
            SSM = SSM + t(fYPost[nc,ni-1,,j]) %*% invMat[j,,] %*% fYPost[nc,ni-1,,j] /
              tauYPost[nc,ni-1,j]
          }
        }

        aStar = a + n*(sum(gammaYPost[nc,ni-1,])+1)/2
        bStar = b + (t(y - designY %*% betaYPost[nc,ni-1,] - apply(fYPost[nc,ni-1,,], 1, sum)) %*%
                       (y - designY %*% betaYPost[nc,ni-1,] - apply(fYPost[nc,ni-1,,], 1, sum)) / 2) +
          SSM/2

        if (length(whichCat) > 0) {
          aStar = a + n*(sum(gammaYPost[nc,ni-1,1:pCont])+1)/2 +
            gammaYPost[nc,ni-1,-c(1:pCont)]*sapply(colsCat, length)/2
          bStar = b + (t(y - designY %*% betaYPost[nc,ni-1,] - apply(fYPost[nc,ni-1,,], 1, sum)) %*%
                         (y - designY %*% betaYPost[nc,ni-1,] - apply(fYPost[nc,ni-1,,], 1, sum)) / 2) +
            (SSM/2) + sum(betaYPost[nc,ni-1,-c(1:(1+dt))]^2)/(2*Group1VarY)
        }

        sigmaYPost[nc,ni] = 1/rgamma(1, aStar, bStar)
        Zy = y
      } else {
        sigmaYPost[nc,ni] = 1
        Zy = rep(NA, n)

        meanZy = designY %*% betaYPost[nc,ni-1,] + apply(fYPost[nc,ni-1,,], 1, sum)

        Zy[y==1] = truncnorm::rtruncnorm(sum(y==1), a=0, mean = meanZy[y==1], sd=1)
        Zy[y==0] = truncnorm::rtruncnorm(sum(y==0), b=0, mean = meanZy[y==0], sd=1)
      }

      ## update tauYPost
      for (j in 1 : pCont) {
        if (gammaYPost[nc,ni-1,j] == 1) {
          tauYPost[nc,ni,j] = 1/rgamma(1, e + n/2, f + t(fYPost[nc,ni-1,,j]) %*%
                                         invMat[j,,] %*% fYPost[nc,ni-1,,j]/(2*sigmaYPost[nc,ni]))
        } else {
          tauYPost[nc,ni,j] = 1/rgamma(1,e,f)
        }
      }

      ## Update w
      wYPost[nc,ni] = rbeta(1, c + sum(gammaYPost[nc,ni-1,] == 1),
                            d + sum(gammaYPost[nc,ni-1,] == 0))

      ## update beta
      tempCols = 1:(dt + 1)
      if (nCatCols == 0) {
        ytilde = Zy - apply(fYPost[nc,ni-1,,], 1, sum)
      } else if (nCatCols == 1) {
        ytilde = Zy - apply(fYPost[nc,ni-1,,], 1, sum) -
          designY[,-tempCols] * betaYPost[nc,ni-1,-tempCols]
      } else {
        ytilde = Zy - apply(fYPost[nc,ni-1,,], 1, sum) -
          designY[,-tempCols] %*% betaYPost[nc,ni-1,-tempCols]
      }

      mu = solve(t(designY[,tempCols]) %*% designY[,tempCols]) %*%
        t(designY[,tempCols]) %*% ytilde
      sig = sigmaYPost[nc,ni]*solve(t(designY[,tempCols]) %*% designY[,tempCols])
      betaYPost[nc,ni,1:(dt+1)] = mvtnorm::rmvnorm(1, mu, sig)

      ## Update categorical covariate parameters if needed
      if (length(whichCat) > 0) {
        sigmaBetaYPost[nc,ni] = 1/rgamma(1, e + sum(gammaYPost[nc,ni-1,-c(1:pCont)]*
                                                      sapply(colsCat, length))/2,
                                         f + sum(betaYPost[nc,ni-1,-c(1:(1+dt))]^2)/
                                           (2*sigmaYPost[nc,ni]))
        Group1VarY = sigmaBetaYPost[nc,ni]



        ## Update regression coefficients and variable inclusion parameters
        tempBeta = betaYPost[nc,ni-1,]
        tempBeta[1:(dt+1)] = betaYPost[nc,ni,1:(dt+1)]

        for (j in 1 : pCat) {
          tempCols = colsCat[[j]]

          PriorSigmaY = Group1VarY*sigmaYPost[nc,ni]*diag(length(tempCols))
          muY = rep(0, length(tempCols))

          yStar = Zy - designY[,-tempCols] %*% tempBeta[-tempCols] -
            apply(fYPost[nc,ni-1,,], 1, sum)

          ## probability of being in group zero
          p0 = log(1 - wYPost[nc,ni])

          ## probability of being in top group
          muVar = solve(t(designY[,tempCols]) %*% designY[,tempCols]  / sigmaYPost[nc,ni] +
                          solve(PriorSigmaY))
          muBeta = muVar %*% (t(designY[,tempCols]) %*% yStar/sigmaYPost[nc,ni] +
                                solve(PriorSigmaY) %*% muY)
          p1 = log(wYPost[nc,ni]) + mvtnorm::dmvnorm(rep(0, length(tempCols)),
                                            mean=muY, sigma=PriorSigmaY, log=TRUE) -
            mvtnorm::dmvnorm(rep(0, length(tempCols)), mean=muBeta, sigma=muVar, log=TRUE)

          maxlog = max(p0,p1)

          p0new = exp(-maxlog + p0)
          p1new = exp(-maxlog + p1)

          gammaYPost[nc,ni,j+pCont] = sample(0:1, size=1, p=c(p0new,p1new))

          tempBeta[tempCols] = rep(0, length(tempCols))
          if (gammaYPost[nc,ni,j+pCont] == 1) tempBeta[tempCols] = mvtnorm::rmvnorm(1, muBeta, sigma=muVar)
        }
        betaYPost[nc,ni,] = tempBeta
      }

      ## update f(x_j)
      tempF = fYPost[nc,ni-1,,]
      for (j in 1 : pCont) {
        ytilde = Zy - designY %*% betaYPost[nc,ni,] - apply(tempF[,-j], 1, sum)

        tempObj1 = tauYPost[nc,ni,j]*eigMat[j,] /
          (1 + tauYPost[nc,ni,j]*eigMat[j,])

        VarMat = t(t(vecMat[j,,]) * tempObj1) %*%
          t(vecMat[j,,])

        VarMatInv = diag(n) + t((t(vecMat[j,,]) *
                                   (1/eigMat[j,]))) %*% t(vecMat[j,,])/tauYPost[nc,ni,j]

        mu = VarMat %*% ytilde

        logNumerator = -(n/2)*log(2*pi*sigmaYPost[nc,ni]*tauYPost[nc,ni,j]) -
          0.5*logDet[j] + log(wYPost[nc,ni])
        logDenominator = -(n/2)*log(2*pi*sigmaYPost[nc,ni]) -
          0.5*sum(log(tauYPost[nc,ni,j]*eigMat[j,] / (1 + tauYPost[nc,ni,j]*eigMat[j,]))) -
          0.5*t(rep(0,n) - mu) %*% VarMatInv %*% (rep(0,n) - mu)/sigmaYPost[nc,ni]

        logP1 = logNumerator - logDenominator
        logP0 = log(1 - wYPost[nc,ni])

        maxP = max(logP1, logP0)

        P1 = exp(logP1 - maxP)
        P0 = exp(logP0 - maxP)

        P = as.numeric(P1 / (P1 + P0))

        gammaYPost[nc,ni,j] = rbinom(1, 1, P)

        if (gammaYPost[nc,ni,j] == 0) {
          tempF[,j] = rep(0, n)
        } else {
          ## update rj first
          rj = rnorm(n, mean=0, sd = sqrt(sigmaYPost[nc,ni]*tempObj1))

          tempF[,j] = VarMat %*% ytilde + vecMat[j,,] %*% as.vector(rj)
        }
      }
      fYPost[nc,ni,,] = tempF
    }
  }

  keep = seq(nBurn + 1, nScans, by=thin)

  l = list(f = fYPost[,keep,,],
           beta = betaYPost[,keep,],
           gamma = gammaYPost[,keep,],
           sigma = sigmaYPost[,keep],
           sigmaBeta = tauYPost[,keep,])

  return(l)

}






GPTreatmentMCMC = function(t, x, band = 1, type="continuous",
                           whichCat,
                           nBurn=1000, nScans=5000,
                           thin=4, nChains=2, a=0.001, b=0.001,
                           c=2, d=dim(x)[2], e=0.5, f=0.5) {

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    pCat = 0
    nCatCols = 0
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]
    pCat = p - pCont

    lengthCat = c()
    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }

    for (j in 1 : pCat) {
      colsCat[[j]] = colsCat[[j]] + 1
    }

    nCatCols = dim(xCat2)[2]
  }

  ## parameters and prior specification

  n = length(t)

  ## add noise to observations with input values too similar
  for (j in 1 : pCont) {
    if (length(unique(xCont[,j])) < n) xCont[,j] = xCont[,j] + rnorm(n, sd=.00001)
  }

  fTPost = array(NA, dim=c(nChains, nScans, n, pCont))
  sigmaTPost = array(NA, dim=c(nChains, nScans))
  gammaTPost = array(NA, dim=c(nChains, nScans, p))
  betaTPost = array(NA, dim=c(nChains, nScans, 1+nCatCols))
  tauTPost = array(NA, dim=c(nChains, nScans, pCont))
  wTPost = array(NA, dim=c(nChains, nScans))
  sigmaBetaTPost = array(NA, dim=c(nChains, nScans))

  Group1VarT = 1
  sigmaBetaTPost[,1] = Group1VarT

  ## starting values

  fTPost[,1,,] = 0
  sigmaTPost[,1] = 1
  gammaTPost[,1,] = 0

  for (nc in 1 : nChains) {
    betaTPost[nc,1,] = c(mean(y), rep(0, nCatCols))
  }

  wTPost[,1] = 0.1

  if (length(whichCat) > 0) {
    designT = cbind(rep(1, n), xCat2)
  } else {
    designT = cbind(rep(1, n))
  }

  tauTPost[,1,] = 1



  ## We need to do SVD outside of MCMC first
  kernMat = array(NA, dim=c(pCont,n,n))
  vecMat = array(NA, dim=c(pCont,n,n))
  invMat = array(NA, dim=c(pCont,n,n))
  eigMat = array(NA, dim=c(pCont,n))
  nEigMat = rep(NA, pCont)
  logDet = rep(NA, pCont)


  for (j in 1 : pCont) {
    kernMat[j,,] = diag(n)
    for (n1 in 1 : n) {
      for (n2 in 1 : n1) {
        #kernMat[j,n1,n2] = exp(-((x[n1,j] - x[n2,j])^2)/band)
        kernMat[j,n1,n2] = exp(-abs(xCont[n1,j] - xCont[n2,j])/band)
        kernMat[j,n2,n1] = kernMat[j,n1,n2]
      }
    }

    eig = eigen(kernMat[j,,])
    eigMat[j,] = eig$values
    vecMat[j,,] = eig$vectors
    logDet[j] = determinant(kernMat[j,,])$modulus
    invMat[j,,] = solve(kernMat[j,,])

  }

  for (j in 1 : pCont) {
    fTPost[1,1,,j] = rep(0,n)
    fTPost[2,1,,j] = rep(0,n)
  }


  for (ni in 2 : nScans) {
    for (nc in 1 : nChains) {
      if (nc == 1 & ni %% 100 == 0) print(ni)

      ################# outcome model ##################

      ## update sigma
      if (type == "continuous") {
        SSM = 0
        for (j in 1 : pCont) {
          if (gammaTPost[nc,ni-1,j] == 1) {
            SSM = SSM + t(fTPost[nc,ni-1,,j]) %*% invMat[j,,] %*% fTPost[nc,ni-1,,j] /
              tauTPost[nc,ni-1,j]
          }
        }

        aStar = a + n*(sum(gammaTPost[nc,ni-1,])+1)/2
        bStar = b + (t(t- designT %*% betaTPost[nc,ni-1,] - apply(fTPost[nc,ni-1,,], 1, sum)) %*%
                       (t- designT %*% betaTPost[nc,ni-1,] - apply(fTPost[nc,ni-1,,], 1, sum)) / 2) +
          SSM/2

        if (length(whichCat) > 0) {
          aStar = a + n*(sum(gammaTPost[nc,ni-1,1:pCont])+1)/2 +
            gammaTPost[nc,ni-1,-c(1:pCont)]*sapply(colsCat, length)/2
          bStar = b + (t(t- designT %*% betaTPost[nc,ni-1,] - apply(fTPost[nc,ni-1,,], 1, sum)) %*%
                         (t- designT %*% betaTPost[nc,ni-1,] - apply(fTPost[nc,ni-1,,], 1, sum)) / 2) +
            (SSM/2) + sum(betaTPost[nc,ni-1,-c(1)]^2)/(2*Group1VarT)
        }

        sigmaTPost[nc,ni] = 1/rgamma(1, aStar, bStar)
        Zt = t
      } else {
        sigmaTPost[nc,ni] = 1
        Zt = rep(NA, n)

        meanZt = designT %*% betaTPost[nc,ni-1,] + apply(fTPost[nc,ni-1,,], 1, sum)

        Zt[t==1] = truncnorm::rtruncnorm(sum(t==1), a=0, mean = meanZt[t==1], sd=1)
        Zt[t==0] = truncnorm::rtruncnorm(sum(t==0), b=0, mean = meanZt[t==0], sd=1)
      }

      ## update tauTPost
      for (j in 1 : pCont) {
        if (gammaTPost[nc,ni-1,j] == 1) {
          tauTPost[nc,ni,j] = 1/rgamma(1, e + n/2, f + t(fTPost[nc,ni-1,,j]) %*%
                                         invMat[j,,] %*% fTPost[nc,ni-1,,j]/(2*sigmaTPost[nc,ni]))
        } else {
          tauTPost[nc,ni,j] = 1/rgamma(1,e,f)
        }
      }

      ## Update w
      wTPost[nc,ni] = rbeta(1, c + sum(gammaTPost[nc,ni-1,] == 1),
                            d + sum(gammaTPost[nc,ni-1,] == 0))

      ## update beta
      tempCols = 1
      if (nCatCols == 0) {
        ttilde = Zt - apply(fTPost[nc,ni-1,,], 1, sum)
      } else if (nCatCols == 1) {
        ttilde = Zt - apply(fTPost[nc,ni-1,,], 1, sum) -
          designT[,-tempCols] * betaTPost[nc,ni-1,-tempCols]
      } else {
        ttilde = Zt - apply(fTPost[nc,ni-1,,], 1, sum) -
          designT[,-tempCols] %*% betaTPost[nc,ni-1,-tempCols]
      }

      mu = solve(t(designT[,tempCols]) %*% designT[,tempCols]) %*%
        t(designT[,tempCols]) %*% ttilde
      sig = sigmaTPost[nc,ni]*solve(t(designT[,tempCols]) %*% designT[,tempCols])
      betaTPost[nc,ni,1] = mvtnorm::rmvnorm(1, mu, sig)

      ## Update categorical covariate parameters if needed
      if (length(whichCat) > 0) {
        sigmaBetaTPost[nc,ni] = 1/rgamma(1, e + sum(gammaTPost[nc,ni-1,-c(1:pCont)]*
                                                      sapply(colsCat, length))/2,
                                         f + sum(betaTPost[nc,ni-1,-c(1)]^2)/
                                           (2*sigmaTPost[nc,ni]))
        Group1VarT = sigmaBetaTPost[nc,ni]



        ## Update regression coefficients and variable inclusion parameters
        tempBeta = betaTPost[nc,ni-1,]
        tempBeta[1] = betaTPost[nc,ni,1]

        for (j in 1 : pCat) {
          tempCols = colsCat[[j]]

          PriorSigmaT = Group1VarT*sigmaTPost[nc,ni]*diag(length(tempCols))
          muT = rep(0, length(tempCols))

          if (length(tempBeta) - length(tempCols) == 1) {
            tStar = Zt - designT[,-tempCols] * tempBeta[-tempCols] -
              apply(fTPost[nc,ni-1,,], 1, sum)
          } else {
            tStar = Zt - designT[,-tempCols] %*% tempBeta[-tempCols] -
              apply(fTPost[nc,ni-1,,], 1, sum)
          }

          ## probability of being in group zero
          p0 = log(1 - wTPost[nc,ni])

          ## probability of being in top group
          muVar = solve(t(designT[,tempCols]) %*% designT[,tempCols]  / sigmaTPost[nc,ni] +
                          solve(PriorSigmaT))
          muBeta = muVar %*% (t(designT[,tempCols]) %*% tStar/sigmaTPost[nc,ni] +
                                solve(PriorSigmaT) %*% muT)
          p1 = log(wTPost[nc,ni]) + mvtnorm::dmvnorm(rep(0, length(tempCols)),
                                            mean=muT, sigma=PriorSigmaT, log=TRUE) -
            mvtnorm::dmvnorm(rep(0, length(tempCols)), mean=muBeta, sigma=muVar, log=TRUE)

          maxlog = max(p0,p1)

          p0new = exp(-maxlog + p0)
          p1new = exp(-maxlog + p1)

          gammaTPost[nc,ni,j+pCont] = sample(0:1, size=1, p=c(p0new,p1new))

          tempBeta[tempCols] = rep(0, length(tempCols))
          if (gammaTPost[nc,ni,j+pCont] == 1) tempBeta[tempCols] = mvtnorm::rmvnorm(1, muBeta, sigma=muVar)
        }
        betaTPost[nc,ni,] = tempBeta
      }

      ## update f(x_j)
      tempF = fTPost[nc,ni-1,,]
      for (j in 1 : pCont) {
        ttilde = Zt - designT %*% betaTPost[nc,ni,] - apply(tempF[,-j], 1, sum)

        tempObj1 = tauTPost[nc,ni,j]*eigMat[j,] /
          (1 + tauTPost[nc,ni,j]*eigMat[j,])

        VarMat = t(t(vecMat[j,,]) * tempObj1) %*%
          t(vecMat[j,,])

        VarMatInv = diag(n) + t((t(vecMat[j,,]) *
                                   (1/eigMat[j,]))) %*% t(vecMat[j,,])/tauTPost[nc,ni,j]

        mu = VarMat %*% ttilde

        logNumerator = -(n/2)*log(2*pi*sigmaTPost[nc,ni]*tauTPost[nc,ni,j]) -
          0.5*logDet[j] + log(wTPost[nc,ni])
        logDenominator = -(n/2)*log(2*pi*sigmaTPost[nc,ni]) -
          0.5*sum(log(tauTPost[nc,ni,j]*eigMat[j,] / (1 + tauTPost[nc,ni,j]*eigMat[j,]))) -
          0.5*t(rep(0,n) - mu) %*% VarMatInv %*% (rep(0,n) - mu)/sigmaTPost[nc,ni]

        logP1 = logNumerator - logDenominator
        logP0 = log(1 - wTPost[nc,ni])

        maxP = max(logP1, logP0)

        P1 = exp(logP1 - maxP)
        P0 = exp(logP0 - maxP)

        P = as.numeric(P1 / (P1 + P0))

        gammaTPost[nc,ni,j] = rbinom(1, 1, P)

        if (gammaTPost[nc,ni,j] == 0) {
          tempF[,j] = rep(0, n)
        } else {
          ## update rj first
          rj = rnorm(n, mean=0, sd = sqrt(sigmaTPost[nc,ni]*tempObj1))

          tempF[,j] = VarMat %*% ttilde + vecMat[j,,] %*% as.vector(rj)
        }
      }
      fTPost[nc,ni,,] = tempF
    }
  }

  keep = seq(nBurn + 1, nScans, by=thin)

  l = list(f = fTPost[,keep,,],
           beta = betaTPost[,keep,],
           gamma = gammaTPost[,keep,],
           sigma = sigmaTPost[,keep],
           sigmaBeta = tauTPost[,keep,])

  return(l)

}




DRmcmcCut = function(y, t, x, whichCat,
                     nChains, totalScans, PostT, PostY,
                     modY = "GP", modT = "GP",
                     dfY = NULL, dfT = NULL, nBoot = 100,
                     lower=0.05, upper=0.95) {

  n = length(t)

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    pCat = 0
    nCatCols = 0
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]
    pCat = p - pCont

    lengthCat = c()
    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }
  }

  if (modT == "GP") {
    betaTPost = PostT$beta
    fTPost = PostT$f
  } else {
    betaTPost = PostT$beta
    designT = cbind(rep(1,n))

    for (j in 1 : pCont) {
      tempT = scale(splines::ns(xCont[,j], dfT))
      designT = cbind(designT, tempT)
    }

    if (length(whichCat) > 0) {
      designT = cbind(designT, xCat2)
    }
  }


  if (modY == "GP") {
    betaYPost = PostY$beta
    fYPost = PostY$f
  } else {
    betaYPost = PostY$beta
    designY = cbind(rep(1,n), t)

    for (j in 1 : pCont) {
      tempY = scale(splines::ns(xCont[,j], dfY))
      designY = cbind(designY, tempY)
    }
    if (length(whichCat) > 0) {
      designY = cbind(designY, xCat2)
    }
  }


  DRestPost = matrix(NA, nChains, totalScans)

  for (nc in 1 : nChains) {
    for (ni in 1 : totalScans) {
      if (modT == "GP") {
        if (length(whichCat) > 0) {
          designT = cbind(rep(1,n), xCat2)
          eHat = pnorm(designT %*% betaTPost[nc,ni,] + apply(fTPost[nc,ni,,], 1, sum))
        } else {
          eHat = pnorm(betaTPost[nc,ni] + apply(fTPost[nc,ni,,], 1, sum))
        }
      } else {
        eHat = pnorm(designT %*% betaTPost[nc,ni,])
      }

      if (modY == "GP") {
        if (length(whichCat) > 0) {
          muHat1 = cbind(rep(1,n), rep(1,n), xCat2) %*% betaYPost[nc,ni,] +
            apply(fYPost[nc,ni,,], 1, sum)
          muHat0 = cbind(rep(1,n), rep(0,n), xCat2) %*% betaYPost[nc,ni,] +
            apply(fYPost[nc,ni,,], 1, sum)
        } else {
          muHat1 = cbind(rep(1,n), rep(1,n)) %*% betaYPost[nc,ni,] +
            apply(fYPost[nc,ni,,], 1, sum)
          muHat0 = cbind(rep(1,n), rep(0,n)) %*% betaYPost[nc,ni,] +
            apply(fYPost[nc,ni,,], 1, sum)
        }
      } else {
        muHat1 = cbind(rep(1,n), rep(1,n), designY[,-c(1,2)]) %*% betaYPost[nc,ni,]
        muHat0 = cbind(rep(1,n), rep(0,n), designY[,-c(1,2)]) %*% betaYPost[nc,ni,]
      }

      eHat[eHat < lower] = lower
      eHat[eHat > upper] = upper

      part1 = mean(t*y/eHat - (t - eHat)*muHat1/eHat)
      part2 = mean((1 - t)*y/(1 - eHat) + (t - eHat)*muHat0/(1 - eHat))

      DRestPost[nc,ni] = part1 - part2

    }
  }


  BootMat = array(NA, dim=c(nBoot, nChains, totalScans))

  for (ii in 1 : nBoot) {
    samp = sample(1:n, n, replace=TRUE)
    for (nc in 1 : nChains) {
      for (ni in 1 : totalScans) {
        if (modT == "GP") {
          if (length(whichCat) > 0) {
            designT = cbind(rep(1,n), xCat2)
            eHat = pnorm(designT[samp,] %*% betaTPost[nc,ni,] +
                           apply(fTPost[nc,ni,,], 1, sum)[samp])
          } else {
            eHat = pnorm(betaTPost[nc,ni] + apply(fTPost[nc,ni,,], 1, sum)[samp])
          }
        } else {
          eHat = pnorm(designT[samp,] %*% betaTPost[nc,ni,])
        }

        if (modY == "GP") {
          if (length(whichCat) > 0) {
            muHat1 = cbind(rep(1,n), rep(1,n), xCat2[samp,]) %*% betaYPost[nc,ni,] +
              apply(fYPost[nc,ni,,], 1, sum)[samp]
            muHat0 = cbind(rep(1,n), rep(0,n), xCat2[samp,]) %*% betaYPost[nc,ni,] +
              apply(fYPost[nc,ni,,], 1, sum)[samp]
          } else {
            muHat1 = cbind(rep(1,n), rep(1,n)) %*% betaYPost[nc,ni,] +
              apply(fYPost[nc,ni,,], 1, sum)[samp]
            muHat0 = cbind(rep(1,n), rep(0,n)) %*% betaYPost[nc,ni,] +
              apply(fYPost[nc,ni,,], 1, sum)[samp]
          }
        } else {
          muHat1 = cbind(rep(1,n), rep(1,n), designY[samp,-c(1,2)]) %*% betaYPost[nc,ni,]
          muHat0 = cbind(rep(1,n), rep(0,n), designY[samp,-c(1,2)]) %*% betaYPost[nc,ni,]
        }

        eHat[eHat < lower] = lower
        eHat[eHat > upper] = upper

        part1 = mean(t[samp]*y[samp]/eHat - (t[samp] - eHat)*muHat1/eHat)
        part2 = mean((1 - t[samp])*y[samp]/(1 - eHat) + (t[samp] - eHat)*muHat0/(1 - eHat))

        BootMat[ii,nc,ni] = part1 - part2

      }
    }
  }

  totalVar = mean(apply(BootMat, 1, sd)^2) + var(apply(BootMat, 1, mean))
  totalSD = sqrt(totalVar)

  l = list(est = mean(DRestPost),
           se = totalSD,
           BootQuantile = quantile(BootMat, c(.025, .975)))
}



DRmcmcContinuousCut = function(y, t, tMat, tMatNew, x, whichCat,
                               locations=seq(quantile(t, .05), quantile(t, .95), length=20),
                               nChains, totalScans, PostT, PostY,
                               modY = "GP", modT = "GP",
                               dfY = NULL, dfT = NULL, nBoot = 100, threshold=0.1) {


  n = length(t)

  if (is.vector(tMat) == TRUE) {
    dt = 1
  } else {
    dt = dim(tMat)[2]
  }

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    pCat = 0
    nCatCols = 0
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]
    pCat = p - pCont

    lengthCat = c()
    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }
  }

  if (modT == "GP") {
    betaTPost = PostT$beta
    fTPost = PostT$f
  } else {
    betaTPost = PostT$beta
    designT = cbind(rep(1,n))

    for (j in 1 : pCont) {
      tempT = scale(splines::ns(xCont[,j], dfT))
      designT = cbind(designT, tempT)
    }

    if (length(whichCat) > 0) {
      designT = cbind(designT, xCat2)
    }
  }


  if (modY == "GP") {
    betaYPost = PostY$beta
    fYPost = PostY$f
  } else {
    betaYPost = PostY$beta
    designY = cbind(rep(1,n), tMat)

    for (j in 1 : pCont) {
      tempY = scale(splines::ns(xCont[,j], dfY))
      designY = cbind(designY, tempY)
    }

    if (length(whichCat) > 0) {
      designY = cbind(designY, xCat2)
    }
  }

  DRestPost = array(NA, dim=c(nChains, totalScans, dim(tMatNew)[1]))

  for (nc in 1 : nChains) {
    for (ni in 1 : totalScans) {
      if (ni %% 100 == 0 & nc==1) print(ni)

      if (modT == "GP") {
        applyT = apply(fTPost[nc,ni,,], 1, sum)
        designT = cbind(rep(1, n))
        if (length(whichCat) > 0) {
          designT = cbind(designT, xCat2)
          pHat = designT %*% betaTPost[nc,ni,] + applyT
        } else {
          pHat = betaTPost[nc,ni] + applyT
        }
      } else {
        pHat = designT %*% betaTPost[nc,ni,]
      }

      if (modY == "GP") {
        applyY = apply(fYPost[nc,ni,,], 1, sum)
        designY = cbind(rep(1,n), tMat)
        if (length(whichCat) > 0) {
          designY = cbind(designY, xCat2)
        }
        muHat = designY %*% betaYPost[nc,ni,] +
          applyY
      } else {
        muHat = cbind(rep(1,n), tMat, designY[,-c(1:(1+dt))]) %*% betaYPost[nc,ni,]
        xBetaY = designY[,-c(1:(1+dt))] %*% betaYPost[nc,ni,-c(1:(1+dt))]
      }

      pSig = sqrt(PostT$sigma[nc,ni])
      pDens = dnorm(t, mean=pHat, sd=pSig)
      pHatA = rep(NA, n)
      muHatA = rep(NA, n)

      if (modY == "GP") {
        for (i in 1 : n) {
          pHatA[i] = mean(dnorm(t[i], mean=pHat, sd=pSig))
          if (length(whichCat) > 0) {
            muHatA[i] = mean(cbind(rep(1,n),
                                   t(matrix(rep(tMat[i,], n), ncol=n)), xCat2) %*%
                               betaYPost[nc,ni,] + applyY)
          } else {
            muHatA[i] = mean(cbind(rep(1,n),
                                   t(matrix(rep(tMat[i,], n), ncol=n))) %*%
                               betaYPost[nc,ni,] + applyY)
          }
        }
      } else {
        for (i in 1 : n) {
          pHatA[i] = mean(dnorm(t[i], mean=pHat, sd=pSig))
          muHatA[i] = mean(cbind(rep(1, n),
                                 t(matrix(rep(tMat[i,], n), ncol=n))) %*%
                             betaYPost[nc,ni,1:(1+dt)] + xBetaY)
        }
      }

      ratio = pHatA / pDens
      ratio[ratio < threshold] = threshold
      ratio[ratio > (1/threshold)] = (1/threshold)

      newOutcome = ((y - muHat) * ratio) + muHatA

      coefs = lm(newOutcome ~ tMat)$coefficients

      DRestPost[nc,ni,] = cbind(rep(1, dim(tMatNew)[1]), tMatNew) %*% coefs
    }
  }

  BootMat = array(NA, dim=c(nBoot, nChains, totalScans, dim(tMatNew)[1]))

  sampList = list()
  for (ii in 1 : nBoot) sampList[[ii]] = sample(1:n, n, replace=TRUE)

  for (nc in 1 : nChains) {
    for (ni in 1 : totalScans) {

      if (modT == "GP") {
        applyT = apply(fTPost[nc,ni,,], 1, sum)
      }

      if (modY == "GP") {
        applyY = apply(fYPost[nc,ni,,], 1, sum)
      } else {
        xBetaY = designY[,-c(1:(1+dt))] %*% betaYPost[nc,ni,-c(1:(1+dt))]
      }

      for (ii in 1 : nBoot) {
        if (ni %% 10 == 0 & nc==1 & ii==1) print(ni)

        if (modT == "GP") {
          if (length(whichCat) > 0) {
            pHat = designT[sampList[[ii]],] %*% betaTPost[nc,ni,] +
              applyT[sampList[[ii]]]
          } else {
            pHat = betaTPost[nc,ni] + applyT[sampList[[ii]]]
          }
        } else {
          pHat = designT[sampList[[ii]],] %*% betaTPost[nc,ni,]
        }

        if (modY == "GP") {
          muHat = designY[sampList[[ii]],] %*% betaYPost[nc,ni,] +
            applyY[sampList[[ii]]]
        } else {
          muHat = cbind(rep(1,n), tMat[sampList[[ii]],],
                        designY[sampList[[ii]],-c(1:(1+dt))]) %*% betaYPost[nc,ni,]
        }

        pSig = sqrt(PostT$sigma[nc,ni])
        pDens = dnorm(t[sampList[[ii]]], mean=pHat, sd=pSig)
        pHatA = rep(NA, n)
        muHatA = rep(NA, n)

        if (modY == "GP") {
          for (i in 1 : n) {
            pHatA[i] = mean(dnorm(t[sampList[[ii]][i]], mean=pHat, sd=pSig))
            if (length(whichCat) > 0) {
              muHatA[i] = mean(cbind(rep(1,n),
                                     t(matrix(rep(tMat[sampList[[ii]][i],], n), ncol=n)),
                                     xCat2[sampList[[ii]],]) %*%
                                 betaYPost[nc,ni,] + applyY[sampList[[ii]]])
            } else {
              muHatA[i] = mean(cbind(rep(1,n),
                                     t(matrix(rep(tMat[sampList[[ii]][i],], n), ncol=n))) %*%
                                 betaYPost[nc,ni,] + applyY[sampList[[ii]]])
            }
          }
        } else {
          for (i in 1 : n) {
            pHatA[i] = mean(dnorm(t[sampList[[ii]][i]], mean=pHat, sd=pSig))
            muHatA[i] = mean(cbind(rep(1, n),
                                   t(matrix(rep(tMat[sampList[[ii]][i],], n), ncol=n))) %*%
                               betaYPost[nc,ni,1:(1+dt)] + xBetaY[sampList[[ii]]])
          }
        }

        ratio = pHatA / pDens
        ratio[ratio < threshold] = threshold
        ratio[ratio > (1/threshold)] = (1/threshold)

        newOutcome = ((y[sampList[[ii]]] - muHat) * ratio) + muHatA

        coefs = lm(newOutcome ~ tMat[sampList[[ii]],])$coefficients

        BootMat[ii,nc,ni,] =  cbind(rep(1, dim(tMatNew)[1]), tMatNew) %*% coefs

      }
    }
  }

  totalVar = CIlower = CIupper = rep(NA, dim(tMatNew)[1])
  for (k in 1 : dim(tMatNew)[1]) {
    totalVar[k] = mean(apply(BootMat[,,,k], 1, sd)^2) + var(apply(BootMat[,,,k], 1, mean))
    CIlower[k] = quantile(BootMat[,,,k], .025)
    CIupper[k] = quantile(BootMat[,,,k], .975)
  }
  totalSD = sqrt(totalVar)

  l = list(est = apply(DRestPost, 3, mean),
           se = totalSD,
           CIlower = CIlower,
           CIupper = CIupper)

  return(l)
}




RegContinuousEst = function(x, locations=seq(quantile(t, .05), quantile(t, .95), length=20),
                            whichCat, nChains, totalScans, PostY,
                            modY = "GP", nBoot = 50,
                            dfY = NULL) {

  df = dfY

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    pCat = 0
    nCatCols = 0
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]
    pCat = p - pCont

    lengthCat = c()
    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }
  }


  if (modY == "GP") {
    betaYPost = PostY$beta
    fYPost = PostY$f
  } else {
    betaYPost = PostY$beta
    designY = scale(splines::ns(xCont[,1], dfY))

    for (j in 2 : pCont) {
      tempY = scale(splines::ns(xCont[,j], dfY))
      designY = cbind(designY, tempY)
    }

    if (length(whichCat) > 0) {
      designY = cbind(designY, xCat2)
    }
  }

  nLocations = length(locations)

  POUT = array(NA, dim=c(nChains, totalScans, nLocations))

  for (nc in 1 : nChains) {
    for (ni in 1 : totalScans) {
      for (nl in 1 : nLocations) {
        tempT = locations[nl]
        tempTmat = t(matrix(rep(c(tempT, tempT^2, tempT^3), n), ncol=n))

        if (modY == "GP") {
          if (length(whichCat) == 0) {
            POUT[nc,ni,nl] = mean(cbind(rep(1, n), tempTmat) %*% betaYPost[nc,ni,] +
                                    apply(fYPost[nc,ni,,], 1, sum))
          } else {
            POUT[nc,ni,nl] = mean(cbind(rep(1, n), tempTmat, xCat2) %*% betaYPost[nc,ni,] +
                                    apply(fYPost[nc,ni,,], 1, sum))
          }
        } else {
          totalMat = cbind(rep(1, n), tempTmat, designY)
          POUT[nc,ni,nl] = mean(totalMat %*% betaYPost[nc,ni,])
        }
      }
    }
  }


  POUTboot = array(NA, dim=c(nBoot, nChains, totalScans, nLocations))

  for (nb in 1 : nBoot) {
    samp = sample(1:n, n, replace=TRUE)
    for (nc in 1 : nChains) {
      for (ni in 1 : totalScans) {
        for (nl in 1 : nLocations) {
          tempT = locations[nl]
          tempTmat = t(matrix(rep(c(tempT, tempT^2, tempT^3), n), ncol=n))

          if (modY == "GP") {
            if (length(whichCat) == 0) {
              POUTboot[nb,nc,ni,nl] = mean(cbind(rep(1, n), tempTmat) %*% betaYPost[nc,ni,] +
                                             apply(fYPost[nc,ni,,], 1, sum)[samp])
            } else {
              POUTboot[nb,nc,ni,nl] = mean(cbind(rep(1, n), tempTmat, xCat2[samp,]) %*% betaYPost[nc,ni,] +
                                             apply(fYPost[nc,ni,,], 1, sum)[samp])
            }
          } else {
            totalMat = cbind(rep(1, n), tempTmat, designY[samp,])
            POUTboot[nb,nc,ni,nl] = mean(totalMat %*% betaYPost[nc,ni,])
          }
        }
      }
    }
  }

  est = apply(POUT, 3, mean)
  se = apply(POUTboot, 4, sd)
  CIlower = apply(POUTboot, 4, quantile, .025)
  CIupper = apply(POUTboot, 4, quantile, .975)

  return(list(est=est, se=se, CIlower=CIlower, CIupper=CIupper))
}







WAICoutcome = function(y, x, tMat, Post, whichCat,
                       type="continuous", modY, dfY = NULL,
                       nChains, totalScans) {

  df = dfY

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    pCat = 0
    nCatCols = 0
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]
    pCat = p - pCont

    lengthCat = c()
    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }
  }

  if (modY == "Splines") {
    designY = cbind(rep(1, n), tMat)

    for (j in 1 : pCont) {
      tempY = scale(splines::ns(xCont[,j], dfY))
      designY = cbind(designY, tempY)
    }

    if (length(whichCat) > 0) {
      designY = cbind(designY, xCat2)
    }
  }

  n = length(y)

  LHood = array(NA, dim=c(nChains, totalScans, n))

  if (type == "continuous") {
    if (modY == "Splines") {
      for (nc in 1 : nChains) {
        for (ni in 1 : totalScans) {
          linPred = designY %*% Post$beta[nc,ni,]
          LHood[nc,ni,] = (1/(sqrt(2*pi*Post$sigma[nc,ni]))) *
            exp((-(y - linPred)^2)/(2*Post$sigma[nc,ni]))
        }
      }
    } else {
      designY = cbind(rep(1, n), tMat)
      if (length(whichCat) > 0) {
        designY = cbind(designY, xCat2)
      }
      for (nc in 1 : nChains) {
        for (ni in 1 : totalScans) {
          linPred = designY %*% Post$beta[nc,ni,] + apply(Post$f[nc,ni,,], 1, sum)
          LHood[nc,ni,] = (1/(sqrt(2*pi*Post$sigma[nc,ni]))) *
            exp((-(y - linPred)^2)/(2*Post$sigma[nc,ni]))
        }
      }
    }
  } else {
    if (modY == "Splines") {
      for (nc in 1 : nChains) {
        for (ni in 1 : totalScans) {
          linPred = pnorm(designY %*% Post$beta[nc,ni,])
          LHood[nc,ni,] = (linPred^y) * (1 - linPred)^(1-y)
        }
      }
    } else {
      designY = cbind(rep(1, n), tMat)
      if (length(whichCat) > 0) {
        designY = cbind(designY, xCat2)
      }
      for (nc in 1 : nChains) {
        for (ni in 1 : totalScans) {
          linPred = pnorm(designY %*% Post$beta[nc,ni,] + apply(Post$f[nc,ni,,], 1, sum))
          LHood[nc,ni,] = (linPred^y) * (1 - linPred)^(1-y)
        }
      }
    }
  }

  return(-2*(sum(log(apply(LHood, 3, mean))) - sum(apply(log(LHood), 3, sd)^2)))
}

WAICtreatment = function(t, x, Post, modT, whichCat,
                         type="binary", dfT = NULL,
                         nChains, totalScans) {

  df = dfT

  if (length(whichCat) == 0) {
    xCont = scale(x)
    p = dim(x)[2]
    pCont = p
    pCat = 0
    nCatCols = 0
  } else {
    p = dim(x)[2]
    xCat = x[,whichCat]
    xCont = scale(x[,-whichCat])
    pCont = dim(xCont)[2]
    pCat = p - pCont

    lengthCat = c()
    if (length(whichCat) == 1) {
      lengthCat[1] = length(unique(xCat))
    } else {
      for (j2 in 1 : length(whichCat)) {
        lengthCat[j2] = length(unique(xCat[,j2]))
      }
    }
    xCat2 = matrix(NA, n, sum(lengthCat) - length(whichCat))

    colsCat = list()
    colsCat[[1]] = 1:(lengthCat[1]-1)

    if (length(whichCat) == 1) {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat == unique(xCat)[j3])
      }
    } else {
      for (j3 in 1 : length(colsCat[[1]])) {
        xCat2[,colsCat[[1]][j3]] = as.numeric(xCat[,1] == unique(xCat[,1])[j3])
      }
    }

    if (length(whichCat) > 1) {
      for (j2 in 2 : length(whichCat)) {
        colsCat[[j2]] = (cumsum(lengthCat-1)[j2-1] + 1) : (cumsum(lengthCat-1)[j2])
        for (j3 in 1 : length(colsCat[[j2]])) {
          xCat2[,colsCat[[j2]][j3]] = as.numeric(xCat[,j2] == unique(xCat[,j2])[j3])
        }
      }
    }
  }

  if (modT == "Splines") {
    designT = cbind(rep(1, n))

    for (j in 1 : pCont) {
      tempT = scale(splines::ns(xCont[,j], dfT))
      designT = cbind(designT, tempT)
    }

    if (length(whichCat) > 0) {
      designT = cbind(designT, xCat2)
    }
  }

  n = length(t)

  LHood = array(NA, dim=c(nChains, totalScans, n))

  if (type == "continuous") {
    if (modT == "Splines") {
      for (nc in 1 : nChains) {
        for (ni in 1 : totalScans) {
          linPred = designT %*% Post$beta[nc,ni,]
          LHood[nc,ni,] = (1/(sqrt(2*pi*Post$sigma[nc,ni]))) *
            exp((-(t - linPred)^2)/(2*Post$sigma[nc,ni]))
        }
      }
    } else {
      designT = cbind(rep(1, n))
      if (length(whichCat) > 0) {
        designT = cbind(designT, xCat2)
      }
      for (nc in 1 : nChains) {
        for (ni in 1 : totalScans) {
          if (length(whichCat) > 0) {
            linPred = designT %*% Post$beta[nc,ni,] + apply(Post$f[nc,ni,,], 1, sum)
          } else {
            linPred = Post$beta[nc,ni] + apply(Post$f[nc,ni,,], 1, sum)
          }
          LHood[nc,ni,] = (1/(sqrt(2*pi*Post$sigma[nc,ni]))) *
            exp((-(t - linPred)^2)/(2*Post$sigma[nc,ni]))
        }
      }
    }
  } else {

    if (modT == "Splines") {
      for (nc in 1 : nChains) {
        for (ni in 1 : totalScans) {
          linPred = pnorm(designT %*% Post$beta[nc,ni,])
          LHood[nc,ni,] = (linPred^t) * (1 - linPred)^(1-t)
        }
      }
    } else {
      designT = cbind(rep(1, n))
      if (length(whichCat) > 0) {
        designT = cbind(designT, xCat2)
      }
      for (nc in 1 : nChains) {
        for (ni in 1 : totalScans) {
          if (length(whichCat) > 0) {
            linPred = pnorm(designT %*% Post$beta[nc,ni,] + apply(Post$f[nc,ni,,], 1, sum))
          } else {
            linPred = pnorm(Post$beta[nc,ni] + apply(Post$f[nc,ni,,], 1, sum))
          }
          LHood[nc,ni,] = (linPred^t) * (1 - linPred)^(1-t)
        }
      }
    }
  }

  return(-2*(sum(log(apply(LHood, 3, mean))) - sum(apply(log(LHood), 3, sd)^2)))
}

