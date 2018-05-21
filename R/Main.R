#' Estimate causal effect using doubly robust estimators in high dimensions
#'
#' This function will take in the observed data and estimate a treatment effect.
#' y,x, and t must all be supplied, though all other parameters have pre-set values
#' the user can proceed with unless they wish to change the prior specification.
#'
#' @param y              The outcome to be analyzed
#' @param t              The treatment whose causal effect is to be estimated
#' @param x              An n by p matrix of potential confounders
#' @param whichCat       A vector of indices that indicate which variables in x are categorical.
#'                       The default is c(), which means the program defaults to assuming all
#'                       covariates are continuous
#' @param y_type         A categorical variable indicating whether y is binary or continuous.
#'                       Possible values are "binary" or "continuous", and the program defaults
#'                       to "continuous"
#' @param dfY            the degrees of freedom of the splines used to model the relationship between
#'                       the covariates and the outcome. If the user wants to use gaussian process
#'                       priors instead of splines then dfY should be set to "GP"
#' @param dfT            the degrees of freedom of the splines used to model the relationship between
#'                       the covariates and the treatment. If the user wants to use gaussian process
#'                       priors instead of splines then dfY should be set to "GP"
#' @param nScans         The number of MCMC scans to run
#' @param nBurn           The number of MCMC scans that will be dropped as a burn-in
#' @param thin           This number represents how many iterations between each scan
#'                       that is kept
#' @param thetaA         The first parameter of the beta prior on the overall sparsity level
#' @param thetaB         The second parameter of the beta prior on the overall sparsity level
#' @param band           The bandwidth parameter for the gaussian process kernel function
#' @param nBoot          The number of resampling iterations to use when estimating the credible
#'                       intervals
#' @param lower          The lowest value the estimated propensity score can take in the DR estimator. 
#'                       This parameter defaults to 0 so that the propensity score is not trimmed.
#' @param upper          The largest value the propensity score can take in the DR estimator. This
#'                       parameter defaults to 1 so that the propensity score is not trimmed.
#'
#'
#' @return A list of values that contain the treatment effect, confidence interval for the
#'         treatment effect, WAIC for the chosen treatment model and outcome model,
#'
#' @export
#' @examples
#'
#' ## p can be larger than n, but we keep the number of covariates small here
#' ## just for illustration so that the code will finish faster
#' n = 200
#' p = 20
#' x = matrix(rnorm(n*p), n, p)
#' t = rbinom(n, 1, p=pnorm(0.7*x[,1] + 0.3*x[,2]))
#' y = rnorm(n, mean=t + 0.3*x[,1] + 0.6*x[,2] + 0.5*x[,3], sd=1)
#'
#' est = DRbayes(y=y, t=t, x=x, nScans=2000, nBurn=1000, thin=2)

DRbayes = function(nScans = 20000, nBurn = 10000, thin = 10,
                   y, x, t, whichCat = c(), y_type = "continuous",
                   dfY = 1, dfT = 1, band = 3,
                   thetaA = 1, thetaB = 0.2*dim(x)[2],
                   nBoot=500, lower=0, upper=1) {

  n = dim(x)[1]
  p = dim(x)[2]

  x = scale(x)
  
  totalScans = floor((nScans - nBurn)/thin)

  if (length(unique(t)) > 2) {
    stop("t must be binary for this function")
  } else {

    if (dfY == "GP") {
      modY = "GP"
      dfY = NULL
      PostY = GPOutcomeMCMC(y=y, tMat=t, x=x, whichCat = whichCat, type=y_type,
                            nScans=nScans, nBurn=nBurn, thin=thin, band=band)
      WAICY = WAICoutcome(y=y, x=x, tMat=t, Post=PostY, modY="GP", type=y_type,
                          whichCat=whichCat, dfY = NULL, nChains = 2, 
                          totalScans = totalScans)
    } else {
      modY = "Splines"
      PostY = SplineOutcomeMCMC(y=y, tMat=t, x=x, whichCat = whichCat, type=y_type,
                                df=dfY, nScans=nScans, nBurn=nBurn, thin=thin)
      WAICY = WAICoutcome(y=y, x=x, tMat=t, Post=PostY, modY="Splines", type=y_type,
                           whichCat=whichCat, dfY = dfY, nChains = 2, 
                           totalScans = totalScans)
    }

    if (dfT == "GP") {
      modT = "GP"
      dfT = NULL
      PostT = GPTreatmentMCMC(t=t, x=x, whichCat = whichCat, type="binary",
                            nScans=nScans, nBurn=nBurn, thin=thin, band=band)
      WAICT = WAICtreatment(t=t, x=x, Post=PostT, modT="GP", type="binary",
                            whichCat=whichCat, dfT = NULL, nChains = 2, 
                            totalScans = totalScans)
    } else {
      modT = "Splines"
      PostT = SplineTreatmentMCMC(t=t, x=x, whichCat = whichCat, type="binary",
                                df=dfT, nScans=nScans, nBurn=nBurn, thin=thin)
      WAICT = WAICtreatment(t=t, x=x, Post=PostT, modT="Splines", type="binary",
                             whichCat=whichCat, dfT = dfT, nChains = 2, 
                            totalScans = totalScans)
    }

    DR = DRmcmcCut(y=y, t=t, x=x, lower=lower, upper=upper,
                   nChains = 2, totalScans = totalScans, whichCat=whichCat,
                   PostT = PostT, PostY = PostY, modY = modY,
                   modT = modT, dfY = dfY, dfT=dfT, nBoot = nBoot,
                   y_type = y_type)
    
    
    l = list(TreatEffect = DR$est,
             TreatEffectSE = DR$se,
             TreatEffectCI = DR$BootQuantile,
             WAICtreatment = WAICT,
             WAICoutcome = WAICY)

  }

  return(l)
}







#' Estimate causal exposure response curve using a doubly robust estimator in high dimensions
#'
#' This function will take in the observed data and estimate the treatment effect curve,
#' which is the average potential outcomes across a range of exposure values.
#' y,x, and t must all be supplied, though all other parameters have pre-set values
#' the user can proceed with unless they wish to change the prior specification.
#'
#' @param y              The outcome to be analyzed
#' @param t              The treatment whose causal effect is to be estimated
#' @param x              An n by p matrix of potential confounders
#' @param locations      The locations for t at which the user wants to estimate E(Y(t)), the average
#'                       potential outcome at level t. The default is a grid across the range of t.
#' @param whichCat       A vector of indices that indicate which variables in x are categorical.
#'                       The default is c(), which means the program defaults to assuming all
#'                       covariates are continuous
#' @param y_type         A categorical variable indicating whether y is binary or continuous.
#'                       Possible values are "binary" or "continuous", and the program defaults
#'                       to "continuous"
#' @param dfY            the degrees of freedom of the splines used to model the relationship between
#'                       the covariates and the outcome. If the user wants to use gaussian process
#'                       priors instead of splines then dfY should be set to "GP"
#' @param dfT            the degrees of freedom of the splines used to model the relationship between
#'                       the covariates and the treatment. If the user wants to use gaussian process
#'                       priors instead of splines then dfY should be set to "GP"
#' @param nScans         The number of MCMC scans to run
#' @param nBurn           The number of MCMC scans that will be dropped as a burn-in
#' @param thin           This number represents how many iterations between each scan
#'                       that is kept
#' @param thetaA         The first parameter of the beta prior on the overall sparsity level
#' @param thetaB         The second parameter of the beta prior on the overall sparsity level
#' @param band           The bandwidth parameter for the gaussian process kernel function
#' @param nBoot          The number of resampling iterations to use when estimating the credible
#'                       intervals
#' @param threshold      The lowest value that the ratio of propensities in the definition of the
#'                       pseudo outcome used in the continuous doubly robust estimator can take.
#'                       This is analagous to trimming the propensity score for binary treatments,
#'                       and the parameter defaults to not trimming.
#'
#'
#' @return A list of values that contain the treatment effect, confidence interval for the
#'         treatment effect, WAIC for the chosen treatment model and outcome model,
#'
#' @export
#' @examples
#'
#' ## p can be larger than n, but we keep the number of covariates small here
#' ## just for illustration so that the code will finish faster
#' n = 200
#' p = 20
#' x = matrix(rnorm(n*p), n, p)
#' t <- 0.6*x[,1] + 0.6*x[,2] + rnorm(n)
#' y <- 5 + 0.05*t^3 - 0.1*t^2 + 0.5*x[,1] + 0.5*x[,2] + rnorm(n)
#'
#' est = DRbayesER(y=y, t=t, x=x, nScans=2000, nBurn=1000, thin=2)

DRbayesER = function(nScans = 20000, nBurn = 10000, thin = 10,
                   y, x, t, locations = seq(quantile(t, .05), 
                                            quantile(t, .95), length=20),
                   whichCat = c(), y_type = "continuous",
                   dfY = 1, dfT = 1, band = 3,
                   thetaA = 1, thetaB = 0.2*dim(x)[2],
                   nBoot=500, threshold=0.00001) {
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  x = scale(x)
  
  tMat = cbind(t, t^2, t^3)
  
  tMatNew = cbind(locations, locations^2, locations^3)
  
  totalScans = floor((nScans - nBurn)/thin)
  
  if (length(unique(t)) <= 2) {
    stop("t must be continuous for this function")
  } else {
    
    if (dfY == "GP") {
      modY = "GP"
      dfY = NULL
      PostY = GPOutcomeMCMC(y=y, tMat=tMat, x=x, whichCat = whichCat, type=y_type,
                            nScans=nScans, nBurn=nBurn, thin=thin, band=band)
      WAICY = WAICoutcome(y=y, x=x, tMat=tMat, Post=PostY, modY="GP", type=y_type,
                          whichCat=whichCat, dfY = NULL, nChains = 2, 
                          totalScans = totalScans)
    } else {
      modY = "Splines"
      PostY = SplineOutcomeMCMC(y=y, tMat=tMat, x=x, whichCat = whichCat, type=y_type,
                                df=dfY, nScans=nScans, nBurn=nBurn, thin=thin)
      WAICY = WAICoutcome(y=y, x=x, tMat=tMat, Post=PostY, modY="Splines", type=y_type,
                          whichCat=whichCat, dfY = dfY, nChains = 2, 
                          totalScans = totalScans)
    }
    
    if (dfT == "GP") {
      modT = "GP"
      dfT = NULL
      PostT = GPTreatmentMCMC(t=t, x=x, whichCat = whichCat, type="continuous",
                              nScans=nScans, nBurn=nBurn, thin=thin, band=band)
      WAICT = WAICtreatment(t=t, x=x, Post=PostT, modT="GP", type="continuous",
                            whichCat=whichCat, dfT = NULL, nChains = 2, 
                            totalScans = totalScans)
    } else {
      modT = "Splines"
      PostT = SplineTreatmentMCMC(t=t, x=x, whichCat = whichCat, type="continuous",
                                  df=dfT, nScans=nScans, nBurn=nBurn, thin=thin)
      WAICT = WAICtreatment(t=t, x=x, Post=PostT, modT="Splines", type="continuous",
                            whichCat=whichCat, dfT = dfT, nChains = 2, 
                            totalScans = totalScans)
    }
    
    DR = DRmcmcContinuousCut(y=y, t=t, tMat=tMat, x=x, tMatNew=tMatNew,
                             nChains = 2, totalScans = totalScans, whichCat=whichCat,
                             PostT = PostT, PostY = PostY, modY = modY,
                             modT = modT, dfY = dfY, dfT=dfT, nBoot = nBoot, 
                             threshold = threshold, y_type = y_type)
    
    
    l = list(TreatEffect = DR$est,
             TreatEffectSE = DR$se,
             TreatEffectCI = cbind(DR$CIlower, DR$CIupper),
             WAICtreatment = WAICT,
             WAICoutcome = WAICY)
    
  }
  
  return(l)
}
