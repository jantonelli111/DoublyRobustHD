#' Estimate causal effect using doubly robust estimators in high-dimensions
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
#' est = DRbayes(y=y, t=t, x=x, nScans=200, nBurn=100, thin=1)

DRbayes = function(nScans = 20000, nBurn = 10000, thin = 10,
                   y, x, t, whichCat = c(), y_type = "continuous",
                   dfY = 1, dfT = 1, band = 3,
                   thetaA = 1, thetaB = 0.2*dim(x)[2]) {

  n = dim(x)[1]
  p = dim(x)[2]

  x = scale(x)

  if (length(unique(t)) > 2) {
    stop("t must be binary for this function")
  } else {

    if (dfY == "GP") {
      PostY = GPOutcomeMCMC(y=y, tMat=t, x=x, whichCat = whichCat, type=y_type,
                            nScans=nScans, nBurn=nBurn, thin=thin, band=band)
    } else {
      PostY = SplineOutcomeMCMC(y=y, tMat=t, x=x, whichCat = whichCat, type=y_type,
                                df=dfY, nScans=nScans, nBurn=nBurn, thin=thin)
    }

    if (dfT == "GP") {
      PostT = GPTreatmentMCMC(t=t, x=x, whichCat = whichCat, type="binary",
                            nScans=nScans, nBurn=nBurn, thin=thin, band=band)
    } else {
      PostT = SplineTreatmentMCMC(t=t, x=x, whichCat = whichCat, type="binary",
                                df=dfT, nScans=nScans, nBurn=nBurn, thin=thin)
    }

    l = list(TreatEffect = apply(PostY$beta, 3, mean)[2])

  }

  return(l)
}
