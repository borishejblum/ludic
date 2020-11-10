#' Association testing using Han & Lahiri estimating equations and jackknife approach
#'
#'@param match_prob matching probabilities matrix (e.g. obtained through \code{\link{recordLink}}) of 
#'dimensions \code{n1 x n2}.
#'
#'@param y response variable of length \code{n1}. Only binary or gaussian phenotypes 
#'are supported at the moment.
#'
#'@param x a \code{matrix} or a \code{data.frame} of predictors of dimensions \code{n2 x p}. 
#'An intercept is automatically added within the function.
#'
#'@param jackknife_nrep the number of jackknife repetitions.
#'Default is 100 (from Han et al.).
#'
#'@param jackknife_blocksize the number of observations to remove in each jackknife.
#'
#'@param dist_family a character string indicating the distribution family for the glm. 
#'Currently, only \code{'gaussian'} and  \code{'binomial'} are supported. Default 
#'is \code{'gaussian'}.
#'
#'@param methods a character vector which must be a subset of \code{("F", "M", "M2")} 
#'indicating which estimator from Han et al. 2018 should be computed. Default is all 3.
#'
#'@importFrom rootSolve multiroot
#'@importFrom stats na.omit rnorm as.formula model.matrix
#'
#'@return a list containing the following for each estimator in \code{methods}:
#'\itemize{
#'   \item \code{beta} a vector containing the \code{p} estimated coefficients
#'   \item \code{varcov} the \code{p x p} variance-covariance \code{matrix} of the \code{beta} coefficients
#'   \item \code{zscores} a vector containing the \code{p} Z-scores
#'   \item \code{pval} the corresponding Gaussian assumption p-values
#'}
#'
#'@references Han, Y., and Lahiri, P. (2019) Statistical Analysis with Linked Data. 
#'International Statistical Review, 87: S139– S157. https://doi.org/10.1111/insr.12295. 
#'
#'@export
#'
#'@examples
#'#rm(list=ls())
#'res <- list()
#'n_sims <- 1#500
#'for(n in 1:n_sims){
#'x <- matrix(ncol=2, nrow=99, stats::rnorm(n=99*2))
#'
#'#plot(density(rbeta(n=1000, 1,2)))
#'match_prob <- matrix(rbeta(n=103*99, 1, 2), nrow=103, ncol=99)
#'
#'
#'y <- rnorm(n=103, 1, 0.5)
#'res[[n]] <- test_han2018(match_prob, y, x, jackknife_blocksize = 10)
#'#y <- rbinom(n=103, 1, prob=0.5)
#'#res[[n]] <- test_han2018(match_prob, y, x, dist_family = "binomial")
#'cat(n, "/", n_sims, "\n", sep="")
#'}
#'pvals_F <- sapply(lapply(res, "[[", "F"), "[[", "pval")
#'pvals_M <- sapply(lapply(res, "[[", "M"), "[[", "pval")
#'pvals_M2 <- sapply(lapply(res, "[[", "M2"), "[[", "pval")
#'rowMeans(pvals_F<0.05)
#'rowMeans(pvals_M<0.05)
#'rowMeans(pvals_M2<0.05)
#'

test_han2018 <- function(match_prob, y, x,
                         jackknife_nrep = 100, 
                         jackknife_blocksize = max(floor(min(length(y), nrow(x))/jackknife_nrep), 1) ,
                         methods = c("F", "M", "M2"), 
                         dist_family = c("gaussian", "binomial")
){
  # sanity checks
  stopifnot(is.matrix(match_prob))
  stopifnot(is.vector(y))
  
  if(length(which(is.na(y))) > 0){
    warning("y contains NA/nan: to be able to provide results associated observations will be removed")
    y_toremove <- which(is.na(y))
    y <- y[-y_toremove]
    match_prob <- match_prob[-y_toremove, ]
  }
  if(length(which(is.na(x))) > 0){
    warning("x contains NA/nan: to be able to provide results associated observations will be removed")
    x_toremove <- unique(which(is.na(x), arr.ind = TRUE)[, "row"])
    x <- x[-x_toremove, ]
    match_prob <- match_prob[, -x_toremove]
  }
  
  if(is.data.frame(x)){
    x <- stats::model.matrix(stats::as.formula(paste0("~", paste(colnames(x), collapse=" + "))), data = x)
  }else if(is.matrix(x)){
    warning("'x' is a matrix, it is assume to be as would be output from stats::model.matrix()")
  }else{
    stop("x is neither a data.frame nor a matrix")
  }
  
  
  n1 <- length(y)
  n2 <- nrow(x)
  p <- ncol(x)
  stopifnot(nrow(match_prob) == n1)
  stopifnot(ncol(match_prob) == n2)
  
  n1 <- length(y)
  stopifnot(nrow(match_prob) == n1)
  
  if(length(dist_family) > 1){
    dist_family <-  dist_family[1]
  }
  if(!(dist_family %in% c("gaussian", "binomial"))){
    stop("'gaussian' or 'binomial' are the only valid values for dist_family currently supported")
  }
  
  
  if(dist_family == "binomial"){ #logistic regression
    est_eq <- function(beta, x, y, match_prob){
      
      tcrossprod(t(x), match_prob) %*% (y - match_prob %*% expit(x %*% beta))
      
    }
  }else if(dist_family == "gaussian"){  #normal regression
    
    est_eq <- function(beta, x, y, match_prob){
      
      tcrossprod(t(x), match_prob) %*% (y - match_prob %*% x %*% beta)
      
    }
    
  }else{
    stop("'gaussian' or 'binomial' are the only valid values for dist_family currently supported")
  }
  
  #samples2jackknife <- cbind(sample(n1, size=jackknife_nrep, replace=TRUE), 
  #                           sample(n2, size=jackknife_nrep, replace=TRUE))
  samples2jackknife <- lapply(seq_len(jackknife_nrep), function(x){
    cbind(sample(n1, jackknife_blocksize, replace=FALSE),
          sample(n2, jackknife_blocksize, replace=FALSE)
    )}
  )
  
  
  
  res <- list()
  
  
  
  if("F" %in% methods){
    betaF <- rootSolve::multiroot(f = function(b){est_eq(beta = b, x = x, y = y, match_prob = match_prob)},
                                  start = rep(0, times = p))$root
    betaF_b <- matrix(NA, nrow = jackknife_nrep, ncol = p)
    for(j in 1:jackknife_nrep){
      betaF_b[j, ] <- rootSolve::multiroot(f = function(b){est_eq(beta = b, 
                                                                  x = x[-samples2jackknife[[j]][, 2], , drop = FALSE], 
                                                                  y = y[-samples2jackknife[[j]][, 1]], 
                                                                  match_prob = match_prob[-samples2jackknife[[j]][, 1],
                                                                                          -samples2jackknife[[j]][, 2]])},
                                           start = rep(0, times = p))$root
    }
    errF <- betaF_b - matrix(colMeans(betaF_b), nrow= jackknife_nrep, ncol=p, byrow = TRUE)
    varcovF <- (jackknife_nrep-1)/jackknife_nrep * crossprod(errF)
    zscoreF <- betaF/sqrt(diag(varcovF))
    
    res[["F"]] <- list("beta" = betaF, "varcov" = varcovF, 
                       "zscore" = zscoreF, "pval" = 2*(1-pnorm(abs(zscoreF))))
    
  }
  
  
  
  if("M" %in% methods){
    match_probM <- matrix(0, nrow = n1, ncol = n2)
    imax <- max.col(match_prob)
    for(i in 1:n1){
      match_probM[i, imax[i]] <- match_prob[i, imax[i]] 
    }
    betaM <- rootSolve::multiroot(f = function(b){est_eq(beta = b, x = x, y = y, match_prob = match_probM)},
                                  start = rep(0, times = p))$root
    betaM_b <- matrix(NA, nrow = jackknife_nrep, ncol = p)
    for(j in seq_len(jackknife_nrep)){
      match_probM <- matrix(0, nrow = n1-jackknife_blocksize, ncol = n2-jackknife_blocksize)
      match_prob_temp = match_prob[-samples2jackknife[[j]][, 1],
                                   -samples2jackknife[[j]][, 2]]
      imax <- max.col(match_prob_temp)
      for(i in seq_len(nrow(match_prob_temp))){
        match_probM[i, imax[i]] <- match_prob_temp[i, imax[i]] 
      }
      betaM_b[j, ] <- rootSolve::multiroot(f = function(b){est_eq(beta = b,
                                                                  x = x[-samples2jackknife[[j]][, 2], , drop = FALSE],
                                                                  y = y[-samples2jackknife[[j]][, 1]], 
                                                                  match_prob = match_probM)},
                                           start = rep(0, times = p))$root
    }
    errM <- betaM_b - matrix(colMeans(betaM_b), nrow= jackknife_nrep, ncol=p, byrow = TRUE)
    varcovM <- (jackknife_nrep-1)/jackknife_nrep * crossprod(errM)
    zscoreM <- betaM/sqrt(diag(varcovM))
    
    res[["M"]] <- list("beta" = betaM, "varcov" = varcovM, 
                       "zscore" = zscoreM, "pval" = 2*(1-pnorm(abs(zscoreM))))
  }
  
  
  
  if("M2" %in% methods){
    
    match_probM2 <- matrix(0, nrow = n1, ncol = n2)
    imax <- max.col(match_prob)
    match_prob_bis <- match_prob
    match_prob_bis[cbind(seq_len(n1), imax)] <- -Inf
    imax2 <- max.col(match_prob_bis)
    for(i in seq_len(n1)){
      match_probM2[i, imax[i]] <- match_prob[i, imax[i]]
      match_probM2[i, imax2[i]] <- match_prob[i, imax2[i]]
    }
    betaM2 <- rootSolve::multiroot(f = function(b){est_eq(beta = b, x = x, y = y, match_prob = match_probM2)},
                                   start = rep(0, times = p))$root
    betaM2_b <- matrix(NA, nrow = jackknife_nrep, ncol = p)
    
    for(j in 1:jackknife_nrep){
      match_probM2 <- matrix(0, nrow = n1-jackknife_blocksize, ncol = n2-jackknife_blocksize)
      match_prob_temp = match_prob[-samples2jackknife[[j]][, 1],
                                   -samples2jackknife[[j]][, 2]]
      imax <- max.col(match_prob_temp)
      match_prob_temp_bis <- match_prob_temp
      match_prob_temp_bis[cbind(seq_len(nrow(match_prob_temp)), imax)] <- -Inf
      imax2 <- max.col(match_prob_temp_bis)
      for(i in seq_len(nrow(match_prob_temp))){
        match_probM2[i, imax[i]] <- match_prob_temp[i, imax[i]]
        match_probM2[i, imax2[i]] <- match_prob_temp[i, imax2[i]]
      }
      betaM2_b[j, ] <- rootSolve::multiroot(f = function(b){est_eq(beta = b, 
                                                                   x = x[-samples2jackknife[[j]][, 2], , drop = FALSE], 
                                                                   y = y[-samples2jackknife[[j]][, 1]], 
                                                                   match_prob = match_probM2)},
                                            start = rep(0, times = p))$root
    }
    errM2 <- betaM2_b - matrix(colMeans(betaM2_b), nrow= jackknife_nrep, ncol=p, byrow = TRUE)
    varcovM2 <- (jackknife_nrep-1)/jackknife_nrep * crossprod(errM2)
    zscoreM2 <- betaM2/sqrt(diag(varcovM2))
    
    res[["M2"]] <- list("beta" = betaM2, "varcov" = varcovM2, 
                        "zscore" = zscoreM2, "pval" = 2*(1-pnorm(abs(zscoreM2))))
  }
  
  
  return(res)
  
}


