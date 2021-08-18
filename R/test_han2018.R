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
#'@param covar_y a \code{matrix} or a \code{data.frame} of predictors of dimensions \code{n1 x q1}. 
#'An intercept is automatically added within the function.
#'
#'@param covar_x a \code{matrix} or a \code{data.frame} of predictors of dimensions \code{n2 x q2}. 
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
#'International Statistical Review, 87: S139– S157. 
#'\doi{10.1111/insr.12295}. 
#'
#'@export
#'
#'@examples
#'# rm(list=ls())
#'# n_sims <- 500
#'# res <- pbapply::pblapply(1:n_sims, function(n){
#'# nx <- 99
#'# ny <- 103
#'# x <- matrix(ncol=2, nrow=ny, stats::rnorm(n=ny*2))
#'# 
#'# #plot(density(rbeta(n=1000, 1,2)))
#'# match_prob <- diag(ny)[, 1:nx]#matrix(rbeta(n=ny*nx, 1, 2), nrow=ny, ncol=99)
#'# 
#'# covar_y <- matrix(rnorm(n=ny, 1, 0.5), ncol=1)
#'# covar_x <- matrix(ncol=3, nrow=ny, stats::rnorm(n=ny*3))
#'# 
#'# #y <- rnorm(n=ny, mean = x %*% c(2,-3)  + covar_x %*% rep(0.2, ncol(covar_x)) + 0.5*covar_y, 0.5)
#'# y <- rbinom(n=ny, 1, prob=expit(x %*% c(2,-3)  + covar_x %*% 
#'#             rep(0.2, ncol(covar_x)) + 0.5*covar_y))
#'# #glm(y~0+x+covar_y+covar_x, family = "binomial")
#'# return(
#'# #test_han2018(match_prob, y, x, jackknife_blocksize = 10, covar_x = NULL, covar_y = NULL)
#'# test_han2018(match_prob, y[1:ny], x[1:nx, ], dist_family = "binomial", 
#'#              jackknife_blocksize = 10, covar_x = covar_x[1:nx, ], 
#'#              covar_y = covar_y[1:ny, , drop=FALSE])
#'# )
#'# }, cl=parallel::detectCores()-1)
#'# pvals_F <- sapply(lapply(res, "[[", "F"), "[[", "beta")
#'# pvals_M <- sapply(lapply(res, "[[", "M"), "[[", "beta")
#'# pvals_M2 <- sapply(lapply(res, "[[", "M2"), "[[", "beta")
#'# quantile(pvals_F)
#'# quantile(pvals_M)
#'# quantile(pvals_M2)
#'# rowMeans(pvals_F<0.05)
#'# rowMeans(pvals_M<0.05)
#'# rowMeans(pvals_M2<0.05)
#'

test_han2018 <- function(match_prob, y, x, covar_y = NULL, covar_x = NULL,
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
  
  if(length(which(is.na(covar_x))) > 0){
    warning("covar_x contains NA/nan: to be able to provide results associated observations will be removed")
    covar_x_toremove <- unique(which(is.na(covar_x), arr.ind = TRUE)[, "row"])
    covar_x <- covar_x[-covar_x_toremove, ]
    match_prob <- match_prob[, -covar_x_toremove]
  }
  
  if(length(which(is.na(covar_y))) > 0){
    warning("covar_x contains NA/nan: to be able to provide results associated observations will be removed")
    covar_y_toremove <- unique(which(is.na(covar_y), arr.ind = TRUE)[, "row"])
    covar_y <- covar_y[-covar_y_toremove, ]
    match_prob <- match_prob[-covar_y_toremove, ]
  }
  
  if(is.data.frame(x)){
    #no intercept according to Han estimating equations...
    x <- stats::model.matrix(stats::as.formula(paste0("~ 0 +", paste(colnames(x), collapse=" + "))), data = x)
  }else if(is.matrix(x)){
    warning("'x' is a matrix, it is assume to be as would be the output from stats::model.matrix() without the intercept. If the intercept is found present, it is removed.")
    if(all(x[,1]==1) | ifelse(!is.null(colnames(x)), colnames(x)[1] == "(Intercept)", FALSE)){
      x <- x[, -1, drop = FALSE]
    }
  }else{
    stop("x is neither a data.frame nor a matrix")
  }
  
  if(!is.null(covar_x) && is.data.frame(covar_x)){
    covar_x <- stats::model.matrix(stats::as.formula(paste0("~ 0 +", paste(colnames(covar_x), collapse=" + "))), data = covar_x)
  }else if(is.matrix(covar_x)){
    warning("'covar_x' is a matrix, it is assume to be as would be the output from stats::model.matrix() without the intercept. If the intercept is found present, it is removed.")
    if(all(covar_x[,1]==1) | ifelse(!is.null(colnames(covar_x)), colnames(covar_x)[1] == "(Intercept)", FALSE)){
      covar_x <- covar_x[, -1, drop = FALSE]
    }
  }else if(!is.null(covar_x)){
    stop("covar_x is neither a data.frame nor a matrix")
  }
  
  if(!is.null(covar_y) && is.data.frame(covar_y)){
    covar_y <- stats::model.matrix(stats::as.formula(paste0("~ 0 +", paste(colnames(covar_y), collapse=" + "))), data = covar_y)
  }else if(is.matrix(covar_y)){
    warning("'covar_y' is a matrix, it is assume to be as would be the output from stats::model.matrix() without the intercept. If the intercept is found present, it is removed.")
    if(all(covar_y[,1]==1) | ifelse(!is.null(colnames(covar_y)), colnames(covar_y)[1] == "(Intercept)", FALSE)){
      covar_y <- covar_y[, -1, drop = FALSE]
    }
  }else if(!is.null(covar_y)){
    stop("covar_y is neither a data.frame nor a matrix")
  }
  
  
  n1 <- length(y)
  n2 <- nrow(x)
  p <- ncol(x)
  stopifnot(nrow(match_prob) == n1)
  stopifnot(ncol(match_prob) == n2)
  if(!is.null(covar_x)){
    stopifnot(nrow(covar_x) == n2)
  }
  if(!is.null(covar_y)){
    stopifnot(nrow(covar_y) == n1)
  }
  n1 <- length(y)
  stopifnot(nrow(match_prob) == n1)
  
  if(length(dist_family) > 1){
    dist_family <-  dist_family[1]
  }
  if(!(dist_family %in% c("gaussian", "binomial"))){
    stop("'gaussian' or 'binomial' are the only valid values for dist_family currently supported")
  }
  
  
  if(dist_family == "binomial"){ #logistic regression
    if(!is.null(covar_x) && !is.null(covar_y)){
      ptot <- p + ncol(covar_x) + ncol(covar_y)
      est_eq <- function(beta, x, y, covar_x=NULL, covar_y=NULL, match_prob){
        xx <- cbind(x, covar_x)
        H <- rbind(tcrossprod(t(xx), match_prob), t(covar_y))
        nxx <- ncol(xx)
        ncy <- ncol(covar_y)
        return(H %*% (y - expit(match_prob %*% xx %*% beta[1:nxx] + covar_y %*% beta[nxx + 1:ncy])))
      }
    }else if(!is.null(covar_x)){
      ptot <- p + ncol(covar_x)
      est_eq <- function(beta, x, y, covar_x=NULL, covar_y=NULL, match_prob){
        xx <- cbind(x, covar_x)
        H <- tcrossprod(t(xx), match_prob)
        return(H %*% (y - expit(match_prob %*% xx %*% beta)))
      }
    }else if(!is.null(covar_y)){
      ptot <- p + ncol(covar_y)
      est_eq <- function(beta, x, y, covar_x=NULL, covar_y=NULL, match_prob){
        nx <- ncol(x)
        ncy <- ncol(covar_y)
        H <- rbind(tcrossprod(t(x), match_prob), t(covar_y))
        return(H %*% (y - expit(match_prob %*% x %*% beta[1:nx] + covar_y %*% beta[nx + 1:ncy])))
      }
    }else{
      ptot <- p
      est_eq <- function(beta, x, y, covar_x=NULL, covar_y=NULL, match_prob){
        H <- tcrossprod(t(x), match_prob)
        return(H %*% (y - expit(match_prob %*% x %*% beta)))
      }
    }
  }else if(dist_family == "gaussian"){  #linear regression
    if(!is.null(covar_x) && !is.null(covar_y)){
      ptot <- p + ncol(covar_x) + ncol(covar_y)
      est_eq <- function(beta, x, y, covar_x=NULL, covar_y=NULL, match_prob){
        xx <- cbind(x, covar_x)
        H <- rbind(tcrossprod(t(xx), match_prob), t(covar_y))
        nxx <- ncol(xx)
        ncy <- ncol(covar_y)
        return(H %*% (y - (match_prob %*% xx %*% beta[1:nxx] + covar_y %*% beta[nxx + 1:ncy])))
      }
    }else if(!is.null(covar_x)){
      ptot <- p + ncol(covar_x)
      est_eq <- function(beta, x, y, covar_x=NULL, covar_y=NULL, match_prob){
        xx <- cbind(x, covar_x)
        H <- tcrossprod(t(xx), match_prob)
        return(H %*% (y - match_prob %*% xx %*% beta))
      }
    }else if(!is.null(covar_y)){
      ptot <- p + ncol(covar_y)
      est_eq <- function(beta, x, y, covar_x=NULL, covar_y=NULL, match_prob){
        nx <- ncol(x)
        ncy <- ncol(covar_y)
        H <- rbind(tcrossprod(t(x), match_prob), t(covar_y))
        return(H %*% (y - (match_prob %*% x %*% beta[1:nx] + covar_y %*% beta[nx + 1:ncy])))
      }
    }else{
      ptot <- p
      est_eq <- function(beta, x, y, covar_x=NULL, covar_y=NULL, match_prob){
        H <- tcrossprod(t(x), match_prob)
        return(H %*% (y - match_prob %*% x %*% beta))
      }
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
    betaF <- rootSolve::multiroot(f = function(b){est_eq(beta = b, x = x, y = y, 
                                                         covar_x = covar_x, covar_y = covar_y, 
                                                         match_prob = match_prob)},
                                  start = rep(0, times = ptot))$root
    betaF_b <- matrix(NA, nrow = jackknife_nrep, ncol = ptot)
    for(j in 1:jackknife_nrep){
      betaF_b[j, ] <- rootSolve::multiroot(f = function(b){est_eq(beta = b, 
                                                                  x = x[-samples2jackknife[[j]][, 2], , drop = FALSE], 
                                                                  y = y[-samples2jackknife[[j]][, 1]], 
                                                                  covar_x = if(!is.null(covar_x)){covar_x[-samples2jackknife[[j]][, 2], , drop = FALSE]}else{NULL},
                                                                  covar_y = if(!is.null(covar_y)){covar_y[-samples2jackknife[[j]][, 1], , drop = FALSE]}else{NULL},
                                                                  match_prob = match_prob[-samples2jackknife[[j]][, 1],
                                                                                          -samples2jackknife[[j]][, 2]])},
                                           start = rep(0, times = ptot))$root
    }
    errF <- betaF_b - matrix(colMeans(betaF_b), nrow = jackknife_nrep, ncol = ptot, byrow = TRUE)
    varcovF <- (jackknife_nrep-1)/jackknife_nrep * crossprod(errF)
    zscoreF <- betaF/sqrt(diag(varcovF))
    
    res[["F"]] <- list("beta" = betaF[1:p], "varcov" = varcovF[1:p,1:p], 
                       "zscore" = zscoreF[1:p], "pval" = 2*(1-pnorm(abs(zscoreF)))[1:p])
    
  }
  
  
  
  if("M" %in% methods){
    match_probM <- matrix(0, nrow = n1, ncol = n2)
    imax <- max.col(match_prob)
    for(i in 1:n1){
      match_probM[i, imax[i]] <- match_prob[i, imax[i]] 
    }
    betaM <- rootSolve::multiroot(f = function(b){est_eq(beta = b, x = x, y = y, 
                                                         covar_x = covar_x, covar_y = covar_y, 
                                                         match_prob = match_probM)},
                                  start = rep(-0, times = ptot))$root
    betaM_b <- matrix(NA, nrow = jackknife_nrep, ncol = ptot)
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
                                                                  covar_x =if(!is.null(covar_x)){covar_x[-samples2jackknife[[j]][, 2], , drop = FALSE]}else{NULL},
                                                                  covar_y =if(!is.null(covar_y)){covar_y[-samples2jackknife[[j]][, 1], , drop = FALSE]}else{NULL},
                                                                  match_prob = match_probM)},
                                           start = rep(0, times = ptot))$root
    }
    errM <- betaM_b - matrix(colMeans(betaM_b), nrow = jackknife_nrep, ncol = ptot, byrow = TRUE)
    varcovM <- (jackknife_nrep-1)/jackknife_nrep * crossprod(errM)
    zscoreM <- betaM/sqrt(diag(varcovM))
    
    res[["M"]] <- list("beta" = betaM[1:p], "varcov" = varcovM[1:p, 1:p], 
                       "zscore" = zscoreM[1:p], "pval" = 2*(1-pnorm(abs(zscoreM)))[1:p])
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
    betaM2 <- rootSolve::multiroot(f = function(b){est_eq(beta = b, x = x, y = y, 
                                                          covar_x = covar_x, covar_y = covar_y, 
                                                          match_prob = match_probM2)},
                                   start = rep(0, times = ptot))$root
    betaM2_b <- matrix(NA, nrow = jackknife_nrep, ncol = ptot)
    
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
                                                                   covar_x =if(!is.null(covar_x)){covar_x[-samples2jackknife[[j]][, 2], , drop = FALSE]}else{NULL},
                                                                   covar_y =if(!is.null(covar_y)){covar_y[-samples2jackknife[[j]][, 1], , drop = FALSE]}else{NULL},
                                                                   match_prob = match_probM2)},
                                            start = rep(0, times = ptot))$root
    }
    errM2 <- betaM2_b - matrix(colMeans(betaM2_b), nrow = jackknife_nrep, ncol = ptot, byrow = TRUE)
    varcovM2 <- (jackknife_nrep-1)/jackknife_nrep * crossprod(errM2)
    zscoreM2 <- betaM2/sqrt(diag(varcovM2))
    
    res[["M2"]] <- list("beta" = betaM2[1:p], "varcov" = varcovM2[1:p, 1:p], 
                        "zscore" = zscoreM2[1:p], "pval" = 2*(1-pnorm(abs(zscoreM2)))[1:p])
  }
  
  
  return(res)
  
}


