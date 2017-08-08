#' Testing the combination of several matching thresholds
#' 
#' Computes p-value for each considered threshold, and computes a p-value 
#' associated with their combination through Fisher's method using perturbation resampling.
#'
#'@param match_prob matching probabilities matrix (e.g. obtained through \code{\link{recordLink}}) of 
#'dimensions \code{n1 x n2}.
#'
#'@param y response variable of length \code{n1}. Only binary phenotypes are supported at the moment.
#'
#'@param x a matrix of predictors of dimensions \code{n2 x p}
#'
#'@param thresholds a vector (possibly of length \code{1}) containing the different threshold 
#'to use to call a match. Default is \code{seq(from = 0.1, to = 0.9, by = 0.2)}.
#'
#'@param nb_perturb the number of perturbation used for the p-value combination.
#'Default is 200.
#'
#'@importFrom landpred VTM
#'@importFrom fGarch dsstd sstdFit
#'@importFrom stats binomial glm na.omit rnorm
#'
#'@return a list containing the following:
#'\itemize{
#'   \item \code{influencefn_pvals} p-values obtained from influence fonction perturbations
#'   with the covariates as columns and the \code{thresholds} as rows, with an additional row 
#'   at the top for the combination 
#'   \item \code{wald_pvals} a matrix containing the p-values obtained from the Wald 
#'   test with the covariates as columns and the \code{thresholds} as rows
#'   \item \code{ptbed_pvals} a list containing, for each covariates, a matrix with
#'   the \code{nb_perturb} perturbed p-values with the different \code{thresholds}
#'   as rows
#'   \item \code{theta_avgImpute} a matrix of the estimated coefficients from the glm when imputing 
#'   the weighted average for covariates (as columns) with the \code{thresholds} as rows
#'   \item \code{sd_theta} a matrix of the estimated SD (from the influence function) of the 
#'   coefficients from the glm when imputing the weighted average for covariates (as columns),
#'   with the \code{thresholds} as rows
#'   \item \code{ptbed_theta_avgImpute} a list containing, for each covariates, a matrix with
#'   the \code{nb_perturb} perturbed estimated coefficients from the glm when imputing 
#'   the weighted average for covariates, with the different \code{thresholds}
#'   as rows
#'}
#'
#'@export
#'
#'@examples
#'res <- list()
#'n_sims <- 1#000
#'for(n in 1:n_sims){
#'y <- rbinom(n=103, 1, prob=0.5)
#'x <- matrix(ncol=2, nrow=99, stats::rnorm(n=99*2))
#'
#'#plot(density(rbeta(n=1000, 1,2)))
#'match_prob <- matrix(rbeta(n=103*99, 1, 2), nrow=103, ncol=99)
#'
#'res[[n]] <- test_combine(match_prob, y, x)$influencefn_pvals
#'cat(n, "/", n_sims, "\n", sep="")
#'}
#'size <- matrix(NA, ncol=nrow(res[[1]]), nrow=ncol(res[[1]])-1)
#'colnames(size) <- rownames(res[[1]])
#'rownames(size) <- colnames(res[[1]])[-ncol(res[[1]])]
#'for(i in 1:(ncol(res[[1]])-1)){
#'  size[i, ] <- rowMeans(sapply(res, function(m){m[, i]<0.05}))
#'}
#'size
#'

test_combine <- function(match_prob, y, x,
                     thresholds = seq(from = 0.1, to = 0.9, by = 0.2), #c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95),
                     nb_perturb = 200){
  
  # sanity checks
  stopifnot(is.matrix(match_prob))
  stopifnot(is.matrix(x))
  stopifnot(is.vector(y))
  
  n1 <- length(y)
  n2 <- nrow(x)
  stopifnot(nrow(match_prob) == n1)
  stopifnot(ncol(match_prob) == n2)
  
  nb_thresholds <- length(thresholds)
  mismatch_avg <- numeric(nb_thresholds)
  names(mismatch_avg) <- thresholds
  
  # initializing results
  betahat_best <- matrix(NA, ncol = 3, nrow = length(thresholds))
  betahat_avg <- matrix(NA, ncol = 3, nrow = length(thresholds))
  rownames(betahat_best) <- as.character(thresholds)
  rownames(betahat_avg) <- as.character(thresholds)
  eta <- list()
  theta_avg <- matrix(NA, ncol = ncol(x) + 1, nrow = length(thresholds))
  rownames(theta_avg) <- as.character(thresholds)
  wald_pvals <- matrix(NA, ncol = ncol(x) + 1, nrow = length(thresholds))
  rownames(wald_pvals) <- as.character(thresholds)
  
  for(i in 1:length(thresholds)){
    cut_p <- thresholds[i]
    prob_sup_cut <- (match_prob > cut_p)
    
    #construct the data-frame for the glm
    xi <- rowSums(prob_sup_cut) > 0
    n_rho <- sum(xi)
    y_match <- y*xi
    match_prob_sel <-  diag(1*xi) %*% match_prob
    x_best <- x[max.col(match_prob_sel), ]
    x_impute <- match_prob_sel %*% x / rowSums(match_prob_sel) #diag(1/rowSums(match_prob_sel)) %*% match_prob_sel %*% x/rowSums(match_prob_sel)
    
    y_match_sub <- y[xi]
    x_best_sub <- x[max.col(match_prob[xi, ]), ]
    x_impute_sub <- match_prob[xi, ]%*%x/rowSums(match_prob[xi, ])
    
    temp_1 <- try(summary(stats::glm(y_match ~ x_best, family = stats::binomial, na.action = stats::na.omit))$coef[2,-3], silent = TRUE)
    if(inherits(temp_1, "try-error")){temp_1 <- rep(NA, 3)}
    temp_avg <- try(summary(stats::glm(y_match ~ x_impute, family = stats::binomial, na.action = stats::na.omit))$coef[2,-3], silent = TRUE)
    if(inherits(temp_avg, "try-error")){temp_avg <- rep(NA, 3)}
    
    theta <- summary(stats::glm(y_match ~ x_impute, family = stats::binomial, na.action = stats::na.omit))$coef[,"Estimate", drop=FALSE]
    wald_pvals[i, ] <- summary(stats::glm(y_match ~ x_impute, family = stats::binomial, na.action = stats::na.omit))$coef[,"Pr(>|z|)", drop=FALSE]
    #summary(stats::glm(y_match_sub ~ x_impute_sub, family = stats::binomial, na.action = stats::na.omit))$coef
    #summary(stats::glm(y_match ~ x_impute, family = stats::binomial, na.action = stats::na.omit))$coef
  
    # Z_sub <- stats::model.matrix( ~ x_impute_sub)
    # I_rho <- 1/n_rho*crossprod(apply(Z_sub, 2, function(colu){colu*sqrt(expit_dev1(Z_sub%*%theta))}))
    # eta[[as.character(cut_p)]] <- 1/n_rho*solve(I_rho)%*%t(Z_sub)%*%diag(x=(y_match_sub - expit(Z_sub%*%theta)[, "Estimate"]))
    # sqrt(apply(eta[[as.character(cut_p)]], 1, crossprod))
    
    x_impute_noNA <-  x_impute 
    x_impute_noNA[is.na(x_impute[, 1])] <- 0
    Z <- diag(xi) %*% stats::model.matrix( ~ x_impute_noNA)
    I_rho <- 1/n_rho*crossprod(apply(Z, 2, function(colu){colu*sqrt(expit_dev1(Z%*%theta))}))
    eta[[as.character(cut_p)]] <- 1/n_rho*solve(I_rho)%*%t(Z)%*%diag(x=(y_match - xi*expit(Z%*%theta)[, "Estimate"]))

    # fit <- stats::glm(y_match ~ x_impute, family = stats::binomial, na.action = stats::na.omit)
    # summary(fit)$coef
    # sqrt(apply(eta[[as.character(cut_p)]], 1, crossprod))
    # solve(vcov(fit))/n_rho
    
    theta_avg[i,] <- summary(stats::glm(y_match ~ x_impute, family = stats::binomial))$coef[,"Estimate", drop=FALSE]
  }
  colnames(theta_avg) <- rownames(theta)
  colnames(wald_pvals) <- rownames(theta)
  colnames(betahat_best) <- names(temp_1)
  colnames(betahat_avg) <- names(temp_avg)
  betahat_avg <- cbind(betahat_avg, thresholds)
  betahat_best <- cbind(betahat_best, thresholds) 
  
  sigma <- list()
  sds <- matrix(ncol=ncol(theta_avg), nrow=length(thresholds))
  colnames(sds) <- colnames(theta_avg)
  rownames(sds) <- thresholds
  eta_mat <- list()
  for(j in 1:ncol(theta_avg)){
    eta_mat[[j]] <- sapply(X=eta, FUN=function(m){m[j,]})
    sigma[[colnames(theta_avg)[j]]] <- crossprod(eta_mat[[j]])  
    sds[, j] <- sqrt(diag(sigma[[colnames(theta_avg)[j]]]))
  }
  names(eta_mat) <- colnames(theta_avg)
  
  
  # Combine p-values
  pvals <- pval_zscore(theta_avg, sds)
  fisher_comb <- apply(pvals, MARGIN = 2, FUN = comb_pvals)
  
  B <- nb_perturb #number of perturbations
  theta_avg_star <- lapply(eta_mat, function(m){crossprod(m, matrix(stats::rnorm(n1*B, mean = 0, sd = 1), nrow=n1, ncol=B))})
  pvals_star <- list()
  for(k in 1:ncol(sds)){
    pvals_star[[k]] <- apply(theta_avg_star[[k]], MARGIN=2, FUN=pval_zscore, sigma=sds[, k])
  }
  fisher_comb_star <- lapply(pvals_star, function(m){apply(m, MARGIN=2, FUN=comb_pvals)})
  pval_combined <- matrix(nrow=1, ncol=ncol(pvals))
  row.names(pval_combined) <- "Combined p-value"
  colnames(pval_combined) <- colnames(pvals)
  for(k in 1:ncol(sds)){
    pval_combined[1, k] <- 1-sum(fisher_comb[k] >= fisher_comb_star[[k]])/B
  }

  return(list("influencefn_pvals" = cbind.data.frame(rbind(pval_combined, pvals), "matching_threshold"=c("combined", rownames(pvals))),
              "wald_pvals" = wald_pvals,
              "ptbed_pvals" = pvals_star,
              "theta_avgImpute" = theta_avg,
              "sd_theta"=sds,
              "ptbed_theta_avgImpute" = theta_avg_star)
  )
  
}


