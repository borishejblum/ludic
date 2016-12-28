#' Compute pvalues from Z-score
#' 
#' @param beta the estimate
#' 
#' @param sigma estimate's estimated variance
#' 
#' @return the p-value
#' 
#' @importFrom stats pnorm
#' 
#' @export

pval_zscore <- function(beta, sigma){
  2*stats::pnorm(-abs(beta/sigma))
}