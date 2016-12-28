#'Pseudo-likelihood computation
#'
#'Compute the loglikelihood of a match between a patient database and an observation, 
#'both of which having the same features measured.
#'
#'@param Amat n x K0 matrix the database into which a match is looked for.
#'
#'@param Bvec a vector of length K0 the observation to be matched.
#'
#'@param eps_p a vector of length \code{K} giving the prior discrepancy rate 
#'expected from A to B for the positives, for each variable.
#'
#'@param eps_n a vector of length \code{K} giving the prior discrepancy rate 
#'expected from A to B for the negatives, for each variable.
#'
#'@param piA a vector of length \code{K} giving the prior probabilities of 
#'observing each variable in A.
#'
#'@param piB a vector of length \code{K} giving the prior probabilities of 
#'observing each variable in B.
#'
#'@seealso loglikC
#'
#'@export

loglik_FUN <- function(Amat, Bvec, eps_p, eps_n, piA, piB){
  #Amat = as.matrix(Amat)
  #Bvec = c(Bvec)
  if(ncol(Amat)>1){
    Amat = t(Amat) # K0 x n matrix
  }
  colSums(Amat*Bvec*log((1-eps_n)/piB) + (1-Amat)*(1-Bvec)*log((1-eps_p)/(1-piB)) + Amat*(1-Bvec)*log(eps_n/(1-piB)) + (1-Amat)*Bvec*log(eps_p/piB))
}
