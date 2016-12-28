#'Pseudo-likelihood computation
#'
#'compute the loglikelihood of a match between a patient database and an observation, 
#'both of which having the same features measured.
#'
#'@param Amat n x K0 matrix the database into which a match is looked for.
#'
#'@param Bvec a vector of length K0 the observation to be matched.
#'
#'@param eps_p
#'
#'@param eps_n
#'
#'@param piA
#'
#'@param piB
#'
#'@export

loglik_cont2diff <- function(Amat, Bmat, d_max){
  Diffmat <- array(NA, dim=c(nrow(Amat), nrow(Bmat), ncol(Amat)))
  for(i in 1:nrow(Amat)){
    for(j in 1:nrow(Bmat)){
      Diffmat[i,j,] <- Amat[i,]-Bmat[j,]
    }
  }
  browser()
  Diffmat_bin <- 1*(Diffmat <d_max)
  pi_same_EB <- colMeans()
  
  if(ncol(Amat)>1){
    Amat = t(Amat) # K0 x n matrix
  }
  colSums(Amat*Bvec*log((1-eps_n)/piB) + (1-Amat)*(1-Bvec)*log((1-eps_p)/(1-piB)) + Amat*(1-Bvec)*log(eps_n/(1-piB)) + (1-Amat)*Bvec*log(eps_p/piB))
}
