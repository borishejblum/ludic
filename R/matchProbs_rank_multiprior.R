#'Compute the posterior probabilities for an observation to match other observations 
#'updating a matrix of prior matching probabilities
#'
#'For a given observation of index id0 (in \code{1:n}), all the posterior 
#'matching probabilities are computed for the complementary \code{p} 
#'observations using a  matrix of prior matching probabilities previously calculated
#'
#'@param id0 the index for the reference observation in \code{1:n} to be match to p other observations
#'@param computed_dist an \code{n x p} matrix of computed distances used for ranking.
#'@param match_priorprobs an \code{n x p} matrix of prior matching probabilities.
#'@param yes_aggregate logical flag for wether the results are sorted according to the
#'original probability or the aggregated probability. Default is \code{TRUE}.
#'@param k0 number of top (highest) matching probabilities to return. 
#'Default is 5.
#'
#'@return a character vector of length \code{k0*4} containing for each k0
#'potential match, inturn, the probability, the aggregated probability, the logLR 
#'abd the corresponding potential matching ID.
#'
#'@export
matchProbs_rank_multiprior <- function(id0, computed_dist, match_priorprobs, yes_aggregate=FALSE, k0=5){
  dist_sorted <- sort(computed_dist[id0,])
  id_sorted <- names(dist_sorted)
  #if(yes.match){
  #  dis.true = tmp.dis[colnames(dist_all)==temp_id]; 
  #  tmpout = c(sum(tmp.dis > dis.true), sum(tmp.dis >= dis.true))
  #}
  logexptrick_const <- max(dist_sorted)
  prob0 <- exp(dist_sorted + log(match_priorprobs[id0,]) - logexptrick_const)/(exp(-logexptrick_const) + sum(exp(dist_sorted - logexptrick_const + log(match_priorprobs[id0,]))))
  prob1 <- tapply(X=prob0, INDEX=dist_sorted, FUN=sum)
  prob1 <- prob1[match(dist_sorted, unique(dist_sorted))]
  res <- cbind(round(cbind("prob"=prob0, "prob_agg"=prob1, "logLR"=dist_sorted),3),
               "ID"=id_sorted)
  if(yes_aggregate){
    res <- res[order(prob1, decreasing=TRUE),]
  }else{
    res <- res[order(prob0, decreasing=TRUE),]
  }
  return(c(res[1:k0,]))
}
