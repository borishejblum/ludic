#'Compute the probabilities for an observation to match other observations
#'
#'For a given observation of index id0 (in \code{1:n}), all the matching probabilities are computed
#'for the other \code{p} observations.
#'
#'@param id0 the index for the reference observation in \code{1:n} to be match to p other observations
#'@param computed_dist an \code{n x p} matrix of computed distances used for ranking.
#'@param prop_match estimated proportion of match ("rho_1")
#'@param yes_aggregate logical flag for wether the results are sorted according to the
#'original probability or the aggregated probability. Default is \code{TRUE}. Used only by \code{matchProbs_rank}.
#'@param k0 number of top (highest) matching probabilities to return. 
#'Default is 5. Used only by \code{matchProbs_rank}.
#'
#'@return \itemize{
#'\item \code{matchProbs_rank} returns a vector of length \code{k0*4} containing for each k0
#'potential match in turn: the probability, the aggregated probability, the logLR 
#'and the corresponding potential matching ID.
#'\item \code{matchProbs_rank_complete} returns a vector of length \code{p} containing the matching 
#'probabilities to each of all \code{p} observations in their original order.
#'}
#'
#'@seealso matchProbs_rank_full_C
#'
#'@export
matchProbs_rank <- function(id0, computed_dist, prop_match, yes_aggregate=TRUE, k0=5){
  dist_sorted <- sort(computed_dist[id0,])
  id_sorted <- names(dist_sorted)
  #if(yes.match){
  #  dis.true = tmp.dis[colnames(dist_all)==temp_id]; 
  #  tmpout = c(sum(tmp.dis > dis.true), sum(tmp.dis >= dis.true))
  #}
  logexptrick_const <- max(dist_sorted)
  prob0 <- exp(dist_sorted + log(prop_match) - logexptrick_const)/(exp(-logexptrick_const) + sum(exp(dist_sorted - logexptrick_const + log(prop_match))))
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

#'@rdname matchProbs_rank
#'@export
matchProbs_rank_complete <- function(id0, computed_dist, prop_match){
  logexptrick_const <- max(computed_dist[id0,])
  prob0 <- exp(computed_dist[id0,] + log(prop_match) - logexptrick_const)/(exp(-logexptrick_const) + sum(exp(computed_dist[id0,] - logexptrick_const + log(prop_match))))
  return(prob0)
}
