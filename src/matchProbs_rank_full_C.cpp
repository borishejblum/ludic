#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//'Compute the probabilities for an observation to match other observations
//'
//'C++ version: for a given observation of index id0 (in \code{1:n}), all the matching probabilities are computed
//'for the other \code{p} observations.
//'
//'@param id0 the index for the reference observation in \code{1:n} to be match to p other observations
//'@param computed_dist an \code{n x p} matrix of computed distances used for ranking.
//'@param prop_match estimated proportion of match ("rho_1")
//'@param yes_aggregate logical flag for wether the results are sorted according to the
//'original probability or the aggregated probability. Default is \code{TRUE}.
//'@param k0 number of top (highest) matching probabilities to return. 
//'Default is 5.
//'
//'@return a character vector of length \code{k0*4} containing for each k0
//'potential match, inturn, the probability, the aggregated probability, the logLR 
//'abd the corresponding potential matching ID.
//'
//'@export
// [[Rcpp::export]]
NumericMatrix matchProbs_rank_full_C(NumericMatrix computed_dist, 
                                     double prop_match){
  
  mat compD = as<mat>(computed_dist);
  int p = compD.n_cols;
  int n = compD.n_rows;
  
  vec logexptrick_const = max(compD, 1);
  mat prob0 = mat(n, p);
  
  
  for(int i=0; i<n; i++){
    double logexptrick_const_rowi = logexptrick_const(i);
    rowvec compD_rowi = compD.row(i);
    double normconst = exp(-logexptrick_const_rowi) + sum(exp(compD_rowi - logexptrick_const_rowi + log(prop_match)));
    prob0.row(i) = exp(compD_rowi + log(prop_match) - logexptrick_const_rowi) / normconst;
  //  exp(computed_dist[id0,] + log(prop_match) - logexptrick_const)/(exp(-logexptrick_const) + sum(exp(computed_dist[id0,] - logexptrick_const + log(prop_match))))
  }
  return(wrap(prob0));
}

//rank_match_1 <- t(sapply(rownames(dist_all), FUN=matchProbs_rank_complete, computed_dist=dist_all, 
//                         prop_match=prop_match_1way))
