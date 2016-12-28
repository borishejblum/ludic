#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]

//'C++ implementation of Winkler's Method E step
//'using a sparse agreement  matrix
//'
//'@keywords internal
// [[Rcpp::export]]
arma::mat estep_C_sparse(arma::sp_mat agreemat, double p, arma::rowvec m, arma::rowvec u, bool Log=false){
  
  int K = agreemat.n_cols;  
  int N = agreemat.n_rows;
  
  mat res(N, 2);
  
  for(int i=0; i<N; i++){
    rowvec agreerow_i(agreemat.row(i));
    double a_log = log(p) + sum(agreerow_i%log(m) + (1.0 - agreerow_i)%log(1.0 - m));
    double b_log = log(1.0 - p) + sum(agreerow_i%log(u) + (1.0 - agreerow_i)%log(1.0 - u));
    double max_log = max(vec(a_log, b_log));
    res(i, 0) = a_log - (max_log + log(exp(a_log - max_log) + exp(b_log - max_log)));
    res(i, 1) = b_log - (max_log + log(exp(a_log - max_log) + exp(b_log - max_log)));
  }
  
  if(!Log){
    res = exp(res);
  }
  
  return(res);
}


//'Fast C++ implementation of Winkler's Method E step
//'using a sparse agreement  matrix
//'
//'@keywords internal
// [[Rcpp::export]]
arma::mat estep_C_vect_sparse(arma::sp_mat agreemat, double p, arma::colvec m, arma::colvec u, bool Log=false){
  
  int K = agreemat.n_cols;  
  int N = agreemat.n_rows;
  
  mat res(N, 2);
  
  sp_mat agreemat_neg = sp_mat(1.0 - mat(agreemat));
  vec a_log = log(p) + agreemat*log(m) + agreemat_neg*log(1.0 - m);
  vec b_log = log(1.0 - p) + agreemat*log(u) + agreemat_neg*log(1.0 - u);
  res.col(0) = - log(1.0 + exp(b_log - a_log));//exp(a_log - (a_log + log(1.0 + exp(b_log - a_log))));
  res.col(1) = b_log - (a_log + log(1.0 + exp(b_log - a_log)));
  
  if(!Log){
    res = exp(res);
  }
  
  return(res);
}

//'Fast C++ implementation of Winkler's Method E step
//'using two sparse matrices 
//'
//'Saving memory while conserving speed by avoiding casting agreemat_neg
//'
//'@keywords internal
// [[Rcpp::export]]
arma::mat estep_C_vect_sparse2(arma::sp_mat agreemat, arma::sp_mat agreemat_neg, double p, arma::colvec m, arma::colvec u, bool Log=false){
  
  int K = agreemat.n_cols;  
  int N = agreemat.n_rows;
  
  mat res(N, 2);
  
  vec a_log = log(p) + agreemat*log(m) + agreemat_neg*log(1.0 - m);
  vec b_log = log(1.0 - p) + agreemat*log(u) + agreemat_neg*log(1.0 - u);
  res.col(0) = - log(1.0 + exp(b_log - a_log));//exp(a_log - (a_log + log(1.0 + exp(b_log - a_log))));
  res.col(1) = b_log - (a_log + log(1.0 + exp(b_log - a_log)));
  
  if(!Log){
    res = exp(res);
  }
  
  return(res);
}


//'C++ implementation of Winkler's Method E step
//'
//'@keywords internal
// [[Rcpp::export]]
arma::mat estep_C(arma::mat agreemat, double p, arma::rowvec m, arma::rowvec u, bool Log=false){
  
  int K = agreemat.n_cols;  
  int N = agreemat.n_rows;
  
  mat res(N, 2);
  
  for(int i=0; i<N; i++){
    rowvec agreerow_i(agreemat.row(i));
    double a_log = log(p) + sum(agreerow_i%log(m) + (1.0 - agreerow_i)%log(1.0 - m));
    double b_log = log(1.0 - p) + sum(agreerow_i%log(u) + (1.0 - agreerow_i)%log(1.0 - u));
    double max_log = max(NumericVector::create(a_log, b_log));
    res(i, 0) = a_log - (max_log + log(exp(a_log - max_log) + exp(b_log - max_log)));
    res(i, 1) = b_log - (max_log + log(exp(a_log - max_log) + exp(b_log - max_log)));
  }
  
  if(!Log){
    res = exp(res);
  }
  
  return(res);
}

//'Fast C++ implementation of Winkler's Method E step
//'
//'@keywords internal
// [[Rcpp::export]]
arma::mat estep_C_vect(arma::mat agreemat, double p, arma::colvec m, arma::colvec u, bool Log=false){
  
  int K = agreemat.n_cols;  
  int N = agreemat.n_rows;
  
  mat res(N, 2);
  
  mat agreemat_neg = 1.0 - agreemat;
  vec a_log = log(p) + agreemat*log(m) + agreemat_neg*log(1.0 - m);
  vec b_log = log(1.0 - p) + agreemat*log(u) + agreemat_neg*log(1.0 - u);
  res.col(0) = - log(1.0 + exp(b_log - a_log));//exp(a_log - (a_log + log(1.0 + exp(b_log - a_log))));
  res.col(1) = b_log - (a_log + log(1.0 + exp(b_log - a_log)));
  
  if(!Log){
    res = exp(res);
  }
  
  return(res);
}



//'C++ implementation of the E and M steps from Winkler's EM algorithm estimating FS method 
//'using sparse matrices for big sample sizes
//'
//'@keywords internal
// [[Rcpp::export]]

List EMstep_C_sparse_big(arma::mat mat_A, arma::mat mat_B, double p, arma::rowvec m, arma::rowvec u){
  
  int K = mat_A.n_cols;  
  int nA = mat_A.n_rows;
  int nB = mat_B.n_rows;
  int N = nA*nB;
  
  vec g_m_log(N);
  vec g_u_log(N);
  
  //E step
  for(int j=0; j<nB; j++){
    for(int i=0; i<nA; i++){
      int l = j*nA + i;
      rowvec agreerow_pairl_neg(abs(mat_A.row(i) - mat_B.row(j)));
      rowvec agreerow_pairl(1.0 - agreerow_pairl_neg);
      double a_log = log(p) + sum(agreerow_pairl%log(m) + (1.0 - agreerow_pairl_neg)%log(1.0 - m));
      double b_log = log(1.0 - p) + sum(agreerow_pairl%log(u) + (1.0 - agreerow_pairl_neg)%log(1.0 - u));
      double max_log = max(NumericVector::create(a_log, b_log));
      g_m_log(l) = a_log - (max_log + log(exp(a_log - max_log) + exp(b_log - max_log)));
      g_u_log(l) = b_log - (max_log + log(exp(a_log - max_log) + exp(b_log - max_log)));
    }
  }
  
  vec m_res(zeros(K));
  vec u_res(zeros(K));
  
  //M step
  double gm_max = max(g_m_log);
  double gu_max = max(g_u_log);
  double gm_summed_log = gm_max + log(exp(sum(g_m_log - gm_max)));
  double gu_summed_log = gu_max + log(exp(sum(g_u_log - gu_max)));
  double p_res = exp(gm_summed_log - log(double(N)));
  for(int j=0; j<nB; j++){
    for(int i=0; i<nA; i++){
      int l = j*nA + i;
      colvec agreerow_pairl(conv_to<colvec>::from(1.0 - abs(mat_A.row(i) - mat_B.row(j))));
      m_res += exp(g_m_log(l) - gm_summed_log)*agreerow_pairl;
      u_res += exp(g_u_log(l) - gu_summed_log)*agreerow_pairl;
    }
  }
  return(Rcpp::List::create(Rcpp::Named("p") = p_res, Rcpp::Named("m") = m_res.t(), 
                            Rcpp::Named("u") = u_res.t()));
}

//'C++ implementation of the E and M steps from Winkler's EM algorithm estimating FS method 
//'using sparse matrices for big sample sizes
//'
//'@keywords internal
//'@export
// [[Rcpp::export]]
arma::mat Estep_C_sparse_big(arma::mat mat_A, arma::mat mat_B, double p, arma::rowvec m, arma::rowvec u, bool Log=false){
  
  int K = mat_A.n_cols;  
  int nA = mat_A.n_rows;
  int nB = mat_B.n_rows;
  int N = nA*nB;
  
  mat g_log(N, 2);
  
  //E step
  for(int j=0; j<nB; j++){
    for(int i=0; i<nA; i++){
      int l = j*nA + i;
      rowvec agreerow_pairl_neg(abs(mat_A.row(i) - mat_B.row(j)));
      rowvec agreerow_pairl(1.0 - agreerow_pairl_neg);
      double a_log = log(p) + sum(agreerow_pairl%log(m) + (1.0 - agreerow_pairl_neg)%log(1.0 - m));
      double b_log = log(1.0 - p) + sum(agreerow_pairl%log(u) + (1.0 - agreerow_pairl_neg)%log(1.0 - u));
      double max_log = max(NumericVector::create(a_log, b_log));
      g_log(l, 0) = a_log - (max_log + log(exp(a_log - max_log) + exp(b_log - max_log)));
      g_log(l, 1) = b_log - (max_log + log(exp(a_log - max_log) + exp(b_log - max_log)));
    }
  }
  
  if(Log){
    return(g_log);
  }else{
    return(exp(g_log));
  }
}



//'C++ implementation of the E and M steps from Winkler's EM algorithm estimating FS method 
//'using sparse matrices for big sample sizes
//'
//'@keywords internal
//'@export
// [[Rcpp::export]]
List Mstep_C_sparse_big(arma::mat mat_A, arma::mat mat_B, arma::vec g_m, arma::vec g_u, bool Log=false){
  
  int K = mat_A.n_cols;  
  int nA = mat_A.n_rows;
  int nB = mat_B.n_rows;
  int N = nA*nB;
  
  vec m_res(zeros(K));
  vec u_res(zeros(K));
  double p_res;
  
  //M step
  if(Log){
    double gm_max = max(g_m);
    double gu_max = max(g_u);
    double gm_summed_log = gm_max + log(exp(sum(g_m - gm_max)));
    double gu_summed_log = gu_max + log(exp(sum(g_u - gu_max)));
    p_res = exp(gm_summed_log - log(double(N)));
    for(int j=0; j<nB; j++){
      for(int i=0; i<nA; i++){
        int l = j*nA + i;
        colvec agreerow_pairl(conv_to<colvec>::from(1.0 - abs(mat_A.row(i) - mat_B.row(j))));
        m_res += exp(g_m(l) - gm_summed_log)*agreerow_pairl;
        u_res += exp(g_u(l) - gu_summed_log)*agreerow_pairl;
      }
    }
  }else{
    double gm_summed = sum(g_m);
    double gu_summed = sum(g_u);
    p_res = gm_summed/N;
    for(int j=0; j<nB; j++){
      for(int i=0; i<nA; i++){
        int l = j*nA + i;
        colvec agreerow_pairl(conv_to<colvec>::from(1.0 - abs(mat_A.row(i) - mat_B.row(j))));
        m_res += g_m(l)*agreerow_pairl / gm_summed;
        u_res += g_u(l)*agreerow_pairl / gu_summed;
      }
    }
  }
  return(Rcpp::List::create(Rcpp::Named("p") = p_res, Rcpp::Named("m") = m_res.t(), 
                            Rcpp::Named("u") = u_res.t()));
}

