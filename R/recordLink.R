#'Probabilistic Patient Record Linkage
#'
#'@param data1 either a binary matrix or dataframe whose rownames are .
#'
#'@param data2 either a binary matrix or a dataframe whose rownames are .
#'
#'@param data1_cont2diff either a matrix or dataframe of continuous features, 
#'such as age, for which the similarity measure uses the difference with 
#'\code{data2_cont2diff}, whose rownames are .
#'
#'@param data2_cont2diff either a matrix or dataframe of continuous features, 
#'such as age, for which the similarity measure uses the difference with 
#'\code{data2_cont1diff}, whose rownames are .
#'
#'@param eps_plus discrepancy rate between \code{data1} and \code{data2}
#'
#'@param eps_minus discrepancy rate between \code{data2} and \code{data1}
#'
#'@param aggreg_2ways a charatcer string indicating how to merge the posterior two 
#'probability matrices obtained for each of the 2 databases. Four possibility are 
#'currently implemented: \code{"maxnorm"}, \code{"max"}, \code{"min"}, \code{"mean"} 
#'and \code{"prod"}. Default is \code{"mean"}.
#'
#'@param min_prev minimum prevalence for the variables used in matching.
#'Default is 1\%.
#'
#'
#'
#'@importFrom landpred VTM
#'@importFrom fGarch dsstd sstdFit
#'
#'@return the posterior probability of matching matrix
#'
#'@export

recordLink <- function(data1, data2, dates1=NULL, dates2=NULL,
                       eps_plus, eps_minus, aggreg_2ways="mean", 
                       min_prev=0.01, 
                       data1_cont2diff, data2_cont2diff,
                       eps_inf, d_max, use_diff=TRUE){
  #browser()
  datesNotNull <- (!is.null(dates1) & !is.null(dates2))
  nb_feat <- ncol(data1)
  if(ncol(data2)!=nb_feat){stop("Number of columns in data2 is different from data1")}
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  if((is.null(data1_cont2diff) || is.null(data2_cont2diff)) && use_diff){
    stop("cannot use_diff when data1_cont2diff and/or data2_cont2diff is NULL")
  }
  if((is.null(dates1) & !is.null(dates2)) | (!is.null(dates1) & is.null(dates2))){
    stop("missing one of the dates tables")
  }
  
  if(use_diff){
    ndiff <- ncol(data1_cont2diff) 
    stopifnot(ncol(data2_cont2diff)==ndiff)
    
    
    if(length(eps_inf)==1 && ndiff>1){
      eps_inf <- rep(eps_inf, ndiff)
    }  
    if(length(d_max)==1 && ndiff>1){
      d_max <- rep(d_max, ndiff)
    }
    stopifnot(length(eps_inf)==ndiff)
    stopifnot(length(d_max)==ndiff)
  }
  
  ind1 <- rownames(data1)
  ind2 <- rownames(data2)
  
  freq_select <- apply(X=(data1==1), MARGIN=2, FUN=mean)>min_prev
  data1_bin <- data1[,freq_select, drop=FALSE]
  data2_bin <- data2[,freq_select, drop=FALSE]
  if(datesNotNull){
    dates1_fs <- dates1[,freq_select, drop=FALSE]
    dates2_fs <- dates2[,freq_select, drop=FALSE]
  }
  #rm(list=c("data1", "data2"))
  
  n_freq10 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.10 | colSums(data1_bin)/nrow(data1_bin)>0.90)
  n_freq5 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.05 | colSums(data1_bin)/nrow(data1_bin)>0.95)
  n_freq1 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.01 | colSums(data1_bin)/nrow(data1_bin)>0.99)
  
  #removing info too rare in any of the 2 database ??
  col_in <- (pmax(apply(X=(data1_bin==0), MARGIN=2, FUN=mean), apply(X=data2_bin==0, MARGIN=2, FUN=mean)) < 0.995 &  
               pmin(apply(data1_bin==0,2,mean),apply(data2_bin==0,2,mean)) > 0.005)
  if(sum(col_in)<1){stop("No features whith adequate prevalence")}
  
  data1_bin <- data1_bin[, col_in, drop=FALSE]
  data2_bin <- data2_bin[, col_in, drop=FALSE]
  if(datesNotNull){
    dates1_fs_ci <- dates1_fs[, col_in, drop=FALSE]
    dates2_fs_ci <- dates2_fs[, col_in, drop=FALSE]
  }
  
  pi1 <- apply(data1_bin==1, MARGIN=2, FUN=mean)
  pi2 <- apply(data2_bin==1, MARGIN=2, FUN=mean)
  K0 <- ncol(data1_bin)
  if(length(eps_plus)==1){
    eps_p <- rep(eps_plus, K0)
  }else if(length(eps_plus)==nb_feat){
    eps_p <- eps_plus[freq_select][col_in]
  }else{
    stop("Length of 'eps_plus' doesn't make sense")
  }
  
  if(length(eps_minus)==1){
    eps_n <- rep(eps_minus, K0)
  }else if(length(eps_minus)==nb_feat){
    eps_n <- eps_minus[freq_select][col_in]
  }else{
    stop("Length of 'eps_minus' doesn't make sense")
  }
  
  
  # need matrices for the C function
  data1_bin <- as.matrix(data1_bin)
  data2_bin <- as.matrix(data2_bin)
  if(datesNotNull){
    dates1_fs_ci <- as.matrix(dates1_fs_ci)
    dates2_fs_ci <- as.matrix(dates2_fs_ci)
  }
  zeros_p <- which(eps_p==0)
  if(length(zeros_p)>0){
    eps_p[zeros_p] <- eps_p[zeros_p] + 1E-6
  }
  zeros_n <- which(eps_n==0)
  if(length(zeros_n)>0){
    eps_n[zeros_n] <- eps_n[zeros_n] + 1E-6
  }
  keep_eps <- which(eps_p!=1 & eps_n!=1)
  #   weights_simil <- cbind("M1-R1"=log((1-eps_n[keep_eps])/pi2[keep_eps]),
  #                          "M0-R0"=log((1-eps_p[keep_eps])/(1-pi2[keep_eps])),
  #                          "M1-R0"=log((eps_n[keep_eps])/(1-pi2[keep_eps])),
  #                          "M0-R1"=log((eps_p[keep_eps])/(pi2[keep_eps])))
  #   weights_simil <- cbind.data.frame("ICD9"=rownames(weights_simil), weights_simil)
  #   write.table(weights_simil, file="weights_simil_AllmedicareXraprod2008.txt", row.names = FALSE, sep="\t")

  if(datesNotNull){
    dist_bin <- loglikC_bin_wDates(Bmat = data2_bin[, keep_eps], Amat = data1_bin[, keep_eps],
                                   Bdates = dates2_fs_ci[, keep_eps], Adates = dates1_fs_ci[, keep_eps],
                                   eps_p = eps_p[keep_eps], eps_n = eps_n[keep_eps], 
                                   piA = pi1[keep_eps], piB = pi2[keep_eps])
  }else{
    dist_bin <- loglikC_bin(Bmat = data2_bin[, keep_eps], Amat = data1_bin[, keep_eps], eps_p = eps_p[keep_eps], 
                             eps_n = eps_n[keep_eps], piA = pi1[keep_eps], piB = pi2[keep_eps])
  }
  
  #browser()
  # colnames(dist_bin) <- rownames(data2_bin)
  # rownames(dist_bin) <- rownames(data1_bin)
  # mean(rowSums(data1[rownames(dist_bin)[-which(rownames(dist_bin) %in% truth$bene_id)][apply(dist_bin[-which(rownames(dist_bin) %in% truth$bene_id), truth$linkage_id], 2, which.max)], ]*data2[truth$linkage_id,]))
  # mean(rowSums(data2[colnames(dist_bin)[-which(colnames(dist_bin) %in% truth$linkage_id)][apply(dist_bin[truth$bene_id, -which(colnames(dist_bin) %in% truth$linkage_id)], 1, which.max)], ]*data1[truth$bene_id,]))
  # mean(rowSums(data2[truth$linkage_id,]*data1[truth$bene_id,]))
  # rm(list=c("data1", "data2"))
  
  if(use_diff){
    pi_inf <- numeric(ndiff)
    for(j in 1:ndiff){
      pi_inf[j] <- mean((unlist(lapply(data1_cont2diff[, j, drop=FALSE], function(x){x-data2_cont2diff[, j, drop=FALSE]}))<d_max[j]))
      #pi_sup <- mean(!(unlist(lapply(data1_cont2diff, function(x){x-data2_cont2diff}))<d_max))
    }
    #     dist_diff <- loglikratioC_diff(Bmat = as.matrix(t(data2_cont2diff)), 
    #                                    Amat =  as.matrix(t(data1_cont2diff)), 
    #                                    d_max = d_max, eps_inf = eps_inf, 
    #                                    pi_inf = pi_inf)
    #arbitrary weight on age
    dist_diff <- loglikratioC_diff_arbitrary(Bmat = as.matrix(t(data2_cont2diff)), 
                                             Amat =  as.matrix(t(data1_cont2diff)), 
                                             d_max = d_max, 
                                             cost = rep(-1E+4, ndiff)
    )
    
    
    dist_all <- dist_bin + dist_diff
  }else{
    dist_all <- dist_bin
  }
  #rm(list=c("data1_bin", "data2_bin"))
  colnames(dist_all) <- ind2
  rownames(dist_all) <- ind1
  #hist(as.vector(dist_all), n=1000)
  
  #   dist_all_2way <- loglikC_mat(Bmat = data1_bin, Amat = data2_bin, eps_p = eps_n*pi1/(1-pi2), 
  #                                eps_n = eps_p*(1-pi1)/pi2, piA = pi2, piB = pi1)
  #   colnames(dist_all_2way) <- row.names(data1_bin)
  #   rownames(dist_all_2way) <- row.names(data2_bin)
  
  
  #loglik_FUN(Bvec = data1_bin[1,], Amat = data2_bin, eps_p = eps_n*pi1/(1-pi2), 
  #           eps_n = eps_p*(1-pi1)/pi2, piA = pi2, piB = pi1)
  
  #sstdFit genrerates warnings - to be ignored
  if(length(dist_all)>10000){
    sstdFit_1way <-  try(fGarch::sstdFit(x=sample(dist_all,10000)))
    if(inherits(sstdFit_1way, "try-error")){sstdFit_1way <-  try(fGarch::sstdFit(x=sample(dist_all,10000)))}
    tmp_est <-  sstdFit_1way$estimate
    for(i in 1:4){
      sstdFit_1way_sub <-  try(fGarch::sstdFit(x=sample(dist_all,10000)))
      if(inherits(sstdFit_1way_sub, "try-error")){sstdFit_1way_sub <-  try(fGarch::sstdFit(x=sample(dist_all,10000)))}
      tmp_est <-  rbind(tmp_est, sstdFit_1way_sub$estimate)
    }
    sstdFit_1way$est <- colMeans(tmp_est, na.rm = TRUE)
  }else{
    sstdFit_1way <-  try(fGarch::sstdFit(x=dist_all))
    if(inherits(sstdFit_1way, "try-error")){sstdFit_1way <-  try(fGarch::sstdFit(x=dist_all))}
    sstdFit_1way$est <-  sstdFit_1way$estimate
  }
  
  #   sstdFit_2way <-  try(sstdFit(x=sample(dist_all_2way,10000)))
  #   if(inherits(sstdFit_2way, "try-error")){sstdFit_2way <-  try(sstdFit(x=sample(dist_all_2way,10000)))}
  #   tmp_est <-  sstdFit_2way$estimate
  #   for(i in 1:4){
  #     sstdFit_2way_sub <-  try(sstdFit(x=sample(dist_all_2way,10000)))
  #     if(inherits(sstdFit_2way_sub, "try-error")){sstdFit_2way_sub <-  try(sstdFit(x=sample(dist_all_2way,10000)))}
  #     tmp_est <-  rbind(tmp_est, sstdFit_2way_sub$estimate)
  #   }
  #   sstdFit_2way$est <-  colMeans(tmp_est)
  
  inter_x <- 0.1
  obs_dist_x_1way <- seq(min(dist_all), max(dist_all), by=inter_x)
  # obs_dist_x_2way <- seq(min(dist_all_2way), max(dist_all_2way), by=inter_x)
  
  ### select cutoff
  # empirical solution
  #single skew-t parametric estimation
  fitskewt_dens_est <- fGarch::dsstd(obs_dist_x_1way,
                                     mean=sstdFit_1way$est[1], 
                                     sd=sstdFit_1way$est[2],
                                     nu=sstdFit_1way$est[3],
                                     xi=sstdFit_1way$est[4])
  peak_fitskewt_dens_est <- max(fitskewt_dens_est)
  rho1_hat_empirical_1way <- min((obs_dist_x_1way[obs_dist_x_1way >= peak_fitskewt_dens_est])[fitskewt_dens_est[obs_dist_x_1way >= peak_fitskewt_dens_est]/peak_fitskewt_dens_est < 1/min(n1,n2)])
  #if(is.infinite(rho1_hat_empirical_1way)){
  #}
  
  #   fitskewt_dens_est <- dsstd(obs_dist_x_2way,
  #                              mean=sstdFit_2way$est[1], 
  #                              sd=sstdFit_2way$est[2],
  #                              nu=sstdFit_2way$est[3],
  #                              xi=sstdFit_2way$est[4])
  #   peak_fitskewt_dens_est <- max(fitskewt_dens_est)
  #   rho1_hat_empirical_2way <- min((obs_dist_x_2way[obs_dist_x_2way >= peak_fitskewt_dens_est])[fitskewt_dens_est[obs_dist_x_2way >= peak_fitskewt_dens_est]/peak_fitskewt_dens_est < 1/min(n1,n2)])
  #   
  
  # analytical solution
  nu <-  sstdFit_1way$est[3]
  s <-  sstdFit_1way$est[2]
  m <-  sstdFit_1way$est[1]
  xi <-  sstdFit_1way$est[4]
  Beta <- exp(lgamma(0.5) - lgamma((nu+1)/2) + lgamma(nu/2))
  m1 <- 2*sqrt(nu-2)/((nu-1)*Beta)
  mu <- m1*(xi-1/xi)
  sigma <- sqrt((1-m1^2)*(xi^2+1/xi^2)+2*m1^2-1)
  d1max_x <- m-(mu-xi*sqrt((nu-2)/(nu+2)))*s/sigma  
  #m-mu*s/sigma
  #m-(mu+1/xi*sqrt((nu-2)/(nu+2)))*s/sigma
  
  #   nu <-  sstdFit_2way$est[3]
  #   s <-  sstdFit_2way$est[2]
  #   m <-  sstdFit_2way$est[1]
  #   xi <-  sstdFit_2way$est[4]
  #   Beta <- gamma(0.5)/gamma((nu+1)/2)*gamma(nu/2)
  #   m1 <- 2*sqrt(nu-2)/((nu-1)*Beta)
  #   mu <- m1*(xi-1/xi)
  #   sigma <- sqrt((1-m1^2)*(xi^2+1/xi^2)+2*m1^2-1)
  #   rho1_hat_2way <- m-(mu-xi*sqrt((nu-2)/(nu+2)))*s/sigma
  
  
  d1 <- function(x){
    #c <- 2*sigma*gamma((nu+1)/2)/(s*(xi+1/xi)*sqrt(pi*nu-2)*gamma(nu/2))
    c <- exp(log(2) + log(sigma) + lgamma((nu+1)/2) - log(s*(xi+1/xi)*sqrt(pi*nu-2)) - lgamma(nu/2))
    d <- -c*(nu+1)*sigma/(xi^2*(nu-2)*s)
    y <- (x-m)/s*sigma+mu
    return(d*y*(1+y^2/(xi^2*(nu-2)))^(-(nu+3)/2))
  }
  
  d2 <- function(x){
    #c <- 2*sigma*gamma((nu+1)/2)/(s*(xi+1/xi)*sqrt(pi*nu-2)*gamma(nu/2))
    c <- exp(log(2) + log(sigma) + lgamma((nu+1)/2) - log(s*(xi+1/xi)*sqrt(pi*nu-2)) - lgamma(nu/2))
    d <- -c*(nu+1)*sigma/(xi^2*(nu-2)*s)
    y <- (x-m)/s*sigma+mu
    return(d*sigma/s*(1+y^2/(xi^2*(nu-2)))^(-(nu+5)/2)*(1-y^2*(nu+2)/xi^2*(nu-2)))
  }
  
  rho1_hat <- obs_dist_x_1way[obs_dist_x_1way>d1max_x][min(which(abs(d1(obs_dist_x_1way[obs_dist_x_1way>d1max_x]))<(1/min(n1,n2)) & abs(d2(obs_dist_x_1way[obs_dist_x_1way>d1max_x]))<(1/min(n1,n2))))]
  if(is.na(rho1_hat)){
    rho1_hat <- obs_dist_x_1way[obs_dist_x_1way>d1max_x][min(which(abs(d1(obs_dist_x_1way[obs_dist_x_1way>d1max_x]))<(1/min(n1,n2))))]
  }
  if(is.na(rho1_hat)){
    rho1_hat <- 0
  }
  rm(list="obs_dist_x_1way")
  
  
  dist_all_select <- dist_all > rho1_hat #rho1_hat_empirical_1way
  # sum(dist_all > d1max_x)
  # dist_all_2way_select <- dist_all_2way > rho1_hat_2way
  nmatch_hat_1way <- min(sum(colSums(dist_all_select)>=1), sum(rowSums(dist_all_select)>=1)) # Empirical Bayes prior on P(M=1)
  rm(list="dist_all_select")
  
  prop_match_1way <- nmatch_hat_1way/(n1*n2)
  # nmatch_hat_2way <- sum(dist_all_2way_select) # Empirical Bayes prior on P(M=1)
  # prop_match_2way <- nmatch_hat_2way/(nrow(data1_bin)*nrow(data2_bin))
  
  rank_match_1 <- matchProbs_rank_full_C(dist_all+100, prop_match=prop_match_1way)
  rownames(rank_match_1) <- ind1
  colnames(rank_match_1) <- ind2
  rank_match_2 <- t(matchProbs_rank_full_C(t(dist_all), prop_match=prop_match_1way))
  rm(list="dist_all")
  rownames(rank_match_2) <- ind1
  colnames(rank_match_2) <- ind2
  mode(rank_match_1) <- "numeric" 
  mode(rank_match_2) <- "numeric"
  
  if(aggreg_2ways == "maxnorm"){
    rank_match_add0_1 <- cbind(rank_match_1, 1-rowSums(rank_match_1))
    probMax_1 <- apply(rank_match_add0_1, MARGIN=1, FUN=max)
    rm(list="rank_match_add0_1")
    rank_match_add0_2 <- rbind(rank_match_2, 1-colSums(rank_match_2))
    probMax_2 <- apply(rank_match_add0_2, MARGIN=2, FUN=max)
    rm(list="rank_match_add0_2")
  }
  
  probMatch_mat <- switch(aggreg_2ways,
                          maxnorm = rank_match_1*rank_match_2/matrix(probMax_1, ncol=1)%*%matrix(probMax_2, nrow=1),
                          min = pmin(rank_match_1, rank_match_2),
                          max = pmax(rank_match_1, rank_match_2),
                          mean = (rank_match_1+rank_match_2)/2,
                          prod = rank_match_1*rank_match_2
  )
  
  
  return(probMatch_mat)
  
}


