#'Kernel for matching binary datasets by Probabilistic Patient Record Linkage
#'
#'@param data a matrix of original (continuous) dataset, either a matrix or a dataframe.
#'
#'
#'@param MAF0 Minor Allele Frequency. Default is 0.2
#'
#'@param OR Odds Ratio. Default is 1.3.
#'
#'@param k0 the number of top matching probabilities to store. Default is 5 top probabilities.
#'
#'@param eps1 probability that a feature is recorded in the first database but 
#'not in the second (mis-match or recording error). Default is 0.01.
#'
#'@param eps2 probability that a feature is recorded in the second database but 
#'not in the first (mis-match or recording error). Default is the same as \code{eps1}.
#'
#'@param aggreg_2ways a charatcer string indicating how to merge the posterior two 
#'probability matrices obtained for each of the 2 databases. Four possibility are 
#'currently implemented: \code{"maxnorm"}, \code{"max"}, \code{"min"}, \code{"mean"} 
#'and \code{"prod"}. Default is \code{"maxnorm"}.
#'
#'@param nb_perturb the number of perturbation used for the pvalue combination.
#'Default is 1000.
#'
#'@param outputfile_name a character string for the name of the .txt file where the 
#'results should be stored. Must be provided without extension. Default is 'Default_output'.
#'
#'@importFrom landpred VTM
#'@importFrom fGarch dsstd sstdFit
#'@importFrom mvnfast rmvn
#'
#'@return append the results in the file '\code{outputfile_name}.txt'.
#'
#'@export

match_2in1_2datasets_inter <- function(data1, data2, ind_match, 
                                       MAF0 = 0.2, OR = 1.3, k0 = 5,
                                       eps1 = 0.01, eps2 = eps1,
                                       aggreg_2ways = c("maxnorm", "max", "min", "mean", "prod"), 
                                       nb_perturb = 1000,
                                       outputfile_name="Default_output"){
  #browser()
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  nb_match <- length(ind_match)
  
  ind1 <- rownames(data1)
  ind2 <- rownames(data2)
  id_match <- as.character(ind_match)
  yes_match0 <- (ind1 == landpred::VTM(ind2,length(ind1)))
  
  
  freq_select <- apply(X=(data1==1), MARGIN=2, FUN=mean)>0.01
  data1_bin <- data1[,freq_select]
  data2_bin <- data2[,freq_select]
  
  tmp_Genom <- rbinom(n1, size=2, prob=MAF0)
  tmp_Pheno <- rbinom(n1, size=1, prob=expit(-0.5+log(OR)*tmp_Genom))
  data1_all <- cbind(data1[,freq_select],"Pheno"=tmp_Pheno,"Genom"=tmp_Genom)
  data1_out <- data1_all[,c("Pheno","Genom")]
  
  tmp_Genom <- c(data1_out[id_match, "Genom"], 
                 rbinom(n2-nb_match, size=2, prob=MAF0))
  tmp_Pheno <- c(data1_out[id_match, "Pheno"], 
                 rbinom(n2-nb_match, size=1, prob=expit(-0.5+log(OR)*tmp_Genom[(nb_match+1):n2])))
  data2_all <- cbind(data2[,freq_select],"Pheno"=tmp_Pheno,"Genom"=tmp_Genom)
  data2_out <- data2_all[,c("Pheno","Genom")]
  
  
  
  epsilon_kminus_truth <- sapply(1:ncol(data1_bin),function(v){
    mean(data2_bin[id_match[data1_bin[id_match,v]==1],v]==0)
  })
  epsilon_kplus_truth <- sapply(1:ncol(data2_bin),function(v){
    mean(data1_bin[id_match[data2_bin[id_match,v]==1],v]==0)
  })
  
  epsilon_truth <- c(mean(epsilon_kminus_truth, na.rm = TRUE), mean(epsilon_kplus_truth, na.rm = TRUE),
                     mean(abs(epsilon_kminus_truth -mean(epsilon_kminus_truth, na.rm = TRUE)), na.rm = TRUE),
                     mean(abs(epsilon_kplus_truth-mean(epsilon_kplus_truth, na.rm = TRUE)), na.rm = TRUE)
  )
  
  
  n_freq10 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.10 | colSums(data1_bin)/nrow(data1_bin)>0.90)
  n_freq5 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.05 | colSums(data1_bin)/nrow(data1_bin)>0.95)
  n_freq1 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.01 | colSums(data1_bin)/nrow(data1_bin)>0.99)
  
  
  
  #removing info too rare in any of the 2 database ??
  col_in <- (pmax(apply(X=(data1_bin==0), MARGIN=2, FUN=mean), apply(X=data2_bin==0, MARGIN=2, FUN=mean)) < 0.995 &  
               pmin(apply(data1_bin==0,2,mean),apply(data2_bin==0,2,mean)) > 0.005)
  data1_bin <- data1_bin[,col_in]
  data2_bin <- data2_bin[,col_in]
  
  pi1 <- apply(data1_bin==1, MARGIN=2, FUN=mean)
  pi2 <- apply(data2_bin==1, MARGIN=2, FUN=mean)
  K0 <- ncol(data1_bin)
  eps_p <- rep(eps1, K0)
  eps_n <- rep(eps2, K0)
  
  
  # need matrices for the C function
  data1_bin <- as.matrix(data1_bin)
  data2_bin <- as.matrix(data2_bin)
  
  dist_all <- loglikC_mat_fast(Bmat = data2_bin, Amat = data1_bin, eps_p = eps_p, 
                                 eps_n = eps_n, piA = pi1, piB = pi2)
  
  
  colnames(dist_all) <- rownames(data2_bin)
  rownames(dist_all) <- rownames(data1_bin)
  
  #   dist_all_2way <- loglikC_mat(Bmat = data1_bin, Amat = data2_bin, eps_p = eps_n*pi1/(1-pi2), 
  #                                eps_n = eps_p*(1-pi1)/pi2, piA = pi2, piB = pi1)
  #   colnames(dist_all_2way) <- row.names(data1_bin)
  #   rownames(dist_all_2way) <- row.names(data2_bin)
  
  
  #loglik_FUN(Bvec = data1_bin[1,], Amat = data2_bin, eps_p = eps_n*pi1/(1-pi2), 
  #           eps_n = eps_p*(1-pi1)/pi2, piA = pi2, piB = pi1)
  
  #sstdFit genrerates warnings - to be ignored
  sstdFit_1way <-  try(fGarch::sstdFit(x=sample(dist_all,10000)))
  if(inherits(sstdFit_1way, "try-error")){sstdFit_1way <-  try(fGarch::sstdFit(x=sample(dist_all,10000)))}
  tmp_est <-  sstdFit_1way$estimate
  for(i in 1:4){
    sstdFit_1way_sub <-  try(fGarch::sstdFit(x=sample(dist_all,10000)))
    if(inherits(sstdFit_1way_sub, "try-error")){sstdFit_1way_sub <-  try(fGarch::sstdFit(x=sample(dist_all,10000)))}
    tmp_est <-  rbind(tmp_est, sstdFit_1way_sub$estimate)
  }
  sstdFit_1way$est <- colMeans(tmp_est)
  
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
  
  dist_all_select <- dist_all > rho1_hat #rho1_hat_empirical_1way
  # sum(dist_all > d1max_x)
  # dist_all_2way_select <- dist_all_2way > rho1_hat_2way
  nmatch_hat_1way <- min(sum(colSums(dist_all_select)>=1), sum(rowSums(dist_all_select)>=1)) # Empirical Bayes prior on P(M=1)
  prop_match_1way <- nmatch_hat_1way/(nrow(data1_bin)*nrow(data2_bin))
  # nmatch_hat_2way <- sum(dist_all_2way_select) # Empirical Bayes prior on P(M=1)
  # prop_match_2way <- nmatch_hat_2way/(nrow(data1_bin)*nrow(data2_bin))
  
  names_p_ori <- paste0("p_ori",1:k0)   
  names_p_agg <- paste0("p_agg",1:k0)
  names_ID <- paste0("ID", 1:k0)
  
  
  
  betahat_all <- NULL 
  accuracy_all <- NULL 
  nmatch_final_all <- NULL

  rank_match_1 <- matchProbs_rank_full_C(dist_all, prop_match=prop_match_1way)
  rownames(rank_match_1) <- rownames(dist_all)
  colnames(rank_match_1) <- colnames(dist_all)
  
  rank_match_2 <- t(matchProbs_rank_full_C(t(dist_all), prop_match=prop_match_1way))
  rownames(rank_match_2) <- rownames(dist_all)
  colnames(rank_match_2) <- colnames(dist_all)
  
  mode(rank_match_1) <- "numeric" 
  mode(rank_match_2) <- "numeric"
  
  rank_match_add0_1 <- cbind(rank_match_1, 1-rowSums(rank_match_1))
  rank_match_add0_2 <- rbind(rank_match_2, 1-colSums(rank_match_2))
  
  rankMatch_global <- matrix(NA, ncol=ncol(rank_match_1), nrow=nrow(rank_match_1))
  probMax_1 <- apply(rank_match_add0_1, MARGIN=1, FUN=max)
  probMax_2 <- apply(rank_match_add0_2, MARGIN=2, FUN=max)
  
  probMax_mat <- switch(aggreg_2ways,
                        maxnorm = rank_match_1*rank_match_2/matrix(probMax_1, ncol=1)%*%matrix(probMax_2, nrow=1),
                        min = pmin(rank_match_1, rank_match_2),
                        max = pmax(rank_match_1, rank_match_2),
                        mean = (rank_match_1+rank_match_2)/2,
                        prod = rank_match_1*rank_match_2
  )
  
  
  #plot(density(diag(probMax_mat[id_match, id_match])))
  
  names_temp <- "p_ori1"
  names_prob <- paste0("p_ori", 1:k0)
  
  temp_id <- sapply(rownames(probMax_mat), function(i){names(sort(probMax_mat[i,], decreasing=TRUE))[1:k0]})
  rownames(temp_id) <- paste0("ID", 1:k0)
  mode(temp_id) <- "character"
  
  rank_match_probs <- sapply(rownames(probMax_mat), function(i){sort(probMax_mat[i,], decreasing=TRUE)[1:k0]})
  
  mydata_match <- as.matrix(data1_out[rownames(probMax_mat),])
  mydata_match <- cbind(mydata_match, t(matrix(data2_out[c(temp_id),"Genom"],nrow=k0)), t(rank_match_probs))
  mode(mydata_match) <- "numeric" 
  colnames(mydata_match)[2+1:k0] = paste0("G",1:k0)
  colnames(mydata_match)[2+k0 +1:k0] <- names_prob
  mydata_match <- cbind(mydata_match,"id"=t(temp_id))
  mode(mydata_match) = "numeric"
  
  mydata0 <- mydata_match[id_match,]
  oracle_summary <- summary(glm(mydata0[,"Pheno"]~mydata0[,"Genom"],family=binomial))$coef[2,-3]
  
  thresholds <- seq(from=0.1, to=0.9, by=0.2)#c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95)#
  nb_thresholds <- length(thresholds)
  
  for(cut_p in thresholds){
    if(aggreg_2ways=="maxnorm"){
      matchClaimed_mat <- ((probMax_mat>cut_p) & (rank_match_1>cut_p) & (rank_match_2>cut_p))
    }else{
      matchClaimed_mat <- (probMax_mat>cut_p)
    }
    matchClaimed_mat_indexes <- which(matchClaimed_mat, arr.ind = TRUE)
    TPR <- sum(rownames(probMax_mat)[matchClaimed_mat_indexes[,1]] == colnames(probMax_mat)[matchClaimed_mat_indexes[,2]])/length(id_match)
    FPR <- sum(rownames(probMax_mat)[matchClaimed_mat_indexes[,1]] != colnames(probMax_mat)[matchClaimed_mat_indexes[,2]])/(nrow(data1_bin)*nrow(data2_bin)-length(id_match))
    PPV <- sum(rownames(probMax_mat)[matchClaimed_mat_indexes[,1]] == colnames(probMax_mat)[matchClaimed_mat_indexes[,2]])/nrow(matchClaimed_mat_indexes)
    NPV <- (nrow(data1_bin)*nrow(data2_bin)-nrow(matchClaimed_mat_indexes)-(length(id_match)-sum(rownames(probMax_mat)[matchClaimed_mat_indexes[,1]] == colnames(probMax_mat)[matchClaimed_mat_indexes[,2]])))/(nrow(data1_bin)*nrow(data2_bin)-nrow(matchClaimed_mat_indexes))
    n_match_claimed <- sum(matchClaimed_mat)
    #cat(cut_p, "\n-----\n TPR:", TPR, "; FPR:", FPR, "; PPV:", PPV, "; NPV:", NPV, "; N_match:", n_match_claimed, "\n)
    
    mydata1 <- mydata_match[mydata_match[,names_prob[1]] > cut_p,]
    #setdiff(which(rownames(mydata_match)%in%id_match), which(mydata_match[,names_prob[1]] > cut_p))
    #setdiff(which(mydata_match[,names_prob[1]] > cut_p), which(rownames(mydata_match)%in%id_match))
    mydata1 <- cbind(mydata1, "G.impute" = apply(mydata1[,paste0("G",1:k0)]*mydata1[,names_prob],1,sum)/
                       apply(mydata1[,names_prob],1,sum))
    temp_1 <- try(summary(glm(mydata1[,"Pheno"]~mydata1[,"G1"],family=binomial))$coef[2,-3], silent = TRUE)
    if(inherits(temp_1, "try-error")){temp_1 <- rep(NA, 3)}
    temp_ave <- try(summary(glm(mydata1[,"Pheno"]~mydata1[,"G.impute"],family=binomial))$coef[2,-3], silent = TRUE)
    if(inherits(temp_ave, "try-error")){temp_ave <- rep(NA, 3)}
    beta.hat <- rbind("Oracle"=oracle_summary,
                      "Best.1"=temp_1,
                      "Best.ave"=temp_ave)
    
    accuracy_hat <- c("TPR"=TPR, "FPR"=FPR, "PPV"=PPV, "NPV"=NPV)
    
    nmatch.hat_final <-  sum(matchClaimed_mat)
    
    betahat_all <- rbind(betahat_all,beta.hat)
    accuracy_all <-  rbind(accuracy_all,accuracy_hat)
    nmatch_final_all <- rbind(nmatch_final_all, nmatch.hat_final)
  }
  
  
  
  
  
  selected_thres <- 1:length(thresholds)#c(3,9)#c(0:4)*2+1##
  
  ## OMNIBUS estimator: optimally combine thresholds with 'otpimal' weights#### 
  ## PART 1 ----
  Z <- cbind("Intercept"=rep(1, nrow(mydata_match)), 
             "G.imput"=(apply(mydata_match[,paste0("G",1:k0)]*mydata_match[,names_prob],1,sum)/
                          apply(mydata_match[,names_prob],1,sum)))
  
  eta <- list()
  theta_thresholds <- NULL
  for(cut_p in thresholds){
    xi <- mydata_match[,names_prob[1]] > cut_p
    n_rho <- sum(xi)
    mydata1 <- mydata_match[xi,]
    mydata1 <- as.data.frame(cbind(mydata1, "G.impute" = apply(mydata1[,paste0("G",1:k0)]*mydata1[,names_prob],1,sum)/
                                     apply(mydata1[,names_prob],1,sum)))
    fit <- glm(Pheno~G.impute, data=mydata1, family=binomial)
    #inv_fishInfo_mat <- summary(fit)$cov.scaled
    theta <- summary(fit)$coef[,"Estimate", drop=FALSE]
    theta_thresholds <-c(theta_thresholds, theta[2])
    I_rho <-  1/n_rho*matrix(apply(Z, MARGIN=1, function(r){tcrossprod(r)*(expit_dev1(r%*%theta)[1,])})%*%xi, ncol=2)
    eta[[as.character(cut_p)]] <- 1/n_rho*rbind((1*xi), (1*xi))*solve(I_rho)%*%t(Z)%*%diag(x=(mydata_match[,"Pheno"] - expit(Z%*%theta)[,1]))
    sqrt(crossprod(sapply(X=eta, FUN=function(m){m[2,]})))
  }
  eta_mat <- sapply(X=eta, FUN=function(m){m[2,]})
  sigma <- crossprod(eta_mat[, selected_thres])
  sds <- sqrt(diag(sigma))
  
  #Combine p-values
  pvals <- pval_zscore(theta_thresholds, sds)
  fisher_comb <- comb_pvals(pvals)
  
  B <- nb_perturb #number of perturbations
  theta_thresholds_star <- crossprod(eta_mat, matrix(rnorm(n1*B, mean = 0, sd = 1), nrow=n1, ncol=B))
  pvals_star <- apply(theta_thresholds_star, MARGIN=2, FUN=pval_zscore, sigma=sds)
  fisher_comb_star <- apply(pvals_star, MARGIN=2, FUN=comb_pvals)
  pval_combined <- 1-sum(fisher_comb>=fisher_comb_star)/B
  
  
  #### Predict estimate at threshold 1 ----
  beta_thresholds_weighted <- betahat_all[seq(3, nrow(betahat_all), by=3),]
  rownames(beta_thresholds_weighted) <- thresholds
  lmfit_threshold_data <- cbind.data.frame("y"=beta_thresholds_weighted[selected_thres,"Estimate"],
                                           "ylog"=log(beta_thresholds_weighted[selected_thres,"Estimate"]),
                                           "x"=thresholds[selected_thres],
                                           "y.se" = beta_thresholds_weighted[selected_thres,"Std. Error"])
  lmfit_threshold <- lm(y~x, data=lmfit_threshold_data)#, weights=(1/lmfit_threshold_data$x.se^2))
  #lmlogfit_threshold <- lm(ylog~x, data=lmfit_threshold_data)
  if(lmfit_threshold$coefficients["x"]<=0){
    est_threshold_pred <- lmfit_threshold$coefficients["(Intercept)"]
  }else{
    est_threshold_pred <- predict.lm(lmfit_threshold, se.fit = T, newdata=cbind.data.frame("x"=1))$fit 
  }
  #est_threshold_pred$se.fit
  # all.equal(sqrt(c(1,1)%*%vcov(lmfit_threshold)%*%t(t(c(1,1)))), est_threshold_pred$se.fit)
  # ss <- sqrt(sum((lmfit_threshold_data$y-mean(lmfit_threshold_data$y))^2)/(length(lmfit_threshold_data$y)-1))
  # ss*sqrt(1/length(lmfit_threshold_data$y) + ((1-mean(lmfit_threshold_data$x))^2)/sum((lmfit_threshold_data$x-mean(lmfit_threshold_data$x))^2))   
  
  
  
  ##### Parametric bootstrap for predicted threshold 1 ----
  B <- 5000
  est_threshold_pred_ptb <- numeric(B)
  for(i in 1:B){
    #     y.pb<-NULL
    #     for(j in 1:nrow(lmfit_threshold_data)){
    #       y.pb <- c(y.pb, rnorm(1, mean=lmfit_threshold_data$y[j], sd=lmfit_threshold_data$y.se[j]))
    #     }
    #     data_ptb <- cbind.data.frame(y.pb, "x"=lmfit_threshold_data$x)
    data_ptb <- cbind.data.frame("y.pb"=mvnfast::rmvn(1, mu=lmfit_threshold_data$y, sigma=sigma)[1,],#rnorm(n=nrow(lmfit_threshold_data), mean=lmfit_threshold_data$y, sd=lmfit_threshold_data$y.se),
                                 "x"=lmfit_threshold_data$x
    )
    lmfit_threshold_ptb <- lm(y.pb~x, data=data_ptb)
    if(lmfit_threshold_ptb$coefficients["x"]<=0){
      est_threshold_pred_ptb[i] <- lmfit_threshold_ptb$coefficients["(Intercept)"]
    }else{
      est_threshold_pred_ptb[i] <- predict.lm(lmfit_threshold_ptb, se.fit = T, newdata=cbind.data.frame("x"=1))$fit 
    }
  }
  
  
  
  
  ## OMNIBUS estimator: optimally combine thresholds with influence functions derived weights#### 
  ## PART 2 ----
  one <- rep(1, nrow(sigma))
  
  #   sigma_inv <- solve(sigma)
  #   weights <- (one%*%sigma_inv/(one%*%sigma_inv%*%one)[1,1])[1,]
  #   
  sd_mat_inv <- diag(1/sds)
  cor_mat <- sd_mat_inv%*%sigma%*%sd_mat_inv
  bias_est <- (theta_thresholds[nb_thresholds]-theta_thresholds)
  sd_plus_bias <- sqrt(sds^2 + bias_est^2)#(est_threshold_pred-theta_thresholds)^2)
  sigma_plus_bias_inv <- solve(diag(sd_plus_bias)%*%cor_mat%*%diag(sd_plus_bias))
  weights <- (one%*%sigma_plus_bias_inv/(one%*%sigma_plus_bias_inv%*%one)[1,1])[1,] # adding bias to the 
  #   
  #   weights <- rep(0, nb_thresholds)
  #   weights[1] <- 0.5
  #   weights[length(weights)] <- 0.5
  
  #cor_mat <- solve(diag(sqrt(diag(sigma))))%*%sigma%*%solve(diag(sqrt(diag(sigma))))
  #sigma_bis <- diag(betahat_all[1:10*3, 2][selected_thres])%*%cor_mat%*%diag(betahat_all[1:10*3, 2][selected_thres])
  
  theta_delta <- theta_thresholds[selected_thres]%*%weights
  se_delta <-   sqrt(weights%*%sigma%*%weights)
  #pval_zscore(theta_delta, se_delta)
  
  
  #replacing standard error by influence functions values in betahat
  se_influ <- sqrt(diag(crossprod(eta_mat)))
  betahat_all[grep("Best.ave",rownames(betahat_all)), "Std. Error"] <- se_influ
  
  #predict(lmlogfit_threshold, se.fit = T, newdata=cbind.data.frame("x"=1))
  #betahat_all <- rbind(betahat_all, est_threshold_pred$fit, est_threshold_pred$se.fit)
  
  #   cut_p <- 0.5
  #   xi_05 <- mydata_match[,names_prob[1]] > cut_p
  #   mydata1 <- mydata_match[xi_05,]
  #   mydata1 <- as.data.frame(cbind(mydata1, "G.impute" = apply(mydata1[,paste0("G",1:k0)]*mydata1[,names_prob],1,sum)/
  #                                    apply(mydata1[,names_prob],1,sum)))
  #   fit_05 <- glm(Pheno~G.impute, data=mydata1, family=binomial)
  #   inv_fishInfo_mat_05 <- summary(fit_05)$cov.scaled
  #   theta_05 <- summary(fit_05)$coef[,"Estimate", drop=FALSE]
  #   eta_05 <- rbind((1*xi_05), (1*xi_05))*inv_fishInfo_mat_05%*%t(Z)%*%diag(x=(mydata_match[,"Pheno"] - expit(Z%*%theta_05)[,1]))
  #   
  #   cut_p <- 0.9
  #   xi_09 <- mydata_match[,names_prob[1]] > cut_p
  #   mydata1 <- mydata_match[xi_09,]
  #   mydata1 <- cbind(mydata1, "G.impute" = apply(mydata1[,paste0("G",1:k0)]*mydata1[,names_prob],1,sum)/
  #                      apply(mydata1[,names_prob],1,sum))
  #   fit_09 <- glm(Pheno~G.impute, data=as.data.frame(mydata1), family=binomial)
  #   inv_fishInfo_mat_09 <- summary(fit_09)$cov.scaled
  #   theta_09 <- summary(fit_09)$coef[,"Estimate", drop=FALSE]
  #   eta_09 <- rbind((1*xi_09), (1*xi_09))*inv_fishInfo_mat_09%*%t(Z)%*%diag(x=(mydata_match[,"Pheno"] - expit(Z%*%theta_09)[,1]))
  #   
  #   sigma <- matrix(c(crossprod(eta_05[2,]), eta_05[2,]%*%eta_09[2,], 
  #                     eta_09[2,]%*%eta_05[2,], crossprod(eta_09[2,])), 
  #                   nrow=2, ncol=2)
  #   #sigma <- cov(cbind(eta_05[2,], eta_09[2,]))
  #   sigma_inv <- solve(sigma)
  #   one <- c(1,1)
  #   weights <- (one%*%sigma_inv/(one%*%sigma_inv%*%one)[1,1])[1,]
  #   theta_delta <- theta_05[2]*weights[1]+theta_09[2]*weights[2]
  #   se_delta <- sqrt(weights[1]^2*sigma[1,1] + weights[2]*sigma[2,2] + 2*weights[1]*weights[2]*sigma[1,2])
  
  
  result <- c(nmatch_hat_1way, accuracy_all, betahat_all, 
              est_threshold_pred, sd(est_threshold_pred_ptb), 
              theta_delta, se_delta, 
              nmatch_final_all, epsilon_truth, n_freq1, n_freq5, n_freq10, sds, bias_est, pval_combined)
  write(result, file=paste0(outputfile_name, ".txt"), append=T, ncol=500)
  
}


