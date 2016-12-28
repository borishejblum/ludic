#'Kernel for matching binary datasets by Probabilistic Patient Record Linkage
#'
#'@param data a matrix of original (continuous) dataset, either a matrix or a dataframe.
#'
#'@param nb_rows the number of row to keep (equal or lower than the number of 
#'row of data). Default is \code{NULL} in which case all the rows are kept
#'
#'@param cut_bin an integer specifying the threshold for dichotomizing observations in data.
#'Default is 0, in which case, any observation bigger than 0 becomes \code{TRUE} 
#'in the simulated data.
#'
#'@param MAF0 Minor Allele Frequency. Default is 0.2
#'
#'@param OR Odds Ratio. Default is 1.3.
#'
#'@param k0 the number of top matching probabilities to store. Default is 5 top probabilities.
#'
#'@param nb_match. The number of true match to be simulated. Default is 3\% of the 
#'simulated pairs of patients. 
#'
#'@param n1 the number of patients in the first database. Default is 3000.
#'
#'@param n2 the number of patients in the second database. Default is 2000.
#'
#'@param eps1 probability that a feature is recorded in the first database but 
#'not in the second (mis-match or recording error). Default is 0.01.
#'
#'@param eps2 probability that a feature is recorded in the second database but 
#'not in the first (mis-match or recording error). Default is the same as \code{eps1}.
#'
#'@param bet_ptb guessed value for the level of pertubation added
#'
#'@param outputfile_name a character string for the name of the .txt file where the 
#'results should be stored. Must be provided without extension. Default is 'Default_output'.
#'
#'@importFrom landpred VTM
#'
#'@importFrom fGarch dsnorm dsstd sstdFit snormFit psstd qsstd
#'@importFrom EMMIXuskew dfmmst
#'@importFrom mclust Mclust mclustBIC
#'@importFrom sn st.mple dst
#'
#'@return append the results in the file '\code{outputfile_name}.txt'.
#'
#'@export

match_2ways <- function(data, nb_rows = NULL, 
                        cut_bin = 0, MAF0 = 0.2, OR=1.3, k0=5,
                        nb_match=NULL, n1=3000, n2=2000, 
                        eps1=0.01, eps2=eps1, bet_ptb=6, 
                        yes_plot=FALSE,
                        type="bin",
                        outputfile_name="Default_output"){
  if(!is.null(nb_rows)){
    data <- data[1:nb_rows,]
  }else{
    nb_rows <- nrow(data)
  }
  row.names(data) <- 1:nb_rows
  
  if(is.null(nb_match)){
    nb_match <- round(0.03*nb_rows)
  }
  
  
  ind1 <- sample(1:nb_rows, n1)
  ind2 <- c(sample(ind1, nb_match), sample(setdiff(1:nb_rows,ind1), n2-nb_match))
  ind_match <- intersect(ind1,ind2) ##print(length(ind_match)); 
  id_match <- as.character(ind_match)
  yes_match0 <- (ind1 == VTM(ind2,length(ind1)))
  
  tmp_Genom <- rbinom(nb_rows, size=2, prob=MAF0)
  tmp_Pheno = c(rbinom(nb_rows, size=1, prob=expit(-0.5+log(OR)*tmp_Genom)))
  data_bin_allp <- 1*(data > cut_bin)
  freq_select <- apply(X=(data_bin_allp==1), MARGIN=2, FUN=mean)>0.01
  data_bin <- data_bin_allp[,freq_select]
  data_all = cbind(data_bin,"Pheno"=tmp_Pheno,"Genom"=tmp_Genom)
  data1_all <- data_all[ind1,]
  data1_out <- data1_all[,c("Pheno","Genom")]
  data1_bin <- data_bin[ind1, ]
  
  if(type=="bin"){
    data2_all <- data_all[ind2,]
    data2_out <- data2_all[,c("Pheno","Genom")]
    data2_bin <- data_bin[ind2, ]
    
    data2_bin <- sapply(1:ncol(data2_bin),function(i){
      Perturb_FUN(data2_bin[, i], type = type, bet = bet_ptb)}
    ) 
    row.names(data2_bin) <- row.names(data2_all)
    colnames(data2_bin) <- colnames(data_bin)
  }
  else if(type=="con"){
    data2 <- data[ind2, freq_select]
    data2_ptb <- sapply(1:ncol(data2),function(i){
      Perturb_FUN(data2[, i], type = "con", bet = bet_ptb)}
    )
    data2_bin <- 1*(data2_ptb > cut_bin)
    data2_all <- cbind(data2_bin,"Pheno"=tmp_Pheno[ind2],"Genom"=tmp_Genom[ind2])    
    data2_out <- data2_all[,c("Pheno","Genom")]
    
    row.names(data2_bin) <- as.character(ind2)
    row.names(data2_out) <- as.character(ind2)
    colnames(data2_bin) <- colnames(data_bin)
  }
  
  
  
  
  epsilon_kminus_truth <- sapply(1:ncol(data2_bin),function(v){
    mean(data2_bin[id_match[data1_bin[id_match,v]==1],v]==0)
  })
  epsilon_kplus_truth <- sapply(1:ncol(data2_bin),function(v){
    mean(data2_bin[id_match[data1_bin[id_match,v]==1],v]==0)
  })
  
  epsilon_truth <- c(mean(epsilon_kminus_truth, na.rm = TRUE), mean(epsilon_kplus_truth, na.rm = TRUE),
                     mean(abs(epsilon_kminus_truth -mean(epsilon_kminus_truth, na.rm = TRUE)), na.rm = TRUE),
                     mean(abs(epsilon_kplus_truth-mean(epsilon_kplus_truth, na.rm = TRUE)), na.rm = TRUE)
  )
  
  
  
  
  #removing info too rare in any of the 2 database ??
  col_in <- (pmax(apply(X=(data1_bin==0), MARGIN=2, FUN=mean), apply(X=data2_bin==0, MARGIN=2, FUN=mean)) < 0.995 &  
               pmin(apply(data1_bin==0,2,mean),apply(data2_bin==0,2,mean)) > 0.005)
  data1_bin <- data1_bin[,col_in]
  data2_bin <- data2_bin[,col_in]
  
  pi1 <- apply(data1_bin==1, MARGIN=2, FUN=mean)
  pi2 <- apply(data2_bin==1, MARGIN=2, FUN=mean)
  K0 <- ncol(data1_bin)
  eps_p <- rep(eps1, K0)
  eps_n <- rep(eps2,K0)
  
  
  
  
  dist_all_1way <- loglikC_mat_faster(Bmat = data2_bin, Amat = data1_bin, eps_p = eps_p, 
                                      eps_n = eps_n, piA = pi2, piB = pi2)
  colnames(dist_all_1way) <- row.names(data2_bin)
  rownames(dist_all_1way) <- row.names(data1_bin)
  
  dist_all_2way <- loglikC_mat(Bmat = data1_bin, Amat = data2_bin, eps_p = eps_n, 
                               eps_n = eps_p, piA = pi2, piB = pi2)
  colnames(dist_all_2way) <- row.names(data1_bin)
  rownames(dist_all_2way) <- row.names(data2_bin)
  
  #sstdFit genrerates warnings - to be ignored
  sstdFit_1way <-  try(sstdFit(x=sample(dist_all_1way,10000)))
  if(inherits(sstdFit_1way, "try-error")){sstdFit_1way <-  try(sstdFit(x=sample(dist_all_1way,10000)))}
  tmp_est <-  sstdFit_1way$estimate
  for(i in 1:4){
    sstdFit_1way_sub <-  try(sstdFit(x=sample(dist_all_1way,10000)))
    if(inherits(sstdFit_1way_sub, "try-error")){sstdFit_1way_sub <-  try(sstdFit(x=sample(dist_all_1way,10000)))}
    tmp_est <-  rbind(tmp_est, sstdFit_1way_sub$estimate)
  }
  sstdFit_1way$est <-  colMeans(tmp_est)
  
  sstdFit_2way <-  try(sstdFit(x=sample(dist_all_2way,10000)))
  if(inherits(sstdFit_2way, "try-error")){sstdFit_2way <-  try(sstdFit(x=sample(dist_all_2way,10000)))}
  tmp_est <-  sstdFit_2way$estimate
  for(i in 1:4){
    sstdFit_2way_sub <-  try(sstdFit(x=sample(dist_all_2way,10000)))
    if(inherits(sstdFit_2way_sub, "try-error")){sstdFit_2way_sub <-  try(sstdFit(x=sample(dist_all_2way,10000)))}
    tmp_est <-  rbind(tmp_est, sstdFit_2way_sub$estimate)
  }
  sstdFit_2way$est <-  colMeans(tmp_est)
  
  
  inter_x <- 0.1
  obs_dist_x_1way <- seq(min(dist_all_1way), max(dist_all_1way), by=inter_x)
  obs_dist_x_2way <- seq(min(dist_all_2way), max(dist_all_2way), by=inter_x)
  
  ### select cutoff
  # empirical solution
  #single skew-t parametric estimation
  fitskewt_dens_est <- dsstd(obs_dist_x_1way,
                             mean=sstdFit_1way$est[1], 
                             sd=sstdFit_1way$est[2],
                             nu=sstdFit_1way$est[3],
                             xi=sstdFit_1way$est[4])
  peak_fitskewt_dens_est <- max(fitskewt_dens_est)
  rho1_hat_empirical_1way <- min((obs_dist_x_1way[obs_dist_x_1way >= peak_fitskewt_dens_est])[fitskewt_dens_est[obs_dist_x_1way >= peak_fitskewt_dens_est]/peak_fitskewt_dens_est < 1/min(n1,n2)])
  fitskewt_dens_est <- dsstd(obs_dist_x_2way,
                             mean=sstdFit_2way$est[1], 
                             sd=sstdFit_2way$est[2],
                             nu=sstdFit_2way$est[3],
                             xi=sstdFit_2way$est[4])
  peak_fitskewt_dens_est <- max(fitskewt_dens_est)
  rho1_hat_empirical_2way <- min((obs_dist_x_2way[obs_dist_x_2way >= peak_fitskewt_dens_est])[fitskewt_dens_est[obs_dist_x_2way >= peak_fitskewt_dens_est]/peak_fitskewt_dens_est < 1/min(n1,n2)])
  
  
  # analytical solution
  nu <-  sstdFit_1way$est[3]
  s <-  sstdFit_1way$est[2]
  m <-  sstdFit_1way$est[1]
  xi <-  sstdFit_1way$est[4]
  Beta <- gamma(0.5)/gamma((nu+1)/2)*gamma(nu/2)
  m1 <- 2*sqrt(nu-2)/((nu-1)*Beta)
  mu <- m1*(xi-1/xi)
  sigma <- sqrt((1-m1^2)*(xi^2+1/xi^2)+2*m1^2-1)
  rho1_hat_1way <- m-(mu-xi*sqrt((nu-2)/(nu+2)))*s/sigma
  
  nu <-  sstdFit_2way$est[3]
  s <-  sstdFit_2way$est[2]
  m <-  sstdFit_2way$est[1]
  xi <-  sstdFit_2way$est[4]
  Beta <- gamma(0.5)/gamma((nu+1)/2)*gamma(nu/2)
  m1 <- 2*sqrt(nu-2)/((nu-1)*Beta)
  mu <- m1*(xi-1/xi)
  sigma <- sqrt((1-m1^2)*(xi^2+1/xi^2)+2*m1^2-1)
  rho1_hat_2way <- m-(mu-xi*sqrt((nu-2)/(nu+2)))*s/sigma
  
  
  
  dist_all_1way_select <- dist_all_1way > rho1_hat_1way
  dist_all_2way_select <- dist_all_2way > rho1_hat_2way
  
  nmatch_hat_1way <- sum(dist_all_1way_select) # Empirical Bayes prior on P(M=1)
  prop_match_1way <- nmatch_hat_1way/(nrow(data1_bin)*nrow(data2_bin))
  nmatch_hat_2way <- sum(dist_all_2way_select) # Empirical Bayes prior on P(M=1)
  prop_match_2way <- nmatch_hat_2way/(nrow(data1_bin)*nrow(data2_bin))
  
  id_nomatch_1 <- setdiff(row.names(dist_all_1way), id_match)
  id_nomatch_2 <- setdiff(row.names(dist_all_2way), id_match)
  names_p_ori <- paste0("p_ori",1:k0)   
  names_p_agg <- paste0("p_agg",1:k0)
  names_ID <- paste0("ID", 1:k0)
  
  
  browser()
  for(yes_aggregate in c(FALSE, TRUE)){
    betahat_all <- NULL 
    accuracy_all <- NULL 
    nmatch_final_all <- NULL
    
    rank_match_y_1way <- sapply(id_match, FUN=matchProbs_rank, computed_dist=dist_all_1way, 
                                prop_match=prop_match_1way, yes_aggregate=yes_aggregate, k0 = k0)
    rank_match_n_1way <- sapply(id_nomatch_1, FUN=matchProbs_rank, computed_dist=dist_all_1way,
                                prop_match=prop_match_1way, yes_aggregate=yes_aggregate, k0 = k0)
    rank_match_y_2way <- sapply(id_match, FUN=matchProbs_rank, computed_dist=dist_all_2way, 
                                prop_match=prop_match_2way, yes_aggregate=yes_aggregate, k0 = k0)
    rank_match_n_2way <- sapply(id_nomatch_2, FUN=matchProbs_rank, computed_dist=dist_all_2way,
                                prop_match=prop_match_2way, yes_aggregate=yes_aggregate, k0 = k0)
    mode(rank_match_y_1way) <- "numeric" 
    mode(rank_match_n_1way) <- "numeric"
    mode(rank_match_y_2way) <- "numeric" 
    mode(rank_match_n_2way) <- "numeric"
    row.names(rank_match_y_1way) <- paste0(rep(c("p_ori","p_agg","logLR","ID"),rep(k0,4)),rep(1:k0,4))
    row.names(rank_match_n_1way) <- row.names(rank_match_y_1way)
    row.names(rank_match_y_2way) <- paste0(rep(c("p_ori","p_agg","logLR","ID"),rep(k0,4)),rep(1:k0,4))
    row.names(rank_match_n_2way) <- row.names(rank_match_y_2way)
    rank_match_y_1way[names_p_agg,] <- rank_match_y_1way[names_p_agg,]*(rank_match_y_1way[names_p_ori,] > 1/k0)
    rank_match_n_1way[names_p_agg,] <- rank_match_n_1way[names_p_agg,]*(rank_match_n_1way[names_p_ori,] > 1/k0)
    rank_match_y_2way[names_p_agg,] <- rank_match_y_2way[names_p_agg,]*(rank_match_y_2way[names_p_ori,] > 1/k0)
    rank_match_n_2way[names_p_agg,] <- rank_match_n_2way[names_p_agg,]*(rank_match_n_2way[names_p_ori,] > 1/k0)
    
    names_temp <- ifelse(yes_aggregate, "p_agg1", "p_ori1")
    #temp_fn_y <- function(kk){
    #  temp_id <- as.character(rank_match_y[names_ID,kk])
    #  return(is.element(id_match[kk],temp_id[rank_match_y[names_p_agg,kk]>0]))
    #}
    ## use genomic from data2_bin, use phenotype from data1_bin
    names_prob <- paste0(ifelse(yes_aggregate ,"p_agg", "p_ori"), 1:k0)
    
    mydata_match_1 <- data1_out[c(id_match,id_nomatch_1),]
    mydata_match_2 <- data2_out[c(id_match,id_nomatch_2),]
    temp_id_1 <- as.matrix(cbind(rank_match_y_1way[names_ID,], rank_match_n_1way[names_ID,]))
    temp_id_2 <- as.matrix(cbind(rank_match_y_2way[names_ID,], rank_match_n_2way[names_ID,]))
    mode(temp_id_1) <- "character"
    mode(temp_id_2) <- "character"
    mydata_match_1 <- cbind(mydata_match_1, t(matrix(data2_out[c(temp_id_1),"Genom"],nrow=k0)),
                            t(cbind(rank_match_y_1way[names_prob,],rank_match_n_1way[names_prob,])))
    mydata_match_2 <- cbind(mydata_match_2, t(matrix(data1_out[c(temp_id_2),"Genom"],nrow=k0)),
                            t(cbind(rank_match_y_2way[names_prob,],rank_match_n_2way[names_prob,])))
    mode(mydata_match_1) <- "numeric" 
    mode(mydata_match_2) <- "numeric"
    colnames(mydata_match_1)[2+1:k0] = paste0("G",1:k0)
    colnames(mydata_match_2)[2+1:k0] = paste0("G",1:k0)
    mydata_match_1 <- cbind(mydata_match_1,"id"=t(temp_id_1))
    mydata_match_2 <- cbind(mydata_match_2,"id"=t(temp_id_2))
    mode(mydata_match_1) = "numeric"
    mode(mydata_match_2) = "numeric"
    
    
    mydata0_1 <- mydata_match_1[id_match,]
    mydata0_2 <- mydata_match_2[id_match,]
    for(cut_p in c(0.5, 0.75, 0.9)){
      mydata1 <- mydata_match[mydata_match[,names_prob[1]] > cut_p,]
      mydata1 <- cbind(mydata1, "G.impute" = apply(mydata1[,paste0("G",1:k0)]*mydata1[,names_prob],1,sum)/
                         apply(mydata1[,names_prob],1,sum))
      temp_1 <- try(summary(glm(mydata1[,"Pheno"]~mydata1[,"G1"],family=binomial))$coef[2,-3], silent = TRUE)
      if(inherits(temp_1, "try-error")){temp_1 <- rep(NA, 3)}
      temp_ave <- try(summary(glm(mydata1[,"Pheno"]~mydata1[,"G.impute"],family=binomial))$coef[2,-3], silent = TRUE)
      if(inherits(temp_ave, "try-error")){temp_ave <- rep(NA, 3)}
      beta.hat <- rbind("Oracle"=summary(glm(mydata0[,"Pheno"]~mydata0[,"Genom"],family=binomial))$coef[2,-3],
                        "Best.1"=temp_1,
                        "Best.ave"=temp_ave)
      temp_ppv <- try(mean(apply(row.names(mydata1)==mydata1[,names_ID],1,sum)), silent = TRUE)
      if(inherits(temp_ppv, "try-error")){temp_ppv <- NA}
      accuracy_hat <- c("TPR"=sum(as.numeric(rank_match_y[names_temp,])>cut_p)/ncol(rank_match_y),
                        "FPR"=sum(as.numeric(rank_match_n[names_temp,])>cut_p)/ncol(rank_match_n),
                        "PPV"=temp_ppv,
                        "NPV"=mean(is.element(setdiff(row.names(data1_bin),row.names(mydata1)),id_nomatch))
      )
      #print(accuracy_hat)
      #print(beta.hat)
      nmatch.hat_final <-  sum(as.numeric(rank_match_y[names_temp,])>cut_p) + sum(as.numeric(rank_match_n[names_temp,])>cut_p)
      
      betahat_all <- rbind(betahat_all,beta.hat)
      accuracy_all <-  rbind(accuracy_all,accuracy_hat)
      nmatch_final_all <- rbind(nmatch_final_all, nmatch.hat_final)
    }
    
    result = c(nmatch_hat, accuracy_all, betahat_all, nmatch_final_all, epsilon_truth)
    write(result, file=ifelse(yes_aggregate, paste0(outputfile_name,".txt"), paste0(outputfile_name,"_noagg.txt")),
          append=T, ncol=500)
  }
}