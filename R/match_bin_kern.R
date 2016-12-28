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
#'@param yes_plot logical flag on whether a plot should be made. Default is FALSE.
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

match_bin_kern <- function(data, nb_rows = NULL, 
                           cut_bin = 0, MAF0 = 0.2, OR=1.3, k0=5,
                           nb_match=NULL, n1=3000, n2=2000, 
                           eps1=0.01, eps2=eps1, bet_ptb=6, 
                           yes_plot=FALSE,
                           type="bin",
                           outputfolder_name=".",
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
  
  n_freq10 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.10 | colSums(data1_bin)/ncol(data1_bin)>0.90)
  n_freq5 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.05 | colSums(data1_bin)/ncol(data1_bin)>0.95)
  n_freq1 <- sum(colSums(data1_bin)/nrow(data1_bin)<0.01 | colSums(data1_bin)/ncol(data1_bin)>0.99)
  
  #mean(epsilon_kminus_truth)
  #mean(epsilon_kplus_truth)
  
  
  #print(cbind(apply(dat.A[id_match,],2,mean),apply(dat.B[id_match,],2,mean)))
  #print(rbind(apply(dat.A[id_match,]*(1-dat.B)[id_match,],2,mean)/apply(dat.A[id_match,],2,mean),
  #            apply((1-dat.A[id_match,])*dat.B[id_match,],2,mean)/apply(1-dat.A[id_match,],2,mean)))
  
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
  
  #lik_all <- rbind("11"=log((1-eps_p)/pi2),"00"=log((1-eps_n)/(1-pi2)),"10"=log(eps_p/(1-pi2)),"01"=log(eps_n/pi2))
  
  #dist_all <- sapply(1:nrow(data2_bin),function(i){loglik_FUN(data1_bin, data2_bin[i,], 
  #                                                            eps_p, eps_n, piA=pi2, piB=pi2)
  #})
  dist_all <- loglikC_mat_faster(Bmat = data2_bin, Amat = data1_bin, eps_p = eps_p, 
                                 eps_n = eps_n, piA = pi2, piB = pi2)
  colnames(dist_all) <- row.names(data2_bin)
  rownames(dist_all) <- row.names(data1_bin)
  
  
  #   dist_all2 <- loglikC_mat(Bmat = data1_bin, Amat = data2_bin, eps_p = eps_n, 
  #                           eps_n = eps_p, piA = pi2, piB = pi2)
  #   dim(dist_all2)
  #   colnames(dist_all2) <- row.names(data1_bin)
  #   rownames(dist_all2) <- row.names(data2_bin)
  
  
  #sstdFit genrerates warnings - to be ignored
  junkfit2 <-  try(sstdFit(x=sample(dist_all,10000)))
  if(inherits(junkfit2, "try-error")){junkfit2 <-  try(sstdFit(x=sample(dist_all,10000)))}
  tmp_est <-  junkfit2$estimate
  for(i in 1:4){
    junkfit2_sub <-  try(sstdFit(x=sample(dist_all,10000)))
    if(inherits(junkfit2_sub, "try-error")){junkfit2_sub <-  try(sstdFit(x=sample(dist_all,10000)))}
    tmp_est <-  rbind(tmp_est, junkfit2_sub$estimate)
  }
  junkfit2$est <-  colMeans(tmp_est)
  
  inter_x <- 0.1
  obs_dist_x <- seq(min(dist_all), max(dist_all), by=inter_x)
  
  
  
  if(yes_plot){
    
    #library(msm)
    #library(fitdistrplus)
    #library(mixsmsn)
    #library(sn)
    
    #junkfit22 = sstdFit(pmin(c(dist_all),20)); print(junkfit22$est)
    #junkfit2 = fitdist(dist_all[dist_all < quantile(dist_all, 1-1/max(n1,n2))], distr="tnorm", method="mle",
    #                    start=list(mean=median(dist_all),sd=sd(dist_all)))
    
    ## ======================================================================== ##
    ## fitting a skewed normal distribution to the observed pair-wise distances ##
    ## ======================================================================== ##
    junkfit1 = snormFit(sample(dist_all,10000)) ##print(junkfit1$par)
    ## =========================================================================== ##
    ## fitting a skewed student-t distribution to the observed pair-wise distances ##
    ## =========================================================================== ##
    ##junkfit2 = sstdFit(c(dist_all)); print(junkfit2$est)
    
    ## =============================================================== ##
    ## fitting a mixture of normal to the observed pair-wise distances ##
    ## =============================================================== ##
    junk = Mclust(sample(dist_all,10000),modelNames="V",G=2); ##print(names(junk)); print(junk$parameters)
    
    ##St.analysis <- smsn.mix(c(dist_all), nu = 12, g = 2, family = "Skew.t")
    
    ## =============================================================== ##
    ## fitting a mixture of skew t to the observed pair-wise distances ##
    ## =============================================================== ##
    source("~/Dropbox/HSPH/matching/BorisMatching/myEMMIXuskew.R")
    mixst <- myfmmst(g=1, dat=matrix(as.vector(dist_all), ncol=1)[sample(1:6000000, size=10000),])
    
    
    
    #pdf(file="Density_Fitted_vs_Observed.pdf", width=9, height=7.5)
    plot(density(dist_all), type="n", main = "Observed vs Fitted Distances", xlab="Observed Distance",ylab="Fitted Density")
    hist(dist_all, prob=T, add=TRUE, col="grey90")
    lines(density(dist_all), lwd=3)
    lines(obs_dist_x, pdf_mix(x=obs_dist_x, pdf=dnorm, w=junk$par$pro, 
                              parameters=cbind(junk$par$mean, sqrt(junk$par$var$sig))),
          col=2,lwd=3)#, lty=2)
    
    lines(obs_dist_x, dst(x=obs_dist_x, xi=res$dp["xi"], omega=res$dp["omega"], alpha=res$dp["alpha"], nu=res$dp["nu"]))
    lines(obs_dist_x, dsnorm(obs_dist_x,mean=junkfit1$par[1],sd=junkfit1$par[2],xi=junkfit1$par[3]),col=3,lwd=3)#, lty=3)
    lines(obs_dist_x, dsstd(obs_dist_x,mean=junkfit2$est[1],sd=junkfit2$est[2],nu=junkfit2$est[3],xi=junkfit2$est[4]),col=4,lwd=3)#, lty=4)
    lines(obs_dist_x,mydfmmst(t(t(obs_dist_x)), known=mixst), col=6, lwd=4)#, lty=5)
    legend("topleft",c("Nonparametric", "Normal mixture","Skew Normal","Skew t", "Skew t mixture"),
           lwd=c(rep(3, 4), 4),col=c(1:4,6),bty="n")#, lty=c(1:4,6))
    #dev.off()
  }
  
  
  #single skew-t parametric estimation
  fitskewt_dens_est <- dsstd(obs_dist_x,
                             mean=junkfit2$est[1], 
                             sd=junkfit2$est[2],
                             nu=junkfit2$est[3],
                             xi=junkfit2$est[4])
  
  #   psstd(35,
  #         mean=junkfit2$est[1], 
  #         sd=junkfit2$est[2],
  #         nu=junkfit2$est[3],
  #         xi=junkfit2$est[4])
  
  #   qsstd(0.95,
  #         mean=junkfit2$est[1], 
  #         sd=junkfit2$est[2],
  #         nu=junkfit2$est[3],
  #         xi=junkfit2$est[4]
  #         )
  
  #   tmp_est <- NULL
  #   for(i in 1:30){
  #     fit_st_est <- st.mple(x=matrix(rep(1, 5000), ncol=1), 
  #                           y=as.vector(sample(dist_all,5000)), 
  #                           opt.method = "Nelder-Mead")
  #     tmp_est <-  rbind(tmp_est, fit_st_est$dp)
  #   }
  #   
  #   fit_st_est <- colMeans(tmp_est)
  #  fitskewt_dens_est <- sn::dst(obs_dist_x, xi=fit_st_est["xi"], omega=fit_st_est["omega"], alpha=fit_st_est["alpha"], nu=fit_st_est["nu"])
  
  
  ##USING THE EMPIRICAL DERIVATIVE
  #   y <- density(dist_all)
  #   d1 <- (y$y[-1]-y$y[-length(y$y)])/(y$x[-1]-y$x[-length(y$x)])
  #   plot(x=y$x[-1], y=d1, type="l", col="blue")
  #   
  #   d2 <- (d1[-1]-d1[-length(d1)])/(y$x[-c(1,2)]-y$x[-c(1,length(y$x))])
  #   lines(x=y$x[-c(1,2)], y=d2, col="red")
  #   y$x[-1][which.min(d1)]
  #   y$x[which.max(abs(d1))]
  
  
  ### select cutoff
  # empirical solution
  peak_fitskewt_dens_est <- max(fitskewt_dens_est)
  which_peak_fitskewt_dens_est <- obs_dist_x[which.max(fitskewt_dens_est)]
  rho1_hat_empirical <- min((obs_dist_x[obs_dist_x >= which_peak_fitskewt_dens_est])[fitskewt_dens_est[obs_dist_x >= which_peak_fitskewt_dens_est]/peak_fitskewt_dens_est < 1/min(n1,n2)])
  
  # empirical derivative of skew-t
  #       d1 <- (fitskewt_dens_est[-1]-fitskewt_dens_est[-length(fitskewt_dens_est)])/(obs_dist_x[-1]-obs_dist_x[-length(obs_dist_x)])
  #       plot(x=obs_dist_x[-1], y=d1, type="l", col="blue", 
  #            xlim=c(25,40), ylim=c(-0.05, 0.05),
  #            xlab="Observed distance",
  #            ylab="Derivative")
  #       obs_dist_x[which.max(abs(d1))]
  #     
  #     nb_d=11
  #     for(j in 2:nb_d){
  #       d1 <- (d1[-1]-d1[-length(d1)])/(obs_dist_x[-c(1:j)]-obs_dist_x[-c(1:(j-1),length(obs_dist_x))])
  #       lines(x=obs_dist_x[-c(1:j)], y=d1, col=j)
  #     }
  #       d2 <- (d1[-1]-d1[-length(d1)])/(obs_dist_x[-c(1,2)]-obs_dist_x[-c(1,length(obs_dist_x))])
  #       lines(x=obs_dist_x[-c(1,2)], y=d2, col="red")
  #       obs_dist_x[-1][which.max(abs(d2))]
  #       d3 <- (d2[-1]-d2[-length(d2)])/(obs_dist_x[-c(1:3)]-obs_dist_x[-c(1:2,length(obs_dist_x))])
  #       lines(x=obs_dist_x[-c(1:3)], y=d3, col="green")
  #       d4 <- (d3[-1]-d3[-length(d3)])/(obs_dist_x[-c(1:4)]-obs_dist_x[-c(1:3,length(obs_dist_x))])
  #       lines(x=obs_dist_x[-c(1:4)], y=d4, col="black")  
  #       legend("topright", paste0("d", 1:4), col=c("blue", "red", "green", "black"), lty=1)
  #     d5 <- (d4[-1]-d4[-length(d4)])/(obs_dist_x[-c(1:5)]-obs_dist_x[-c(1:4,length(obs_dist_x))])
  #     lines(x=obs_dist_x[-c(1:5)], y=d5, col="brown")
  #     d6 <- (d5[-1]-d5[-length(d5)])/(obs_dist_x[-c(1:6)]-obs_dist_x[-c(1:5,length(obs_dist_x))])
  #     lines(x=obs_dist_x[-c(1:6)], y=d6, col="yellow")
  
  # analytical solution
  nu <-  junkfit2$est[3]
  s <-  junkfit2$est[2]
  m <-  junkfit2$est[1]
  xi <-  junkfit2$est[4]
  Beta <- gamma(0.5)/gamma((nu+1)/2)*gamma(nu/2)
  m1 <- 2*sqrt(nu-2)/((nu-1)*Beta)
  mu <- m1*(xi-1/xi)
  sigma <- sqrt((1-m1^2)*(xi^2+1/xi^2)+2*m1^2-1)
  d1max_x <- m-(mu-xi*sqrt((nu-2)/(nu+2)))*s/sigma
  #m-mu*s/sigma
  #m-(mu+1/xi*sqrt((nu-2)/(nu+2)))*s/sigma
  
  d1 <- function(x){
    c <- 2*sigma*gamma((nu+1)/2)/(s*(xi+1/xi)*sqrt(pi*nu-2)*gamma(nu/2))
    d <- -c*(nu+1)*sigma/(xi^2*(nu-2)*s)
    y <- (x-m)/s*sigma+mu
    return(d*y*(1+y^2/(xi^2*(nu-2)))^(-(nu+3)/2))
  }
  
  d2 <- function(x){
    c <- 2*sigma*gamma((nu+1)/2)/(s*(xi+1/xi)*sqrt(pi*nu-2)*gamma(nu/2))
    d <- -c*(nu+1)*sigma/(xi^2*(nu-2)*s)
    y <- (x-m)/s*sigma+mu
    return(d*sigma/s*(1+y^2/(xi^2*(nu-2)))^(-(nu+5)/2)*(1-y^2*(nu+2)/xi^2*(nu-2)))
  }
  
  rho1_hat <- obs_dist_x[obs_dist_x>d1max_x][min(which(abs(d1(obs_dist_x[obs_dist_x>d1max_x]))<(1/min(n1,n2)) & abs(d2(obs_dist_x[obs_dist_x>d1max_x]))<(1/min(n1,n2))))]
  
  dist_all_select <- dist_all > rho1_hat #rho1_hat_empirical
  
  nmatch_hat <- min(sum(colSums(dist_all_select)>=1), sum(rowSums(dist_all_select)>=1))# sum(dist_all_select) # Empirical Bayes prior on P(M=1)
  prop_match <- nmatch_hat/(nrow(data1_bin)*nrow(data2_bin))
  
  
  id_nomatch <- setdiff(row.names(dist_all), id_match)
  names_p_ori <- paste0("p_ori",1:k0)   
  names_p_agg <- paste0("p_agg",1:k0)
  names_ID <- paste0("ID", 1:k0)
  
  for(yes_aggregate in c(FALSE, TRUE)){
    betahat_all <- NULL 
    accuracy_all <- NULL 
    nmatch_final_all <- NULL
    
    rank_match_y <- sapply(id_match, FUN=matchProbs_rank, computed_dist=dist_all, 
                           prop_match=prop_match, yes_aggregate=yes_aggregate, k0 = k0)
    rank_match_n <- sapply(id_nomatch, FUN=matchProbs_rank, computed_dist=dist_all,
                           prop_match=prop_match, yes_aggregate=yes_aggregate, k0 = k0)
    mode(rank_match_y) <- "numeric" 
    mode(rank_match_n) <- "numeric"
    row.names(rank_match_y) <- paste0(rep(c("p_ori","p_agg","logLR","ID"),rep(k0,4)),rep(1:k0,4))
    row.names(rank_match_n) <- row.names(rank_match_y)
    rank_match_y[names_p_agg,] <- rank_match_y[names_p_agg,]*(rank_match_y[names_p_ori,] > 1/k0)
    rank_match_n[names_p_agg,] <- rank_match_n[names_p_agg,]*(rank_match_n[names_p_ori,] > 1/k0)
    
    names_temp <- ifelse(yes_aggregate, "p_agg1", "p_ori1")
    #temp_fn_y <- function(kk){
    #  temp_id <- as.character(rank_match_y[names_ID,kk])
    #  return(is.element(id_match[kk],temp_id[rank_match_y[names_p_agg,kk]>0]))
    #}
    ## use genomic from data2_bin, use phenotype from data1_bin
    
    names_prob <- paste0(ifelse(yes_aggregate ,"p_agg", "p_ori"), 1:k0)
    mydata_match <- data1_out[c(id_match,id_nomatch),]
    temp_id <- as.matrix(cbind(rank_match_y[names_ID,],rank_match_n[names_ID,]))
    mode(temp_id) <- "character"
    mydata_match <- cbind(mydata_match, t(matrix(data2_out[c(temp_id),"Genom"],nrow=k0)),
                          t(cbind(rank_match_y[names_prob,],rank_match_n[names_prob,])))
    mode(mydata_match) <- "numeric" 
    colnames(mydata_match)[2+1:k0] = paste0("G",1:k0)
    mydata_match <- cbind(mydata_match,"id"=t(temp_id))
    mode(mydata_match) = "numeric"
    
    
    mydata0 <- mydata_match[id_match,]
    for(cut_p in c(0.5, 0.75, 0.9)){
      mydata1 <- mydata_match[mydata_match[,names_prob[1]] > cut_p,]
      browser()
      mydata1 <- cbind(mydata1, "G.impute" = (apply(X=mydata1[,paste0("G",1:k0)]*mydata1[,names_prob], MARGIN=1, FUN=sum)/
                                                apply(X=mydata1[,names_prob], MARGIN=1, FUN=sum)))
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
      nmatch.hat_final <-  sum(as.numeric(rank_match_y[names_temp,])>cut_p) + sum(as.numeric(rank_match_n[names_temp,])>cut_p)
      
      betahat_all <- rbind(betahat_all,beta.hat)
      accuracy_all <-  rbind(accuracy_all,accuracy_hat)
      nmatch_final_all <- rbind(nmatch_final_all, nmatch.hat_final)
    }
    
    #browser()
    
    result = c(nmatch_hat, accuracy_all, betahat_all, nmatch_final_all, epsilon_truth, n_freq1, n_freq5, n_freq10)
    final_filename <- paste0(outputfolder_name, ifelse(yes_aggregate, "/", "_noagg/"), outputfile_name, ifelse(yes_aggregate, ".txt", "_noagg.txt"))
    write(result, file=final_filename, append=T, ncol=500)
  }
}