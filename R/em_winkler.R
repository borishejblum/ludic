#'Implementation of Winkler's EM algorithmf or Fellegi-Sunter matching method
#'
#'@param data1 either a binary matrix or dataframe whose rownames are .
#'
#'@param data2 either a binary matrix or a dataframe whose rownames are .
#'
#'@param tol tolerance for the EM algorithm convergence
#'
#'@param maxit maximum number of iterations for the EM algorithm
#'
#'@param do_plot a logical flag indicating whether a plot should be drawn for the EM convergence. 
#'Default is \code{TRUE}
#'
#'@param verbose a logical flag indicating wehther intermediate values from the EM algorithm should 
#'be printed. Useful for debugging. Default is \code{FALSE}. 
#'
#'
#'@return a matrix of size \code{n1 x n2} with the matching score for each \code{n1*n2} pair.
#'
#'@references
#'Winkler WE. Using the EM Algorithm for Weight Computation in the Fellegi-Sunter Model of Record Linkage. Proc Sect Surv Res Methods, Am Stat Assoc 1988: 667-71.
#'
#'Grannis SJ, Overhage JM, Hui S, et al. Analysis of a probabilistic record linkage technique without human review. AMIA 2003 Symp Proc 2003: 259-63.
#'
#'@examples
#'#mat1 <- matrix(round(rnorm(n=1000, sd=1.2)), ncol=10, nrow=100)
#'#mat2 <- rbind(mat1[1:10, ],
#'#              matrix(round(rnorm(n=900, sd=1.2)), ncol=10, nrow=90)
#'#              )
#'mat1 <- matrix(round(rnorm(n=1000, sd=12)), ncol=10, nrow=100)
#'mat2 <- rbind(mat1[1:10, ],
#'              matrix(round(rnorm(n=900, sd=12)), ncol=10, nrow=90)
#'              )
#'rownames(mat1) <- paste0("A", 1:nrow(mat1))
#'rownames(mat1) <- paste0("B", 1:nrow(mat1))
#'mat1 <- 1*(mat1>1)
#'mat2 <- 1*(mat2>1)
#'res <- em_winkler(mat1, mat2)
#'res <- em_winkler(mat1, mat1)
#'plot(density(as.vector(res$matchingScore)))
#'
#'
#'@export
em_winkler<-function(data1, data2, tol=0.01, maxit=15, do_plot=TRUE, verbose=FALSE){#}, sens_thres=0.99, spec_thres=0.99){
  
  stopifnot(class(data1)=="matrix")
  stopifnot(class(data2)=="matrix")
  stopifnot(data1%in%c(0,1))
  stopifnot(data2%in%c(0,1))
  
  n1 <- nrow(data1)
  rownames_data1 <- rownames(data1)
  n2 <- nrow(data2)
  rownames_data2 <- rownames(data2)
  K <- ncol(data1)
  stopifnot(ncol(data2)==K)
  
  #agreement
  cat("Computing agreement... ")
  ag_score <- try(agree_C(data1, data2))
  if(inherits(ag_score, "try-error")){
    stop("Data may be to big, try using the 'em_winkler_big' function instead...")
  }
  #rm(data1, data2)
  cat("DONE !\n")
  N <- nrow(ag_score)
  
  #temp
  #ag_score <- 1*(ag_score==1)
  
  #EM algorithm
  #initialization
  p <- 0.01 
  m <- rep(0.95, K)
  u <- rep(0.05, K)
  
  #ag_score_neg <- Matrix::Matrix(1-ag_score, sparse=TRUE)
  for(i in 1:maxit){
    cat("it:", i, "\n")
    p_old <- p
    m_old <- m
    u_old <- u
    
    # E step
    # temp <- estep_C_vect_sparse(ag_score, p=p, m=m, u=u) # estep_C_vect(ag_score, p=p, m=m, u=u) # estep_vect(ag_score, p=p, m=m, u=u)
    temp_log <- estep_C_vect(ag_score, p=p, m=m, u=u, Log=TRUE)
    g_m_log <- temp_log[, 1]
    g_u_log <- temp_log[, 2]
    #rm(temp_log)
    # M step
    g_m_log_max <- max(g_m_log)
    g_u_log_max <- max(g_u_log)
    gm_summed_log <- g_m_log_max + log(sum(exp(g_m_log - g_m_log_max)))
    gu_summed_log <- g_u_log_max + log(sum(exp(g_u_log - g_u_log_max)))
    p <- exp(gm_summed_log - log(N))
    m <- (exp(g_m_log - gm_summed_log)%*%ag_score)[1,]
    u <- (exp(g_u_log - gu_summed_log)%*%ag_score)[1,]
    
    if(length(which(m > (1-tol)))>0){
      m[which(m > (1-tol))] <- (1-tol)
    }
    if(length(which(m < tol))>0){
      m[which(m < tol)] <- tol
    }
    if(length(which(u > (1-tol)))>0){
      u[which(u > (1-tol))] <- (1-tol)
    }
    if(length(which(u < tol))>0){
      u[which(u < tol)] <- tol
    }
    if(length(which(p > (1-tol)))>0){
      p[which(p > (1-tol))] <- (1-tol)
    }
    if(length(which(p < tol))>0){
      p[which(p < tol)] <- tol
    }
    
    if(verbose){
      cat("p:", p, "\n")
      cat("m:", m, "\n")
      cat("u:", u, "\n\n")
    }
    
    conv <- (abs(p_old-p) < tol) && all(abs(m_old-m) < tol) && all(abs(u_old-u) < tol)
    if(is.na(conv)){
      p <- p_old
      m <- m_old
      u <- u_old
      conv_flag <- 1
      break()
    }
    if(conv){
      conv_flag <- 0
      break()
    }
  }
  
  if(i==maxit){
    conv_flag <- 2
  }
  cat("EM finished in", i, "iterations.\n p:", p,"\n\n")
  
  # Final matching score
  ms <- matchingScore_C(ag_score, m, u, nA=n1, nB=n2)   #ms_temp <- apply(ag_score, 1, matching_score, m, u)
  rm(ag_score)
  ms_vec <- as.vector(ms)
  o <- order(ms_vec)
  n_truematches <- round(N*p)
  threshold <- ms_vec[rev(o)][n_truematches]
  rownames(ms) <- rownames_data1
  colnames(ms) <- rownames_data2
  
  # TPR_Grannis <- sapply(ag_score, function(x){prod(m^x*(1-m)^(1-x))})#exp(ag_score%*%log(m) + (1-ag_score)%*%log(1-m))
  # TNR_Grannis <- sapply(ag_score, function(x){prod(u^x*(1-u)^(1-x))})#exp(ag_score%*%log(u) + (1-ag_score)%*%log(1-u))
  # if(do_plot){
  #   plot(y=TPR_Grannis[o], x=ms_vec[o], type="l", col="blue", 
  #        ylab="Estimated rate", xlab="Observed matching score",
  #        ylim=c(0,1))
  #   lines(y=TNR_Grannis[o], x=ms_vec[o], type="l", col="red")
  #   legend("right", c("TPR", "TNR"), col=c("blue", "red"), lty=c(1, 1))
  # }
  # browser()
  # 
  # #which(ms>0, arr.ind=TRUE)
  # #ag_score[which(ms_vec>0),]
  # thresholds <- c(ms_vec[o][which(TPR_Grannis[o]>=sens_thres)[1]],
  #                 ms_vec[o][which(TNR_Grannis[o]>=spec_thres)[1]])
  
  
  
  return(list("matchingScore"=ms, "threshold_ms"=threshold, "n_truematches"=n_truematches, "convergence_status"=conv_flag))
}





#'Implementation of Winkler's EM algorithm or Fellegi-Sunter matching method when the data are too 
#'big to compute the agreement matrix.
#'
#'Agreement is recomputed on the fly each time it is needed. The EM steps are completely done in C++
#'This decreases the RAM usage (still important though), at the cost of increasing computationnal time.
#'
#'@export
em_winkler_big<-function(data1, data2, tol=0.01, maxit=50, do_plot=TRUE, verbose=FALSE){#}, sens_thres=0.99, spec_thres=0.99){
  
  stopifnot(class(data1)=="matrix")
  stopifnot(class(data2)=="matrix")
  stopifnot(data1%in%c(0,1))
  stopifnot(data2%in%c(0,1))
  
  cat("\n**********\nATTENTION: this is the version of Winkler's EM algorithm implemented for 'large' datasets, when the agreement matrix cannot be stored in memory (matrix of size n1*n2 x K)...\n*********\n")
  
  n1 <- nrow(data1)
  rownames_data1 <- rownames(data1)
  n2 <- nrow(data2)
  rownames_data2 <- rownames(data2)
  K <- ncol(data1)
  stopifnot(ncol(data2)==K)
  
  N <- n1*n2
  
  #EM algorithm
  #initialization
  p <- 0.5 
  m <- rep(0.9, K)
  u <- rep(0.1, K)
  
  for(i in 1:maxit){
    cat("it:", i, "\n")
    p_old <- p
    m_old <- m
    u_old <- u
    
    # EM steps in C
    EM_res <- EMstep_C_sparse_big(data1, data2, p=p, m=m, u=u)
    # E_res <- Estep_C_sparse_big(data1, data2, p=p, m=m, u=u)
    # M_res <- Mstep_C_sparse_big(data1, data2, g_m=E_res[, 1], g_u=E_res[, 2])
    p <- EM_res[["p"]]
    m <- EM_res[["m"]]
    u <- EM_res[["u"]]

    
    if(length(which(m > (1-tol)))>0){
      m[which(m > (1-tol))] <- (1-tol)
    }
    if(length(which(m < tol))>0){
      m[which(m < tol)] <- tol
    }
    if(length(which(u > (1-tol)))>0){
      u[which(u > (1-tol))] <- (1-tol)
    }
    if(length(which(u < tol))>0){
      u[which(u < tol)] <- tol
    }
    if(length(which(p > (1-tol)))>0){
      p[which(p > (1-tol))] <- (1-tol)
    }
    if(length(which(p < tol))>0){
      p[which(p < tol)] <- tol
    }
    
    if(verbose){
      cat("p:", p, "\n")
      cat("m:", m, "\n")
      cat("u:", u, "\n\n")
    }
    
    conv <- (abs(p_old-p) < tol) && all(abs(m_old-m) < tol) && all(abs(u_old-u) < tol)
    if(is.na(conv)){
      p <- p_old
      m <- m_old
      u <- u_old
      conv_flag <- 1
      break()
    }
    if(conv){
      conv_flag <- 0
      break()
    }
  }
  
  if(i==maxit){
    conv_flag <- 2
  }
  cat("EM finished.\n p:", p,"\n\n")
  
  # Final matching score
  ms <- matchingScore_C_sparse_big(data1, data2, m, u)   #ms_temp <- apply(ag_score, 1, matching_score, m, u)
  
  ms_vec <- as.vector(ms)
  o <- order(ms_vec)
  n_truematches <- round(N*p)
  threshold <- ms_vec[rev(o)][n_truematches]
  rownames(ms) <- rownames_data1
  colnames(ms) <- rownames_data2
  
  return(list("matchingScore"=ms, "threshold_ms"=threshold, "n_truematches"=n_truematches, "convergence_status"=conv_flag))
}





#' #'@keywords internal
#' estep_vect <- function(ag_score, p, m, u){
#'   a <-exp(log(p) + ag_score%*%log(m) + (1-ag_score)%*%log(1-m))
#'   b <- exp(log(1-p) + ag_score%*%log(u) + (1-ag_score)%*%log(1-u))
#'   return(cbind(a/(a+b), b/(a+b)))
#' }
#' 
#' #'@keywords internal
matching_score <- function(agree, m, u){
  sum((log(m)-log(u))^agree*(log(1-m)-log(1-u))^(1-agree))
}