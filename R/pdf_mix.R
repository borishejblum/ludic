#'Mixture of probability density functions
#'
#'@param x quantile where the probability density is to be evaluated.
#'
#'@param pdf probability density function to be used that takes x as its first argument.
#'
#'@param w a vector of lenght K containing the mixture component weights. Default is \code{NULL}, in which case equal 
#'proportions are assumed.
#'
#'@param parameters a matrix of dimension K x p where p is the number of paramaters 
#'to be passed to the pdf function. The columns must be in the same order as the parameters
#'
#'@export

pdf_mix <- function(x, pdf, w=NULL, parameters){
  
  K <- nrow(parameters)
 
  if(is.null(w)){
    w <- rep(1, K)
  }else{
    if(length(w)!=K){
      stop("length of 'w' doesn't match the number of rows of parameters")
    }
  }
  
  w <- w/sum(w)
  
  param_list <- apply(parameters, MARGIN=1, FUN=function(par){lapply(par, "[")})
  
  out <- 0
  for(k in 1:K){
    out <- out + do.call(pdf, c(list(x), param_list[[k]]))*w[k]
  }
  
  return(out)
}


