#'Adding noise around data
#'
#'@param x original values
#'@param type a character string indicating which type of data is 'x'. Can be 
#'either binary (\code{"bin"} for binomial or \code{"bin_thresh"} for thresholding
#'a multivariate normal), continuous (\code{"con"}), or ... (\code{"dis"})
#'Default ids \code{"bin"}.
#'@param bet guessed value for the level of pertubation added
#'
#' @export
Perturb_FUN = function(x, type=NULL, bet=3){
  mu <- mean(x>0)
  n <- length(x)
  if(is.null(type)){
    K.uni <- length(unique(x))
    if(K.uni==2){
      type <- "bin"
    }
    else if(K.uni/n > 0.05){
      type <- "con"
    }
    else{
      type <- "dis"
    }
  }
  switch(type,
         "bin"={x_perturb <- rbinom(n, size=1, prob=expit(bet*(x-0.5)/sqrt(mu)+0.1*log(mu/(1-mu))))},#1.2*log(mu/(1-mu))))},
         "bin_thresh"={x_perturb <- 1*(rnorm(length(x), x-0.5, sd=0.5/qnorm(1-bet/100))>0)},
         "dis"={x_perturb <- round(exp( (x==0)*log(1-mu+1/n) + log(x+1) + rnorm(length(x))*sd(log(x+1))/20))},
         "con"={x_perturb <- exp(log(abs(x)+0.02*sd(x)) + rnorm(n)*sd(log(abs(x)+0.01*sd(x)))*0.1)*sign(x)}
  )
  return(x_perturb)
}