#'Kernel for matching binary datasets by Probabilistic Patient Record Linkage
#'
#'@param data1 a matrix of  dataset, either a matrix or a dataframe.
#'
#'@param data2 a matrix of original (continuous) dataset, either a matrix or a dataframe.
#'
#'@param ind_match a matrix with 2 columns containing the rownames of the true matches.
#'First column is for \code{data1} rownames and the second corresponds to \code{data2}.
#'
#'@param MAF0 Minor Allele Frequency. Default is 0.2
#'
#'@param OR Odds Ratio. Default is 1.5
#'
#'@return a list containing the 2 datasets with the simulated phenotype and genotype appended 
#'
#'@export

simu_asso <- function(data1, data2, ind_match, 
                      MAF0 = 0.2, OR = 1.5
){
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  nb_match <- nrow(ind_match)
  
  ind1 <- rownames(data1)
  ind2 <- rownames(data2)
  
  tmp_Genom <- rbinom(n1, size=2, prob=MAF0)
  tmp_Pheno <- rbinom(n1, size=1, prob=expit(-0.5+log(OR)*tmp_Genom))
  data1_all <- cbind(data1, "Pheno"=tmp_Pheno, "Genom"=tmp_Genom)
  data1_out <- data1_all[, c("Pheno","Genom")]
  
  tmp_Genom <- c(data1_out[as.character(ind_match[, "data1"]), "Genom"], 
                 rbinom(n2-nb_match, size=2, prob=MAF0))
  tmp_Pheno <- c(data1_out[as.character(ind_match[, "data1"]), "Pheno"], 
                 rbinom(n2-nb_match, size=1, prob=expit(-0.5+log(OR)*tmp_Genom[(nb_match+1):n2])))
  
  notmatch2 <- rownames(data2)[!(rownames(data2) %in% as.character(ind_match[, "data2"]))]
  data2 <- data2[c(as.character(ind_match[, "data2"]), notmatch2), ]
  data2_all <- cbind(data2, "Pheno"=tmp_Pheno, "Genom" = tmp_Genom)
  data2_out <- data2_all[, c("Pheno","Genom")]
  
  return(list("data1" = data1_out, "data2" = data2_out))
}


