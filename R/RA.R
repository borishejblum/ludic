#'Anonymized binarized diagnosis codes from RA study.
#'
#'An anonymized version of the binarized diagnosis code data from the RA1 and RA2 
#'datasets, over both 6-year and 11-year time span.
#'
#'@name RA
#'@rdname RA
#'@aliases RA1_6y RA2_6y RA1_11y RA2_11y silverstandard_truematches
#'
#'@usage data(RA)
#'
#'@format 5 objects\itemize{
#'      \item{\code{RA1_6y}:} an integer matrix of 0s and 1s containing 4,936 
#'      renamed diagnosis codes for 26,681 patients from the dataset RA1 recorded 
#'      over a 6-year time span.
#'      
#'      \item{\code{RA2_6y}:} an integer matrix of 0s and 1s containing 4,936 
#'      renamed diagnosis codes for 5,707 patients from the dataset RA2 recorded 
#'      over a 6-year time span.
#'
#'      \item{\code{RA1_11y}:} an integer matrix of 0s and 1s containing 5,593 
#'      renamed diagnosis codes for 26,687 patients from the dataset RA1 recorded 
#'      over a 11-year time span.
#'      
#'      \item{\code{RA2_11y}:} an integer matrix of 0s and 1s containing 5,593 
#'      renamed diagnosis codes for 6,394 patients from the dataset RA2 recorded 
#'      over a 11-year time span.
#'      
#'      \item{\code{silverstandard_truematches}:} a character matrix with two 
#'      columns containing the identifiers of the 3,831 pairs of silver-standard 
#'      matches.
#'      
#'}
#'
#'@details
#'The ICD-9 diagnosis codes have also been masked and randomly reordered, replaced by 
#'meaningless names. Finally, the silver-standard matching pairs are also provided to 
#'allow the benchmarking of methods for probabilistic record linkage using diagnosis codes.
#'
#'@references Hejblum BP, Weber G, Liao KP, Palmer N, Churchill S, Szolovits P, 
#'Murphy S, Kohane I and Cai T, Probabilistic Record Linkage of De-Identified 
#'Research Datasets Using Diagnosis Codes, \emph{Scientific Data}, 6:180298 (2019). 
#'\doi{10.1038/sdata.2018.298}.

#'@references Liao, K. P. et al. Electronic medical records for discovery research 
#'in rheumatoid arthritis. \emph{Arthritis Care & Research} 62, 1120-1127 (2010). 
#'\doi{10.1002/acr.20184}
#'
#'@references Liao, K. P. et al. Methods to Develop an Electronic Medical Record 
#'Phenotype Algorithm to Compare the Risk of Coronary Artery Disease across 3 
#'Chronic Disease Cohorts. \emph{PLoS ONE} 10, e0136651 (2015). 
#'\doi{10.1371/journal.pone.0136651}
#'
#'@examples
#'
#'if(interactive()){
#'rm(list=ls())
#'library(ludic)
#'data(RA)
#'res_match_6y <- recordLink(data1 = RA1_6y, data2 = RA2_6y, 
#'                           eps_plus = 0.01, eps_minus = 0.01,
#'                           aggreg_2ways ="mean",
#'                           min_prev = 0,
#'                           use_diff = FALSE)
#'
#'res_match_11y <- recordLink(data1 = RA1_11y, data2 = RA2_11y, 
#'                            eps_plus = 0.01, eps_minus = 0.01,
#'                            aggreg_2ways ="mean",
#'                            min_prev = 0,
#'                            use_diff = FALSE)
#'
#'
#'print.res_matching <- function(res, threshold=0.9, ref=silverstandard_truematches){
#'  have_match_row <- rowSums(res>threshold)
#'  have_match_col <- colSums(res>threshold)
#'  bestmatched_pairs_all <- cbind.data.frame(
#'    "D1"=rownames(res)[apply(res[,which(have_match_col>0), drop=FALSE], 2, which.max)],
#'    "D2"=names(have_match_col)[which(have_match_col>0)]
#'  )
#'  nTM_all <- nrow(ref)
#'  nP_all <- nrow(bestmatched_pairs_all)
#'  TPR_all <- sum(apply(bestmatched_pairs_all, 1, paste0, collapse="") 
#'                 %in% apply(ref, 1, paste0, collapse=""))/nTM_all
#'  PPV_all <- sum(apply(bestmatched_pairs_all, 1, paste0, collapse="") 
#'                 %in% apply(ref, 1, paste0, collapse=""))/nP_all
#'  cat("threshold: ", threshold, 
#'      "\nnb matched: ", nP_all,"; nb true matches: ", nTM_all, 
#'      "\nTPR: ", TPR_all, ";   PPV: ", PPV_all, "\n\n", sep="")
#'}
#'print.res_matching(res_match_6y)
#'print.res_matching(res_match_11y)
#'
#'}
#'
#' @keywords datasets
#' @docType data
NULL