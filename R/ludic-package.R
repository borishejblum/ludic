#' ludic
#' 
#' Linkage Using Diagnosis Codes
#' 
#' This package implements probabilistic record linkage methods that relies on the use of diagnosis 
#' codes only, in the absence of direct identifiers .
#' 
#' \tabular{ll}{
#' Package: \tab ludic\cr
#' Type: \tab Package\cr
#' Version: \tab ludic 0.2.0\cr
#' Date: \tab 2021-08-18\cr
#' License: \tab \href{https://cran.r-project.org/web/licenses/MIT}{The "MIT License" (MIT)}\cr
#' }
#' The main function of \code{ludic} is \code{\link{recordLink}}.
#' 
#' @author Boris P. Hejblum, Tianxi Cai
#' --- Maintainer: Boris P. Hejblum
#' 
#'@references Hejblum BP, Weber G, Liao KP, Palmer N, Churchill S, Szolovits P, 
#'Murphy S, Kohane I and Cai T, Probabilistic Record Linkage of De-Identified 
#'Research Datasets Using Diagnosis Codes, \emph{Scientific Data}, 6:180298 (2019). 
#'\doi{10.1038/sdata.2018.298}.
#'
#'@references Zhang HG, Hejblum BP, Weber G, Palmer N, Churchill S, Szolovits P, 
#'Murphy S, Liao KP, Kohane I and Cai T, ATLAS: An automated association test using 
#'probabilistically linked health records with application to genetic studies, 
#'\emph{JAMIA}, in press (2021). 
#'\doi{10.1101/2021.05.02.21256490}.
#'
#' @docType package
#' @name ludic-package
#' @aliases ludic
#' 
#' @useDynLib ludic, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' 
NULL