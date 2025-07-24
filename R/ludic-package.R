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
#' Version: \tab ludic 0.2.1\cr
#' Date: \tab 2025-07-24\cr
#' License: \tab \href{https://cran.r-project.org/web/licenses/MIT}{The "MIT License" (MIT)}\cr
#' }
#' The main function of \code{ludic} is \code{\link{recordLink}}.
#' 
#' @author Boris P. Hejblum, Harrison G. Zhang, Tianxi Cai
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
#'\emph{JAMIA}, 28(12):2582-2592, (2021). 
#'\doi{10.1093/jamia/ocab187}.
#'
#' @name ludic-package
#' @aliases ludic
#' 
#' @useDynLib ludic, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
