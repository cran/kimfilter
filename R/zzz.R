#' Kim Filter
#' 
#' \emph{kimfilter} Rcpp implementation of the multivariate Kim filter, which 
#' combines the Kalman and Hamilton filters for state probability inference. 
#' The filter is designed for state space models and can handle missing values 
#' and exogenous data in the observation and state equations.
#' \code{browseVignettes("kimfilter")} to view it in your browser.
#'  
#' @docType package
#' @author Alex Hubbard
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib kimfilter, .registration=TRUE
#' @name kimfilter
NULL