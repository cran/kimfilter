# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' R's implementation of the Moore-Penrose pseudo matrix inverse
#' 
#' @param m matrix
#' @return matrix inverse of m
Rginv <- function(m) {
    .Call(`_kimfilter_Rginv`, m)
}

#' Generalized matrix inverse
#' 
#' @param m matrix
#' @return matrix inverse of m
gen_inv <- function(m) {
    .Call(`_kimfilter_gen_inv`, m)
}

#' Finds the steady state probabilities from a transition matrix
#' mat = |p_11 p_21 ... p_m1|
#'       |p_12 p_22 ... p_m2|
#'       |...            ...|
#'       |p_1m p_2m ... p_mm|
#' where the columns sum to 1
#' 
#' @param mat square SxS matrix of probabilities with column sums of 1. S
#' represents the number of states
#' @return matrix of dimensions Sx1 with steady state probabilities
#' @examples
#' library(kimfilter)
#' Pm = rbind(c(0.8406, 0.0304), 
#'            c(0.1594, 0.9696))
#' ss_prob(Pm)
#' @export
ss_prob <- function(mat) {
    .Call(`_kimfilter_ss_prob`, mat)
}

#' Matrix self rowbind
#' 
#' @param mat matrix
#' @param times integer
#' @return matrix
self_rbind <- function(mat, times) {
    .Call(`_kimfilter_self_rbind`, mat, times)
}

#' Check if list contains a name
#' 
#' @param s a string name
#' @param L a list object
#' @return boolean
contains <- function(s, L) {
    .Call(`_kimfilter_contains`, s, L)
}

#' Kim Filter
#' 
#' @param ssm list describing the state space model, must include names
#' B0 - N_b x 1 matrix, initial guess for the unobserved components 
#' P0 - N_b x N_b matrix, initial guess for the covariance matrix of the unobserved components
#' Dm - N_b x 1 matrix, constant matrix for the state equation
#' Am - N_y x 1 matrix, constant matrix for the observation equation
#' Fm - N_b X p matrix, state transition matrix
#' Hm - N_y x N_b matrix, observation matrix
#' Qm - N_b x N_b matrix, state error covariance matrix
#' Rm - N_y x N_y matrix, state error covariance matrix
#' betaO - N_y x N_o matrix, coefficient matrix for the observation exogenous data
#' betaS - N_b x N_s matrix, coefficient matrix for the state exogenous data
#' Pm - n_state x n_state matrix, state transition probability matrix
#' @param yt N x T matrix of data
#' @param Xo N_o x T matrix of exogenous observation data
#' @param Xs N_s x T matrix of exogenous state 
#' @param weight column matrix of weights, T x 1
#' @param smooth boolean indication whether to run the backwards smoother
#' @return list of cubes and matrices output by the Kim filter
#' @examples
#' #Stock and Watson Markov switching dynamic common factor
#' library(kimfilter)
#' library(data.table)
#' data(sw_dcf)
#' data = sw_dcf[, colnames(sw_dcf) != "dcoinc", with = FALSE]
#' vars = colnames(data)[colnames(data) != "date"]
#' 
#' #Set up the state space model
#' ssm = list()
#' ssm[["Fm"]] = rbind(c(0.8760, -0.2171, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#'                   c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#'                   c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                   c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#'                   c(0, 0, 0, 0, 0.0364, -0.0008, 0, 0, 0, 0, 0, 0), 
#'                   c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
#'                   c(0, 0, 0, 0, 0, 0, -0.2965, -0.0657, 0, 0, 0, 0), 
#'                   c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
#'                   c(0, 0, 0, 0, 0, 0, 0, 0, -0.3959, -0.1903, 0, 0),
#'                   c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0), 
#'                   c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2436, 0.1281), 
#'                   c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
#' ssm[["Fm"]] = array(ssm[["Fm"]], dim = c(dim(ssm[["Fm"]]), 2))
#' ssm[["Dm"]] = matrix(c(-1.5700, rep(0, 11)), nrow = nrow(ssm[["Fm"]]), ncol = 1)
#' ssm[["Dm"]] = array(ssm[["Dm"]], dim = c(dim(ssm[["Dm"]]), 2))
#' ssm[["Dm"]][1,, 2] = 0.2802
#' ssm[["Qm"]] = diag(c(1, 0, 0, 0, 0.0001, 0, 0.0001, 0, 0.0001, 0, 0.0001, 0))
#' ssm[["Qm"]] = array(ssm[["Qm"]], dim = c(dim(ssm[["Qm"]]), 2))
#' ssm[["Hm"]] = rbind(c(0.0058, -0.0033, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
#'                   c(0.0011, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
#'                   c(0.0051, -0.0033, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0), 
#'                   c(0.0012, -0.0005, 0.0001, 0.0002, 0, 0, 0, 0, 0, 0, 1, 0))
#' ssm[["Hm"]] = array(ssm[["Hm"]], dim = c(dim(ssm[["Hm"]]), 2))
#' ssm[["Am"]] = matrix(0, nrow = nrow(ssm[["Hm"]]), ncol = 1)
#' ssm[["Am"]] = array(ssm[["Am"]], dim = c(dim(ssm[["Am"]]), 2))
#' ssm[["Rm"]] = matrix(0, nrow = nrow(ssm[["Am"]]), ncol = nrow(ssm[["Am"]]))
#' ssm[["Rm"]] = array(ssm[["Rm"]], dim = c(dim(ssm[["Rm"]]), 2))
#' ssm[["B0"]] = matrix(c(rep(-4.60278, 4), 0, 0, 0, 0, 0, 0, 0, 0)) 
#' ssm[["B0"]] = array(ssm[["B0"]], dim = c(dim(ssm[["B0"]]), 2))
#' ssm[["B0"]][1:4,, 2] = rep(0.82146, 4)
#' ssm[["P0"]] = rbind(c(2.1775, 1.5672, 0.9002, 0.4483, 0, 0, 0, 0, 0, 0, 0, 0), 
#'                     c(1.5672, 2.1775, 1.5672, 0.9002, 0, 0, 0, 0, 0, 0, 0, 0), 
#'                     c(0.9002, 1.5672, 2.1775, 1.5672, 0, 0, 0, 0, 0, 0, 0, 0), 
#'                     c(0.4483, 0.9002, 1.5672, 2.1775, 0, 0, 0, 0, 0, 0, 0, 0), 
#'                     c(0, 0, 0, 0, 0.0001, 0, 0, 0, 0, 0, 0, 0), 
#'                     c(0, 0, 0, 0, 0, 0.0001,  0, 0, 0, 0, 0, 0), 
#'                     c(0, 0, 0, 0, 0, 0, 0.0001, -0.0001, 0, 0, 0, 0),
#'                     c(0, 0, 0, 0, 0, 0, -0.0001, 0.0001, 0, 0, 0, 0), 
#'                     c(0, 0, 0, 0, 0, 0, 0, 0, 0.0001, -0.0001, 0, 0), 
#'                     c(0, 0, 0, 0, 0, 0, 0, 0, -0.0001, 0.0001, 0, 0), 
#'                     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0001, -0.0001), 
#'                     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0001, 0.0001))
#' ssm[["P0"]] = array(ssm[["P0"]], dim = c(dim(ssm[["P0"]]), 2))
#' ssm[["Pm"]] = rbind(c(0.8406, 0.0304), 
#'                     c(0.1594, 0.9696))
#'   
#' #Log, difference and standardize the data
#' data[, c(vars) := lapply(.SD, log), .SDcols = c(vars)]
#' data[, c(vars) := lapply(.SD, function(x){
#'   x - shift(x, type = "lag", n = 1)
#' }), .SDcols = c(vars)]
#' data[, c(vars) := lapply(.SD, scale), .SDcols = c(vars)]
#'   
#' #Convert the data to an NxT matrix
#' yt = t(data[, c(vars), with = FALSE])
#' kf = kim_filter(ssm, yt, smooth = TRUE)
#' @export
kim_filter <- function(ssm, yt, Xo = NULL, Xs = NULL, weight = NULL, smooth = FALSE) {
    .Call(`_kimfilter_kim_filter`, ssm, yt, Xo, Xs, weight, smooth)
}

