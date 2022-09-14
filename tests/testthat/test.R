library(kimfilter)
data(sw_dcf)
data = sw_dcf[, colnames(sw_dcf) != "dcoinc"]
vars = colnames(data)[colnames(data) != "date"]

#Set up the state space model
ssm = list()
ssm[["Fm"]] = rbind(c(0.8760, -0.2171, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0.0364, -0.0008, 0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, -0.2965, -0.0657, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, -0.3959, -0.1903, 0, 0),
                    c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2436, 0.1281), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
ssm[["Fm"]] = array(ssm[["Fm"]], dim = c(dim(ssm[["Fm"]]), 2))
ssm[["Dm"]] = matrix(c(-1.5700, rep(0, 11)), nrow = nrow(ssm[["Fm"]]), ncol = 1)
ssm[["Dm"]] = array(ssm[["Dm"]], dim = c(dim(ssm[["Dm"]]), 2))
ssm[["Dm"]][1,, 2] = 0.2802
ssm[["Qm"]] = diag(c(1, 0, 0, 0, 0.0001, 0, 0.0001, 0, 0.0001, 0, 0.0001, 0))
ssm[["Qm"]] = array(ssm[["Qm"]], dim = c(dim(ssm[["Qm"]]), 2))
ssm[["Hm"]] = rbind(c(0.0058, -0.0033, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                    c(0.0011, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                    c(0.0051, -0.0033, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0), 
                    c(0.0012, -0.0005, 0.0001, 0.0002, 0, 0, 0, 0, 0, 0, 1, 0))
ssm[["Hm"]] = array(ssm[["Hm"]], dim = c(dim(ssm[["Hm"]]), 2))
ssm[["Am"]] = matrix(0, nrow = nrow(ssm[["Hm"]]), ncol = 1)
ssm[["Am"]] = array(ssm[["Am"]], dim = c(dim(ssm[["Am"]]), 2))
ssm[["Rm"]] = matrix(0, nrow = nrow(ssm[["Am"]]), ncol = nrow(ssm[["Am"]]))
ssm[["Rm"]] = array(ssm[["Rm"]], dim = c(dim(ssm[["Rm"]]), 2))
ssm[["B0"]] = matrix(c(rep(-4.60278, 4), 0, 0, 0, 0, 0, 0, 0, 0)) 
ssm[["B0"]] = array(ssm[["B0"]], dim = c(dim(ssm[["B0"]]), 2))
ssm[["B0"]][1:4,, 2] = rep(0.82146, 4)
ssm[["P0"]] = rbind(c(2.1775, 1.5672, 0.9002, 0.4483, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(1.5672, 2.1775, 1.5672, 0.9002, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(0.9002, 1.5672, 2.1775, 1.5672, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(0.4483, 0.9002, 1.5672, 2.1775, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0.0001, 0, 0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0.0001,  0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0.0001, -0.0001, 0, 0, 0, 0),
                    c(0, 0, 0, 0, 0, 0, -0.0001, 0.0001, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, 0.0001, -0.0001, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, -0.0001, 0.0001, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0001, -0.0001), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0001, 0.0001))
ssm[["P0"]] = array(ssm[["P0"]], dim = c(dim(ssm[["P0"]]), 2))
ssm[["Pm"]] = rbind(c(0.8406, 0.0304), 
                    c(0.1594, 0.9696))                        

#Log, difference, and standardize the data
for(i in vars){
  data[, c(i)] = scale(c(NA, diff(log(data[, c(i)]))))
}

# #testthat doesn't recognize data.table
# data[, c(vars) := lapply(.SD, log), .SDcols = c(vars)]
# data[, c(vars) := lapply(.SD, function(x){
#   x - shift(x, type = "lag", n = 1)
# }), .SDcols = c(vars)]
# data[, c(vars) := lapply(.SD, scale), .SDcols = c(vars)]

#Convert the data to an NxT matrix
yt = t(data[, c(vars)])

test_that("kim filter", {
  expect_equal(class(kim_filter(ssm, yt)), "list")
})

ssm[["betaO"]] = matrix(0, nrow = nrow(yt), ncol = 1)
ssm[["betaO"]] = array(ssm[["betaO"]], dim = c(dim(ssm[["betaO"]]), 2))
ssm[["betaS"]] = matrix(0, nrow = nrow(ssm[["Fm"]]), ncol = 1)
ssm[["betaS"]] = array(ssm[["betaS"]], dim = c(dim(ssm[["betaS"]]), 2))
Xo = matrix(rnorm(ncol(yt)), ncol = ncol(yt), nrow = 1)
Xs = matrix(rnorm(ncol(yt)), ncol = ncol(yt), nrow = 1)

test_that("kim filter exo", {
  expect_equal(class(kim_filter(ssm, yt, Xo, Xs)), "list")
})

weights = matrix(rnorm(ncol(yt)), ncol = 1, nrow = ncol(yt))

test_that("kim filter weighted", {
  expect_equal(class(kim_filter(ssm, yt, weight = weights)), "list")
})

test_that("steady state probabilities", {
  expect_equal(dim(ss_prob(ssm[["Pm"]])), c(dim(ssm[["Pm"]])[1], 1))
})
