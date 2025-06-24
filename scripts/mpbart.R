# MPBART with hard splits (Xu)
library(GcompBART)

nb = 1500 # number of burn-in samples
nd = 1500 # number of post-burn-in posterior draws 
nt = 200 # number of trees

Prior_mult = function(p, ntree){
  return(list(nu = p-1+1, 
              V = diag(p-1),
              ntrees = ntree,
              kfac = 2,
              pswap = 0.1, 
              pbd = 0.5, 
              pb = 0.25,
              alpha = 0.95, 
              beta = 2.0,
              nc = 100, 
              minobsnode = 10))
}

Mcmc_mult = function(p,w0 = NULL,sig0 = NULL,nb, nd){
  res = list(sigma0 = diag(p-1), burn = nb, ndraws = nd, nSigDr = 50)  
  if(!is.null(w0)){
    res = append(res, list("w0" = w0), length(res))
  } 
  if(!is.null(sig0)){
    res$sigma0 = sig0
  }
  return(res)
}

fml = "y ~ X1 + X2 + X3 + X4 + X5 + X6"

baseyl = 2 # reference level

p=length(unique(trd$y)) # number of outcome levels

bfit_mult = model_bart(as.formula(fml), data = trd, type = "multinomial",
                       base = baseyl,
                       Prior = Prior_mult(p = p, ntree = nt),
                       Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                       correction = FALSE, Kindo = FALSE, do_alpha2_prior = FALSE)