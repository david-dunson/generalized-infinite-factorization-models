#---------------------------------------------------------------------------------------------#
# ---- SIS for Gaussian data ----
# structured factorization of a Gaussian data matrix Y 
# Copyright: Lorenzo Schiavon 2021,
# Licensed under Creative Commons Attribution-NonCommercial 4.0
#---------------------------------------------------------------------------------------------#

Mcmc_SIS = function(Y, X=NULL, as, bs, alpha, a_theta, b_theta, b_beta = 1,  
                    continuous_X = NULL ,  kinit = NULL, kmax = NULL,
                    b0=1, b1=5*10^(-4), start_adapt = 50, 
                    nrun, burn=round(nrun/4), thin = 1, 
                    output = "all", my_seed = 6784, std=T){
  
  # set seed
  set.seed(my_seed)
  
  # dimensions
  n = dim(Y)[1]
  p = dim(Y)[2]
  if(is.null(X)){X=matrix(1, nrow = p, ncol = 1)}
  q = dim(X)[2]
  
  if(is.null(continuous_X)){continuous_X = c(F,rep(T, ncol(X)-1))}
  
  # useful stuffs
  if(is.null(kinit)) kinit = min(floor(log(p)*8), p)
  if(is.null(kmax)) kmax = p+1
  if(kmax<kinit){stop("kmax must be greater or equal than kinit.")}
  
  if(any(output %in% "all")){output = c("covMean", "covSamples", "loadSamples", "loadMean", "loadAbsMean",
                                        "numFactors", "time", "factSamples", "factMean", "locMean", 
                                        "locSamples", "coefMean", "coefSamples")}
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  # scale X continuous
  if(sum(continuous_X)>0){
    X[,continuous_X] = scale(X[,continuous_X])
  }
  
  if(std){
    VY= apply(Y, 2, var)                  # explicitly preserve scale
    Y = scale(Y)
  }else{
    VY= rep(1, p)
  }
  scaleMat = sqrt((VY) %*% t(VY))
  
  
  k = kinit                             # no. of factors to start with
  kstar = k-1
  
  
  # --- Adaptive prob --- #
  prob = 1/exp(b0 + b1*seq(1,nrun))               
  uu = runif(nrun)
  
  
  # --- Initial values --- #
  p_constant =2*exp(1)*log(p)/p          # factor probability constant 
  
  ps = rgamma(p, as, bs )  
  Sigma=diag(1/ps)                                # Sigma = diagonal residual covariance
  
  Lambda_star = matrix(rnorm(p*k), nrow = p, ncol = k)
  eta = matrix(rnorm(n*k), nrow = n, ncol = k)    # factor loadings & latent factors
  
  Beta = matrix(rnorm(q*k), nrow = q, ncol = k)   # local shrinkage coefficients
  pred = X%*%Beta
  logit = plogis(pred)
  Phi = matrix(rbinom(p*k, size = 1, prob = p_constant), nrow = p, ncol = k)
  
  v = c( rbeta(k-1, shape1 = 1, shape2 = alpha), 1)
  w = v*c(1,cumprod(1-v[-k]))                     # weights
  z = rep(k,k)
  rho = rep(1, k)
  
  Lambda = t(t(Lambda_star)* sqrt(rho)) * Phi
  
  Plam = diag( rgamma(k, a_theta, b_theta) )    # precision matrix of lambda star
  
  
  # --- Allocate output object memory --- #
  if(any(output %in% "covMean")) COVMEAN = matrix(0, nrow = p, ncol = p)
  if(any(output %in% "covSamples")) OMEGA = array(dim = c(p, p, sp))
  if(any(output %in% "loadMean")) LOADMEAN = matrix(0, nrow = p, ncol = kmax)
  if(any(output %in% "loadAbsMean")) LOADABSMEAN = matrix(0, nrow = p, ncol = kmax)
  if(any(output %in% "loadSamples")) LAMBDA = list()
  if(any(output %in% "locMean")) LOCMEAN = matrix(0, nrow = p, ncol = kmax)
  if(any(output %in% "locSamples")) PHI = list()
  if(any(output %in% "coefMean")) COEFMEAN = matrix(0, nrow = q, ncol = kmax)
  if(any(output %in% "coefSamples")) BETA = list()
  if(any(output %in% "numFactors")) K = rep(NA, sp)
  if(any(output %in% "numFactors")) runtime = NULL
  if(any(output %in% "factMean")) FACTMEAN = matrix(0, nrow = n, ncol = kmax)
  if(any(output %in% "factSamples")) ETA = list()
  ind = 1
  
  
  #------start gibbs sampling-----#
  
  #cat("Start\n")
  
  t0 = proc.time()
  for (i in 1:nrun){
    
    
    # -- 1. Update eta -- #
    Lmsg = Lambda * ps
    Veta1 = diag(k) + t(Lmsg) %*% Lambda
    eigs = eigen(Veta1)
    if(all(eigs$values > 1e-6)) {
      Tmat = sqrt(eigs$values) * t(eigs$vectors)
    } else {
      Tmat = chol(Veta1)
    }
    R = qr.R(qr(Tmat))
    S = solve(R)
    Veta = S %*% t(S)                                               # Veta = inv(Veta1)
    Meta = Y %*% Lmsg %*% Veta                                      # n x k 
    eta = Meta + matrix(rnorm(n*k), nrow = n, ncol = k) %*% t(S)    # update eta in a block
    
    # -- 2. Update Sigma -- # 
    Ytil = Y - tcrossprod(eta,Lambda) 
    ps = rgamma(p, as + 0.5*n, bs+0.5*colSums(Ytil^2))
    Sigma = diag(1/ps)
    
    # -- 3. Update beta -- #
    pred = X%*%Beta
    logit = plogis(pred)
    
    Phi_L = matrix(1, nrow = p, ncol = k)
    logit_phi0 = logit[which(Phi==0)]
    which_zero = which(runif(length(logit_phi0))<
                         ((1-logit_phi0)/(1- logit_phi0 *p_constant)))
    Phi_L[ which(Phi==0)[which_zero] ] = 0
    
    Dt = matrix( pgdraw(1, pred), nrow = p, ncol = k )
    
    for(h in 1:k){
      Bh_1 = diag(rep(b_beta,q)^{-1})
      
      Qbe = t(X)%*%diag(Dt[,h])%*%X +Bh_1
      bbe = crossprod(X, (Phi_L[,h]-0.5)) 
      Lbe = t(chol(Qbe))
      # mean
      vbe = forwardsolve(Lbe, bbe)
      mbe = backsolve(t(Lbe), vbe)
      # var
      zbe = rnorm(q)
      ybe = backsolve(t(Lbe), zbe)
      Beta[,h] = t(ybe + mbe)
    }
    
    # -- 4. Update Lambda -- #
    for(j in 1:p) {
      
      etaj = t(t(eta)*rho*Phi[j,] ) 
      eta2 = crossprod(etaj)
      
      Qlam = Plam + ps[j]*eta2
      blam = ps[j]*crossprod(etaj,Y[,j])
      Llam = t(chol(Qlam))
      # mean
      vlam = forwardsolve(Llam,blam)
      mlam = backsolve(t(Llam),vlam)
      # var
      zlam = rnorm(k)
      ylam = backsolve(t(Llam), zlam)
      Lambda_star[j,] = t(ylam + mlam)
    }
    Lambda = t(t(Lambda_star)* sqrt(rho)) * Phi
    
    # -- 5.1 Update z -- #
    sdy = matrix(rep(sqrt(1/ps), n), n, p, byrow=T)
    phi_lam = Phi * Lambda_star
    index(phi_lam) = c("j","h")
    index(eta) = c("i", "h")
    eta_phi_lam = einstein(eta, phi_lam, drop = F)  # n x p x k
    
    mu = apply( aperm( aperm(eta_phi_lam, c(3,1,2))*rho, c(2,3,1)),
                c(1,2), sum)
    
    for(h in 1:k){
      
      mu_0 = mu - rho[h]*eta_phi_lam[,,h]
      mu_1 = mu_0 + eta_phi_lam[,,h]
      
      f0 = sum(dnorm(Y, mean= mu_0, sd=sdy, log=T))
      f1 = sum(dnorm(Y, mean= mu_1, sd=sdy, log=T))
      mf = max(c(f0,f1))
      f0 = f0 - mf
      f1 = f1 - mf
      
      prob_h = exp( c(rep(f0, h), rep(f1, k-h)) +log(w))
      
      if (sum(prob_h)==0){
        prob_h = c(rep(0,k-1), 1)
      } else{
        prob_h = prob_h/sum(prob_h)
      }
      z[h] = which(rmultinom(n=1, size=1, prob=prob_h)==1)
    }
    
    # rho
    rho = rep(1,k)
    rho[z <= seq(1,k)]=0
    
    
    # -- 5.2 Update precision matrix of Lambda_star -- #
    Plam = diag(rgamma(k, a_theta+0.5*p, b_theta+0.5*colSums(Lambda_star^2)))
    
    
    # -- 5.3 Update v and w -- #
    for(h in 1:(k-1)){
      v[h] = rbeta(1, shape1 = 1+sum(z==h), shape2 = alpha+sum(z>h))
    }
    v[k] = 1
    w = v*c(1,cumprod(1-v[-k]))
    
    #-- 6. Update Phi -- #
    pred = X%*%Beta 
    logit = plogis(pred)
    
    rho_lam = t(t(Lambda_star)*rho )
    index(rho_lam) = c("j","h")
    eta_rho_lam = einstein(eta, rho_lam, drop = F)  # n x p x k
    
    mu =  apply( aperm( aperm(eta_rho_lam, c(2,3,1))*as.vector(Phi), c(3,1,2)),
                 c(1,2), sum)
    
    for(h in 1:k){
      mu_0 = mu - t(t(eta_rho_lam[,,h])*Phi[,h])  
      mu_1 = mu_0 + eta_rho_lam[,,h]
      
      f0 = colSums(dnorm(Y, mean= mu_0, sd=sdy, log=T))
      f1 = colSums(dnorm(Y, mean= mu_1, sd=sdy, log=T))
      mf = max(c(f0,f1))
      f0 = f0 - mf
      f1 = f1 - mf
      
      lp_phi0 = f0 + log(1-logit[,h]*p_constant)
      lp_phi1 = f1 + log(logit[,h]*p_constant)
      
      sumlog = apply(cbind(lp_phi0, lp_phi1),1, logSumExp)
      
      Phi[,h] = round( runif(p) < exp(lp_phi1-sumlog) )
    }
    
    # -- save sampled values (after thinning) -- #
    if((i %% thin == 0) & (i > burn)) {
      if(k>=kmax){ Lambda_mean = Lambda[,1:kmax]; Phi_mean = Phi[,1:kmax]; eta_mean=eta[,1:kmax]
      Beta_mean = Beta[,1:kmax]
      } else{ Lambda_mean = cbind(Lambda, matrix(0, p, kmax-k))
      Phi_mean = cbind(Phi, matrix(0, p, kmax-k))
      eta_mean= cbind(eta, matrix(0, n, kmax-k))
      Beta_mean=cbind(Beta, matrix(0, q, kmax-k))
      }
      
      Omega = (Lambda %*% t(Lambda) + Sigma) * scaleMat
      if(any(output %in% "covMean")) COVMEAN = COVMEAN + Omega / sp
      if(any(output %in% "covSamples")) OMEGA[,,ind] = Omega
      if(any(output %in% "loadMean"))  LOADMEAN = LOADMEAN + Lambda_mean * sqrt(VY) / sp
      if(any(output %in% "loadAbsMean"))  LOADABSMEAN = LOADABSMEAN + abs(Lambda_mean) * sqrt(VY) / sp
      if(any(output %in% "loadSamples")) LAMBDA[[ind]] = Lambda * sqrt(VY)
      if(any(output %in% "locMean")) LOCMEAN = LOCMEAN + Phi_mean / sp
      if(any(output %in% "locSamples")) PHI[[ind]] = Phi
      if(any(output %in% "coefMean")) COEFMEAN = COEFMEAN + Beta_mean / sp
      if(any(output %in% "coefSamples")) BETA[[ind]] = Beta_mean
      if(any(output %in% "numFactors")) K[ind] = kstar
      if(any(output %in% "factMean")) FACTMEAN = FACTMEAN + eta_mean/sp
      if(any(output %in% "factSamples")) ETA[[ind]] = eta 
      
      ind = ind + 1
    }
    
    # --- Adaptation --- #
    
    
    if((uu[i] < prob[i])&(i> start_adapt)){
      
      #cat("Adaptation at iteration ", i," - rho: ", rho,"\n")
      active = which(z>seq(1,k))
      kstar = length(active)
      
      
      if (kstar< k-1){
        # set truncation to kstar and subset all variables, keeping only active columns
        k = kstar+1
        
        eta = cbind( eta[, active, drop = F], rnorm(n))
        vartheta_k = rgamma(1, a_theta, b_theta)
        Plam = diag( c(diag(Plam)[active, drop = F], vartheta_k))
        Lambda_star = cbind( Lambda_star[, active, drop = F], rnorm(p, 0, sd=sqrt(vartheta_k)))
        Phi = cbind( Phi[, active, drop = F], rbinom(p, size = 1, prob = p_constant))
        rho = c(rho[active, drop = F], 1)
        Lambda = cbind( Lambda[, active, drop = F], 
                        Lambda_star[,k]*sqrt(rho[k])*Phi[,k] )
        Beta = cbind( Beta[, active, drop = F], rnorm(q, 0, sd=sqrt(b_beta)))
        w = c(w[active, drop = F], 1-sum(w[active, drop = F]) )
        v = c(v[active, drop = F], 1)  # just to allocate memory
        z = c(z[active, drop=F], k)
        
      } else if (k<kmax) {
        # increase truncation by 1 and extend all variables, sampling from the prior/model
        k = k+1
        eta = cbind( eta, rnorm(n))
        vartheta_k = rgamma(1, a_theta, b_theta)
        Plam = diag(c(diag(Plam), vartheta_k) )
        Lambda_star = cbind( Lambda_star, rnorm(p, 0, sd=sqrt(vartheta_k)))
        Phi = cbind( Phi, rbinom(p, size = 1, prob = p_constant))
        rho = c(rho, 1)
        Lambda = cbind( Lambda, Lambda_star[,k]*sqrt(rho[k])*Phi[,k] )
        Beta = cbind( Beta, rnorm(q, 0, sd=sqrt(b_beta)))
        v[k-1] = rbeta(1, shape1 = 1, shape2 = alpha)
        v = c(v,1)
        w = v*c(1,cumprod(1-v[-k]))
        z = c(z, k)
      }
      
    }
    
    
    if((i %% 1000) == 0) {
      #cat("---", i,"\n", k, "\n" )
    }
  }
  
  runtime = proc.time()-t0
  out = lapply(output, function(x) {
    if(x == "covMean") return(COVMEAN)
    if(x == "covSamples") return(OMEGA)
    if(x == "loadMean") return(LOADMEAN)
    if(x == "loadAbsMean") return(LOADABSMEAN)
    if(x == "loadSamples") return(LAMBDA)
    if(x == "locMean") return(LOCMEAN)
    if(x == "locSamples") return(PHI)
    if(x == "coefMean") return(COEFMEAN)
    if(x == "coefSamples") return(BETA)
    if(x == "numFactors") return(K)
    if(x == "time") return(runtime[1])
    if(x == "factMean") return(FACTMEAN)
    if(x == "factSamples") return(ETA)
  })
  names(out) = output
  out[["model_prior"]] ="SIS"
  out[["data"]] =Y
  out[["covariates"]]=X
  out[["hyperparameters"]]=list(a_sigma=as, b_sigma=bs, alpha=alpha, a_theta=a_theta,
                                b_theta = b_theta, b_beta = b_beta) 
  return(out)
}




postprocessing_SIS=function(out_MCMC, likelihood="gaussian", 
                            model_prior = "SIS", frac_sampled=1, 
                            parallel=F, ncpu=ifelse(parallel,2,1) ){
  
  # check parallel
  parallel = ifelse(ncpu>1, T, F)
  
  # check model prior
  if(!(model_prior %in% c("SIS","CUSP"))){ 
    stop("argument 'model_prior' is not valid: it can be either SIS or CUSP")
  }
  
  # dimension of the data
  p = ncol(out_MCMC$data)
  n = nrow(out_MCMC$data)
  
  # number of iteration
  t = length(out_MCMC$numFactors)
  # finding final number of columns
  K = max(out_MCMC$numFactors)
  
  # sampling
  if((frac_sampled<=0)|(frac_sampled>1)){
    stop("'frac_sampled' not valid: it must be a number in (0,1]")
  }
  sampled = sample(1:t, ceiling(t*frac_sampled)) 
  
  
  lposterior_function = function(ind){
    hyperpar = out_MCMC$hyperparameters
    
    # log-likelihood
    if(likelihood=="gaussian"){
      ll = sum(dmvn(out_MCMC$data, mu=rep(0,p), Sigma = out_MCMC$covSamples[,,ind], log=T))
    }else{
      ## we compute the prob to have that order simulating from mvn(mu=0, sigma=Omega)
    }
    # complete Lambda with correct number of factors
    kl = ncol(out_MCMC$loadSamples[[ind]])
    if(kl>=K){
      Lambda = out_MCMC$loadSamples[[ind]][,1:K]
    }else{  Lambda = cbind(out_MCMC$loadSamples[[ind]], rep( rep(0,p), K-kl)) }
    
    
    # prior of Lambda 
    if(model_prior=="SIS"){
      # posterior of lambda, Sigma, beta
      if(kl>=K){
        Beta = out_MCMC$coefSamples[[ind]][,1:K]
      }else{  Beta = cbind(out_MCMC$coefSamples[[ind]], rep( rep(0,p), K-kl)) }
      prob = plogis(out_MCMC$covariates %*% Beta)
      p_phi1 =  2*prob*exp(1)*log(p)/p  
      
      # log prior of beta
      lp_beta = sum(dnorm(Beta, mean=0, sd=sqrt(hyperpar$b_beta), log = T))
      
      lp_col_lambda = rep(0, K)
      for(h in 1:K){
        w0 = which(Lambda[,h] == 0)
        w0l = length(w0)
        if(w0l==p){
          lp1 = log(p_gamma1[h])+sum(log(1-p_phi1[,h]))
          lp_col_lambda[h] = logSumExp(c( log(1-p_gamma1[h]),  lp1))
        }else if(w0l>0){
          lp_col_nonzero = LaplacesDemon::dmvt(x=Lambda[-w0,h], mu=rep(0,p-w0l),
                                               S=(hyperpar$b_theta/hyperpar$a_theta)*diag(p-w0l),
                                               df=2*hyperpar$a_theta, log=T) 
          lp_col_lambda[h] = log(p_gamma1[h])+sum(log(1-p_phi1[w0,h]))+
            sum(log(p_phi1[-w0,h]))+lp_col_nonzero
        }else{
          lp_col_nonzero = LaplacesDemon::dmvt(x=Lambda[,h], mu=rep(0,p),
                                               S=(hyperpar$b_theta/hyperpar$a_theta)*diag(p),
                                               df=2*hyperpar$a_theta, log=T) 
          lp_col_lambda[h] = log(p_gamma1[h])+sum(log(p_phi1[,h]))+lp_col_nonzero
        }
      }
      
      lp_lambda = sum(lp_col_lambda) 
      
    }else if(model_prior=="CUSP"){
      
      # log_prior beta
      lp_beta= 0
      
      lp_nonzero = log(p_gamma1) +  LaplacesDemon::dmvt(x=t(Lambda), mu=rep(0,p),
                                                        S=(hyperpar$b_theta/hyperpar$a_theta)*diag(p),
                                                        df=2*hyperpar$a_theta, log=T) 
      lp_zero = log(1-p_gamma1) + 
        colSums(dnorm(Lambda, mean = 0, sd = hyperpar$theta_inf^(1/2), log = T) )
      lp_lambda = sum(apply(cbind(lp_nonzero, lp_zero), 1, logSumExp))
    }
    
    # Sigma
    sigma = diag(out_MCMC$covSamples[,,ind])-diag(tcrossprod(out_MCMC$loadSamples[[ind]]))
    
    # prior of Sigma
    lp_sigma = sum( dgamma(1/sigma, shape=hyperpar$a_sigma, rate = hyperpar$b_sigma, log =T ) )
    
    # log posterior
    lpost = ll + lp_lambda + lp_sigma + lp_beta
    return(lpost)
  }
  
  if(parallel){
    require(snowfall)
    
    sfInit(parallel=TRUE, cpus=ncpu, type="SOCK")
    
    # library in cluster
    sfLibrary(gtools)
    sfLibrary(LaplacesDemon)
    
    # load posterior function and other variables
    sfExport(lposterior_function, structure_function, out_MCMC,
             p,n, t,K,likelihood, model_prior)
    lposterior = unlist(sfLapply(sampled, lposterior_function))
    sfStop() 
  }else{
    require(gtools)
    require(LaplacesDemon)
    lposterior = unlist(lapply(sampled, lposterior_function))
  }
  
  iteration_max = sampled[which.max(lposterior)]
  loadings_max = out_MCMC$loadSamples[[iteration_max]][,1:out_MCMC$numFactors[iteration_max]]
  postprocess_output = list( iteration_max = iteration_max, loadings_max = loadings_max)
  
  return(postprocess_output)
}

