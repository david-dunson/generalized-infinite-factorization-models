#---------------------------------------------------------------------------------------------#
# ---- SIS with Gaussian and missing data for covariate matrix ----
# structured factorization of a covariate matrix Y in the model y = b Y + eps
# Copyright: Lorenzo Schiavon 2021,
# Licensed under Creative Commons Attribution-NonCommercial 4.0
#---------------------------------------------------------------------------------------------#


# y response variable
# Y covariates
# X meta-covariates

Mcmc_SIS_actions = function(y, Y, X=NULL, as, bs, alpha, a_theta, b_theta, sigma_beta = 1,  
                    a_nu=1, b_nu = 2, sigma_b = 1, sigma_b0 = 0,
                    continuous_X = NULL ,  kinit = NULL, kmax = NULL,
                    adapp0=1, adapp1=5*10^(-4), start_adapt = 50, 
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
                                        "locSamples", "coefMean", "coefSamples","shrinkCoefMean",
                                        "shrinkCoefSamples", "varMean", "varSamples")}
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  which.avail = which(!is.na(y))      # which elements of response variable are available
  n.avail = length(which.avail)
  
  # scale X continuous
  if(sum(continuous_X)>0){
    X[,continuous_X] = scale(X[,continuous_X])
  }
  
  # scale Y
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
  prob = 1/exp(adapp0 + adapp1*seq(1,nrun))               
  uu = runif(nrun)
  
  
  # --- Initial values --- #
  p_constant =2*exp(1)*log(p)/p          # factor probability constant 
  
  ps = rgamma(p, as, bs )  
  pnu = rgamma(1,a_nu, b_nu)
  ps.star = c(ps, pnu)
  Sigma.star=diag(1/ps.star)                  # Sigma = diagonal residual covariance
  
  Lambda_star = matrix(rnorm(p*k), nrow = p, ncol = k)
  eta = matrix(rnorm(n*k), nrow = n, ncol = k)    # factor loadings & latent factors
  b = rnorm(k)
  b0 = rnorm(1)
  
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
  if(any(output %in% "coefMean")) COEFMEAN = rep(0, kmax+1)
  if(any(output %in% "coefSamples")) B = list()
  if(any(output %in% "varMean")) VARMEAN = 0
  if(any(output %in% "varSamples")) SIGMANU = rep(0,sp)
  if(any(output %in% "covMean")) COVMEAN = matrix(0, nrow = p, ncol = p)
  if(any(output %in% "covSamples")) OMEGA = array(dim = c(p, p, sp))
  if(any(output %in% "loadMean")) LOADMEAN = matrix(0, nrow = p, ncol = kmax)
  if(any(output %in% "loadAbsMean")) LOADABSMEAN = matrix(0, nrow = p, ncol = kmax)
  if(any(output %in% "loadSamples")) LAMBDA = list()
  if(any(output %in% "locMean")) LOCMEAN = matrix(0, nrow = p, ncol = kmax)
  if(any(output %in% "locSamples")) PHI = list()
  if(any(output %in% "shrinkCoefMean")) SHRINKCOEFMEAN = matrix(0, nrow = q, ncol = kmax)
  if(any(output %in% "shrinkCoefSamples")) BETA = list()
  if(any(output %in% "numFactors")) K = rep(NA, sp)
  if(any(output %in% "time")) runtime = NULL
  if(any(output %in% "factMean")) FACTMEAN = matrix(0, nrow = n, ncol = kmax)
  if(any(output %in% "factSamples")) ETA = list()
  ind = 1
  
  
  #------start gibbs sampling-----#
  
  #cat("Start\n")
  
  t0 = proc.time()
  for (i in 1:nrun){
    
    
    # -- 0. Update intercept -- # 
    
    b0 = rnorm(1, mean=pnu*sum(y-eta%*%b, na.rm = T)/(sigma_b0^(-2)+(n.avail)*pnu),
          sd = (sigma_b0^{-2}+n.avail*pnu)^(-0.5))
    
    
    # -- 1. Update eta -- #
    # y available 
    Lmsg1 = rbind(Lambda, b)* ps.star
    Veta1 = diag(k) + t(Lmsg1) %*% rbind(Lambda, b)
    eigs1 = eigen(Veta1)
    if(all(eigs1$values > 1e-6)) {
      Tmat1 = sqrt(eigs1$values) * t(eigs1$vectors)
    } else {
      Tmat1 = chol(Veta1)
    }
    R1 = qr.R(qr(Tmat1))
    S1 = solve(R1)
    Meta1 = cbind(Y, y-b0)[which.avail,] %*% Lmsg1 %*% S1 %*% t(S1)           # n.avail x k 
    eta[which.avail,] = Meta1 + matrix(rnorm(n.avail*k), nrow = n.avail, ncol = k) %*% t(S1)    # update eta in a block
    
    # with y missing if they exist
    if(n-n.avail>0){
      Lmsg2 = Lambda* ps
      Veta2 = diag(k) + t(Lmsg2) %*% Lambda
      eigs2 = eigen(Veta2)
      if(all(eigs2$values > 1e-6)) {
        Tmat2 = sqrt(eigs2$values) * t(eigs2$vectors)
      } else {
        Tmat2 = chol(Veta2)
      }
      R2 = qr.R(qr(Tmat2))
      S2 = solve(R2)
      Meta2 = Y[-which.avail,] %*% Lmsg2 %*% S2 %*% t(S2)           # n-n.avail x k 
      eta[-which.avail,] = Meta2 + matrix(rnorm((n-n.avail)*k), nrow = n-n.avail, ncol = k) %*% t(S2)    # update eta in a block
    }
 
    
    # -- 2. Update Sigma -- # 
    Ytil = Y - tcrossprod(eta,Lambda) 
    ps = rgamma(p, as + 0.5*n, bs+0.5*colSums(Ytil^2))
    pnu = rgamma(1, a_nu+ 0.5*n.avail, b_nu+ 0.5*sum((y- eta%*%b-b0)^2, na.rm = T))
    ps.star = c(ps,pnu)
    Sigma.star = diag(1/ps.star)
    
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
      Bh_1 = diag(rep(sigma_beta,q)^{-2})
      
      Qbe = t(X)%*%diag(Dt[,h])%*%X +Bh_1
      bbe = sigma_beta^(-2)*crossprod(X, (Phi_L[,h]-0.5)) 
      Lbe = t(chol(Qbe))
      # mean
      vbe = forwardsolve(Lbe, bbe)
      mbe = backsolve(t(Lbe), vbe)
      # var
      zbe = rnorm(q)
      ybe = backsolve(t(Lbe), zbe)
      Beta[,h] = t(ybe + mbe)
    }
    
    # -- 4.1 Update Lambda -- #
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
    
    
    # -- 4.2 Update b -- #
    eta2 = crossprod(eta)
    
    Qb = sigma_b^(-2) * base::diag(k) + pnu*eta2
    bb = pnu*crossprod(eta[which.avail,], y[which.avail])
    Lb = t(chol(Qb))
    # mean
    vb = forwardsolve(Lb,bb)
    mb = backsolve(t(Lb),vb)
    # var
    zb = rnorm(k)
    yb = backsolve(t(Lb), zb)
    b = drop(t(yb + mb))
    
    
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
      Beta_mean = Beta[,1:kmax]; b_mean = b[1:kmax]
      } else{ Lambda_mean = cbind(Lambda, matrix(0, p, kmax-k))
      Phi_mean = cbind(Phi, matrix(0, p, kmax-k))
      eta_mean= cbind(eta, matrix(0, n, kmax-k))
      Beta_mean=cbind(Beta, matrix(0, q, kmax-k))
      b_mean = c(b, rep(0,kmax-k))
      }
      
      Omega = (Lambda %*% t(Lambda) + diag(1/ps)) * scaleMat
      if(any(output %in% "coefMean")) COEFMEAN = COEFMEAN +c(b0,b_mean)/sp
      if(any(output %in% "coefSamples")) B[[ind]] = c(b0,b)
      if(any(output %in% "varMean")) VARMEAN = VARMEAN + 1/(pnu*sp)
      if(any(output %in% "varSamples")) SIGMANU[ind] = 1/pnu
      if(any(output %in% "covMean")) COVMEAN = COVMEAN + Omega / sp
      if(any(output %in% "covSamples")) OMEGA[,,ind] = Omega
      if(any(output %in% "loadMean"))  LOADMEAN = LOADMEAN + Lambda_mean * sqrt(VY) / sp
      if(any(output %in% "loadAbsMean"))  LOADABSMEAN = LOADABSMEAN + abs(Lambda_mean) * sqrt(VY) / sp
      if(any(output %in% "loadSamples")) LAMBDA[[ind]] = Lambda * sqrt(VY)
      if(any(output %in% "locMean")) LOCMEAN = LOCMEAN + Phi_mean / sp
      if(any(output %in% "locSamples")) PHI[[ind]] = Phi
      if(any(output %in% "shrinkCoefMean")) SHRINKCOEFMEAN = SHRINKCOEFMEAN + Beta_mean / sp
      if(any(output %in% "shrinkCoefSamples")) BETA[[ind]] = Beta_mean
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
        b = c(b[active, drop = F], 0)
        Lambda = cbind( Lambda[, active, drop = F], 
                        Lambda_star[,k]*sqrt(rho[k])*Phi[,k] )
        Beta = cbind( Beta[, active, drop = F], rnorm(q, 0, sd=sqrt(sigma_beta)))
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
        b = c(b,0)
        Lambda = cbind( Lambda, Lambda_star[,k]*sqrt(rho[k])*Phi[,k] )
        Beta = cbind( Beta, rnorm(q, 0, sd=sqrt(sigma_beta)))
        v[k-1] = rbeta(1, shape1 = 1, shape2 = alpha)
        v = c(v,1)
        w = v*c(1,cumprod(1-v[-k]))
        z = c(z, k)
      }
      
    }
    
    
    if((i %% 1000) == 0) {
      cat("---", i,"\n", k,"\n", sqrt(mean((y- eta%*%b-b0)^2, na.rm = T)),"\n")
    }
  }
  
  runtime = proc.time()-t0
  out = lapply(output, function(x) {
    if(x == "coefMean") return(COEFMEAN)
    if(x == "coefSamples") return(B)
    if(x == "varMean") return(VARMEAN)
    if(x == "varSamples") return(SIGMANU)
    if(x == "covMean") return(COVMEAN)
    if(x == "covSamples") return(OMEGA)
    if(x == "loadMean") return(LOADMEAN)
    if(x == "loadAbsMean") return(LOADABSMEAN)
    if(x == "loadSamples") return(LAMBDA)
    if(x == "locMean") return(LOCMEAN)
    if(x == "locSamples") return(PHI)
    if(x == "shrinkCoefMean") return(SHRINKCOEFMEAN)
    if(x == "shrinkCoefSamples") return(BETA)
    if(x == "numFactors") return(K)
    if(x == "time") return(runtime[1])
    if(x == "factMean") return(FACTMEAN)
    if(x == "factSamples") return(ETA)
  })
  names(out) = output
  out[["model_prior"]] ="SIS"
  out[["response"]] = y
  out[["covariates"]] =Y
  out[["metacovariates"]]=X
  out[["hyperparameters"]]=list(a_sigma=as, b_sigma=bs, a_nu = a_nu, b_nu= b_nu,
                                alpha=alpha, a_theta=a_theta,
                                b_theta = b_theta, sigma_beta = sigma_beta, sigma_b = sigma_b) 
  return(out)
}




# posterior mode
mode_SIS_actions=function(out_MCMC){
  
  
  # dimension of the data
  p = ncol(out_MCMC$covariates)
  n = nrow(out_MCMC$covariates)
  q = ncol(out_MCMC$metacovariates)
  
  # number of iteration
  t = length(out_MCMC$numFactors)
  # finding final number of columns
  K = max(out_MCMC$numFactors)
  
  # which y available
  which.avail = which(!is.na(out_MCMC$response))

  # posterior  
  lposterior_function = function(ind){
    hyperpar = out_MCMC$hyperparameters

    
    cov_lambdabeta =out_MCMC$loadSamples[[ind]]%*%out_MCMC$coefSamples[[ind]][-1]
    Var.tmp = cbind(rbind(out_MCMC$covSamples[,,ind], t(cov_lambdabeta)),
                c(cov_lambdabeta, crossprod(out_MCMC$coefSamples[[ind]][-1])+ out_MCMC$varSamples[ind] ) )
    
    Variance = Var.tmp
    Variance[p+1,p+1] = Variance[p+1,p+1] *(1+hyperpar$sigma_b * solve(Var.tmp)[p+1,p+1])
    
    # log-likelihood
    ll.avail = sum(dmvn(cbind(out_MCMC$covariates,out_MCMC$response)[which.avail,],
                  mu=rep(0,p+1), Sigma = Variance, log=T))
    if(length(out_MCMC$response)-length(which.avail)>0){
      ll = ll.avail + sum(dmvn( out_MCMC$covariates[-which.avail,],
                          mu=rep(0,p), Sigma = out_MCMC$covSamples[,,ind], log=T))
    }else{ll =ll.avail }
    

    # complete Lambda, beta and b with correct number of factors
    kl = ncol(out_MCMC$loadSamples[[ind]])
    if(kl>=K){
      Lambda = out_MCMC$loadSamples[[ind]][,1:K]
      b = out_MCMC$coefSamples[[ind]][2:(K+1)]
      Beta = out_MCMC$shrinkCoefSamples[[ind]][,1:K]
    }else{  
      Lambda = cbind(out_MCMC$loadSamples[[ind]], matrix( 0, p, K-kl))
      b = c(out_MCMC$coefSamples[[ind]][-1], rep(0, K-kl))
      Beta = cbind(out_MCMC$shrinkCoefSamples[[ind]], rep( 0,q, K-kl))
    }
    
    # posterior of [lambda,b], Sigma, beta
    p_gamma1 = (hyperpar$alpha/(1+hyperpar$alpha))^seq(1,K)
    prob = plogis(out_MCMC$metacovariates %*% Beta)
    p_phi1 =  2*prob*exp(1)*log(p)/p  
    
    # log prior of beta
    lp_beta = sum(dnorm(Beta, mean=0, sd=sqrt(hyperpar$sigma_beta), log = T))
  
    # prior of Lambda     
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
    
    # prior of b
    lp_b = sum(dnorm(b, mean=0, sd=sqrt(hyperpar$sigma_b), log = T))
      
    # prior of sigma_nu
    lp_sigmanu = dgamma(1/out_MCMC$varSample[ind], shape = hyperpar$a_nu, rate = hyperpar$b_nu, log=T)
    
    # Sigma
    sigma = diag(out_MCMC$covSamples[,,ind])-diag(tcrossprod(out_MCMC$loadSamples[[ind]]))
    
    # prior of Sigma
    lp_sigma = sum( dgamma(1/sigma, shape=hyperpar$a_sigma, rate = hyperpar$b_sigma, log =T ) )
    
    # log posterior
    lpost = ll + lp_lambda + lp_sigma + lp_beta + lp_b + lp_sigmanu
    return(lpost)
  }
  
  set = 1:t
  
  require(gtools)
  require(LaplacesDemon)
  lposterior = unlist(lapply(set, lposterior_function))
  
  print(summary(lposterior))
  iteration_max = which.max(lposterior)
  loadings_max = out_MCMC$loadSamples[[iteration_max]][,1:out_MCMC$numFactors[iteration_max]]
  coefficient_max = out_MCMC$coefSamples[[iteration_max]][1:(out_MCMC$numFactors[iteration_max]+1)]
  shrinkCoefficient_max = out_MCMC$shrinkCoefSamples[[iteration_max]][,1:out_MCMC$numFactors[iteration_max]]
  postprocess_output = list( iteration_max = iteration_max, loadings_max = loadings_max,
                             coefficient_max = coefficient_max, shrinkCoefficient_max = shrinkCoefficient_max)
  
  return(postprocess_output)
}




# Multivariate regression Gibbs 
Mcmc_Normal = function(y, X, a_nu, b_nu, sigma_b = 1,
                       nrun, burn=round(nrun/4), thin = 1, 
                       output = "all", my_seed = 6784, std=T){
  
  # set seed
  set.seed(my_seed)
  
  # dimensions
  n = dim(X)[1]
  p = dim(X)[2]

  if(any(output %in% "all")){output = c("time", "coefMean", "coefSamples", "varMean", "varSamples")}
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  # scale Y
  if(std){
    VX= apply(X, 2, var)                  # explicitly preserve scale
    X = scale(X)
  }else{
    VX= rep(1, p)
  }
  
  # --- Initial values --- #
  pnu = rgamma(1,a_nu, b_nu)

  b = rnorm(p+1, sd = sigma_b)
  
  # --- Allocate output object memory --- #
  if(any(output %in% "coefMean")) COEFMEAN = rep(0, p+1)
  if(any(output %in% "coefSamples")) B = matrix(0,p+1,sp)
  if(any(output %in% "varMean")) VARMEAN = 0
  if(any(output %in% "varSamples")) SIGMANU = rep(0,sp)
  ind = 1
  
  
  #------start gibbs sampling-----#
  
  #cat("Start\n")
  
  t0 = proc.time()
  for (i in 1:nrun){
    
    # -- 1. Update coefficients -- # 
    
    X.int = cbind(rep(1,n), X)

    varb_1 = diag(rep(sigma_b,p+1)^{-2})
    
    Qb = pnu*crossprod(X.int) + varb_1
    bb = pnu*crossprod(X.int, y) 
    Lb = t(chol(Qb))
    # mean
    vb = forwardsolve(Lb, bb)
    mb = backsolve(t(Lb), vb)
    # var
    zb = rnorm(p+1)
    yb = backsolve(t(Lb), zb)
    b = drop(t(yb + mb))

    # -- 2. Update Sigma -- # 
    ytil = y - X.int%*%b 
    pnu = rgamma(1, a_nu+ 0.5*n, b_nu+ 0.5*sum(ytil^2))
    
    # -- save sampled values (after thinning) -- #
    if((i %% thin == 0) & (i > burn)) {
      
      if(any(output %in% "coefMean")) COEFMEAN = COEFMEAN +b/sp
      if(any(output %in% "coefSamples")) B[,ind] = b
      if(any(output %in% "varMean")) VARMEAN = VARMEAN + 1/(pnu*sp)
      if(any(output %in% "varSamples")) SIGMANU[ind] = 1/pnu

      ind = ind + 1
    }
    
    if((i %% 1000) == 0) {
      cat("---", i,"\n")
    }
  
  }
    
  runtime = proc.time()-t0
  out = lapply(output, function(x) {
    if(x == "coefMean") return(COEFMEAN)
    if(x == "coefSamples") return(B)
    if(x == "varMean") return(VARMEAN)
    if(x == "varSamples") return(SIGMANU)
  })
  names(out) = output
  out[["model_prior"]] ="Normal-linear"
  out[["response"]] = y
  out[["covariates"]] =X
  out[["hyperparameters"]]=list(a_nu = a_nu, b_nu= b_nu,
                                sigma_b = sigma_b) 
  return(out)
  
}
