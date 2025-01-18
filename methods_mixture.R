Mixture <- function(data, covar){

  # Create binary variable for point mass vs continuous component
  # pm = 0 if in point mass
  pm <- ifelse(data > 0,1,0)
  
  cont <- data[data > 0]
  cont.cov <- covar[data > 0, ,drop=FALSE]
  cens <- data[data==0]
  cens.cov <- covar[data==0, , drop=FALSE]
  # Getting starting parameter values
  data[data==0] <- 0.5*min(cont)
  propTrunc <- pnorm(min(log(cont)),mean=mean(log(data)), sd=sd(log(data)))
  propAll <- length(cens)/length(data)
  # OneMinusTau = 1-tau which is the proportion of true zeros and tau is the proportion
  # from the continuous component
  OneMinusTau <- propAll-propTrunc
  # Can get OneMinusTau less than 1. Check for this and substitute one fifth of
  # the observed proportion missing (this is arbitrary)
  if (OneMinusTau < 0) OneMinusTau <- propAll/5
  Tau <- 1-OneMinusTau
  beta0 <- log(Tau/(1-Tau))
  
  
  # FIT FULL MODEL 
  linear <- lm(log(cont)~cont.cov)
  Params <- c(coefficients(linear), beta0, rep(0,ncol(covar)), summary(linear)$sigma)
  out <- try(optim(Params, mixtureLogL, method="L-BFGS-B",obs=cont, cen=cens, obs.covars=cont.cov,
                   cen.covars=cens.cov, lower=c(rep(-Inf,length(Params)-1),0.05), control=list(maxit=1000)),
             silent=TRUE)
  
  # Check for convergence error
  if (class(out)!="try-error"){
    if (out$convergence==1){
      warning("Algorithm did not converge for full likelihood.")
      full.logl <- NA
    } else {
      full.logl <- out$value
      cntl.beta <- out$par[3]
      cntl.mean <- out$par[1]
      cntl.tau <- exp(out$par[3])/(1+exp(out$par[3]))
      case.mean <- out$par[1]+out$par[2]
      case.beta <- out$par[3]+out$par[4]
      case.tau <- exp(case.beta)/(1+exp(case.beta))
      full.logl <- -1*out$value[1]
      std <- out$par[5]
    }
  } else {
    # If there is a boundary error then optimize without the lower bound specified
    out <- nlminb(Params, mixtureLogL, obs=cont, cen=cens, obs.covars=cont.cov,
                  cen.covars=cens.cov, lower=c(rep(-Inf,length(Params)-1),0.05))
    cntl.mean <- out$par[1]
    cntl.beta <- out$par[3]
    cntl.tau <- exp(out$par[3])/(1+exp(out$par[3]))
    case.mean <- out$par[1]+out$par[2]
    case.beta <- out$par[3]+out$par[4]
    case.tau <- exp(case.beta)/(1+exp(case.beta))
    full.logl <- -1*out$objective[1]
    std <- out$par[5]
  }
  
  
  
  
  # FIT NULL MODEL
  linear.n <- lm(log(cont)~1)
  ParamsN <- c(coefficients(linear.n), beta0, summary(linear.n)$sigma)
  outN <- try(optim(ParamsN, mixtureLogL, method="L-BFGS-B",obs=cont, cen=cens,
                    obs.covars=NULL,
                    cen.covars=NULL, lower=c(rep(-Inf,length(ParamsN)-1),0.05), control=list(maxit=1000)),
              silent=TRUE)
  if (class(outN)!="try-error"){
    if (outN$convergence==1){
      warning("Algorithm did not converge for null likelihood.")
      null.logl <- NA
    } else {
      null.logl <- -outN$value
      mean.n <- outN$par[1]
      tau.n <- exp(outN$par[2])/(1+exp(outN$par[2]))
      std.n<- outN$par[3]
    }
  } else {
    # If there is a boundary error then optimize with nlminb
    outN <- nlminb(ParamsN, mixtureLogL, obs=cont, cen=cens, obs.covars=NULL,
                   cen.covars=NULL, lower=c(rep(-Inf,length(ParamsN)-1),0.05))
    null.logl <- -outN$objective
    mean.n <- outN$par[1]
    tau.n <- exp(outN$par[2])/(1+exp(outN$par[2]))
    std.n<- outN$par[3]
  }
  
  
  #FIT MEAN MODEL
  Params.m <- c(coefficients(linear.n), beta0, rep(0,ncol(covar)), summary(linear)$sigma)
  out.m <- try(optim(Params.m, mixtureLogL.m, method="L-BFGS-B",obs=cont, cen=cens, obs.covars=cont.cov,
                     cen.covars=cens.cov, lower=c(rep(-Inf,length(Params)-1),0.05), control=list(maxit=1000)),
               silent=TRUE)
  # Check for convergence error
  if (class(out.m)!="try-error"){
    if (out.m$convergence==1){
      warning("Algorithm did not converge for full likelihood.")
      mean.logl <- NA
    } else {
      mean.logl <- out.m$value
      cntl.beta.m <- out.m$par[2]
      cntl.mean.m <- out.m$par[1]
      cntl.tau.m <- exp(out.m$par[2])/(1+exp(out.m$par[2]))
      case.mean.m <- out.m$par[1]
      case.beta.m <- out.m$par[2]+out.m$par[3]
      case.tau.m <- exp(case.beta.m)/(1+exp(case.beta.m))
      mean.logl <- -1*out.m$value[1]
      std.m <- out.m$par[4]
    }
  } else {
    # If there is a boundary error then optimize without the lower bound specified
    out.m <- nlminb(Params.m, mixtureLogL.m, obs=cont, cen=cens, obs.covars=cont.cov,
                    cen.covars=cens.cov, lower=c(rep(-Inf,length(Params)-1),0.05))
    cntl.beta.m <- out.m$par[2]
    cntl.mean.m <- out.m$par[1]
    cntl.tau.m <- exp(out.m$par[2])/(1+exp(out.m$par[2]))
    case.mean.m <- out.m$par[1]
    case.beta.m <- out.m$par[2]+out.m$par[3]
    case.tau.m <- exp(case.beta.m)/(1+exp(case.beta.m))
    mean.logl <- -1*out.m$objective[1]
    std.m <- out.m$par[4]
  }
  
  
  #FIT zero MODEL
  Params.p <- c(coefficients(linear), beta0, summary(linear)$sigma)
  out.p <- try(optim(Params.p, mixtureLogL.p, method="L-BFGS-B",obs=cont, cen=cens, obs.covars=cont.cov,
                     cen.covars=cens.cov, lower=c(rep(-Inf,length(Params)-1),0.05), control=list(maxit=1000)),
               silent=TRUE)
  # Check for convergence error
  if (class(out.p)!="try-error"){
    if (out.p$convergence==1){
      warning("Algorithm did not converge for full likelihood.")
      pzero.logl <- NA
    } else {
      pzero.logl <- out.p$value
      cntl.beta.p <- out.p$par[3]
      cntl.mean.p <- out.p$par[1]
      cntl.tau.p <- exp(out.p$par[3])/(1+exp(out.p$par[3]))
      case.mean.p <- out.p$par[1]+out.p$par[2]
      case.beta.p <- out.p$par[3]
      case.tau.p <- exp(case.beta.p)/(1+exp(case.beta.p))
      pzero.logl <- -1*out.p$value[1]
      std.p <- out.p$par[4]
    }
  } else {
    # If there is a boundary error then optimize without the lower bound specified
    out.p <- nlminb(Params.p, mixtureLogL.p, obs=cont, cen=cens, obs.covars=cont.cov,
                    cen.covars=cens.cov, lower=c(rep(-Inf,length(Params)-1),0.05))
    cntl.beta.p <- out.p$par[3]
    cntl.mean.p <- out.p$par[1]
    cntl.tau.p <- exp(out.p$par[3])/(1+exp(out.p$par[3]))
    case.mean.p <- out.p$par[1]+out.p$par[2]
    case.beta.p <- out.p$par[3]
    case.tau.p <- exp(case.beta.p)/(1+exp(case.beta.p))
    pzero.logl <- -1*out.p$objective[1]
    std.p <- out.p$par[4]
  }
  
  
  
  Pvalues <- NULL
  Pvalues.m <- NULL
  Pvalues.p <- NULL
  Stats <- NULL
  Stats.m <- NULL
  Stats.p <- NULL
  cov.names <- colnames(covar)
  if (is.null(colnames(covar))){
    cov.names <- as.character(seq(1,ncol(covar),1))
  }
  
  X2 <- -2*(null.logl-full.logl)
  X2.m <- -2*(mean.logl-full.logl)
  X2.p <- -2*(pzero.logl-full.logl)
  p.value <- 1-pchisq(X2, df=2)
  p.value.m <- 1-pchisq(X2.m,df=1)
  p.value.p <- 1-pchisq(X2.p,df=1)
  Pvalues <- c(Pvalues,p.value)
  Pvalues.m <- c(Pvalues.m,p.value.m)
  Pvalues.p <- c(Pvalues.p,p.value.p)
  Stats <- c(Stats,X2)
  Stats.m <- c(Stats.m,X2.m)
  Stats.p <- c(Stats.p,X2.p)
  ans <- list(Factor=cov.names, Statistic=Stats, Statistic.m=Stats.m, Statistic.p=Stats.p,  
              P.value=Pvalues, P.value.m=Pvalues.m, P.value.p=Pvalues.p, 
              FullLogl=full.logl, NullLogl=null.logl,
              MeanLogl=mean.logl, PzeroLogl=pzero.logl,
              Estimates=c(CntlMean=as.numeric(cntl.mean), CntlTau=cntl.tau,
                          CaseMean=as.numeric(case.mean), CaseTau=case.tau, STD=std),
              Estimates.m=c(CntlMean.m=as.numeric(cntl.mean.m), CntlTau.m=cntl.tau.m,
                            CaseMean.m=as.numeric(case.mean.m), CaseTau.m=case.tau.m, STD.m=std.m),
              Estimates.p=c(CntlMean.p=as.numeric(cntl.mean.p), CntlTau.p=cntl.tau.p,
                            CaseMean.p=as.numeric(case.mean.p), CaseTau.p=case.tau.p, STD.p=std.p),
              Estimates.n=c(CntlMean.n=as.numeric(mean.n),CntlTau.n=as.numeric(tau.n), STD.n=std.n)
  )
  return(ans)
}


################################################################################
# mixtureLogl is the log likelihood for the mixture model that is maximized
# params is a vector of the model parameters with the linear coefficients first
# followed by the logistic coefficients and last is the standard deviation of the
# the log normal distirbution
# This function allows evaluation of multiple covariates
# obs is a vector of the observed (non-zero values)
# cen is a vector of the censored or 0 values (all values are 0)
# obs.covars is a matrix n.obs X p of the covariates of the non-zero observations
# cen.covars is a n.cen X p matrix of the covariates of the zero observations
################################################################################
mixtureLogL <- function(params, obs, cen, obs.covars, cen.covars){
  d <- min(obs)
  Obs.Covars <- cbind(rep(1,length(obs)), obs.covars)
  Cen.Covars <- cbind(rep(1,length(cen)), cen.covars)
  n.param <- ncol(Obs.Covars)
  betas <- 1:n.param
  gammas <- c((n.param+1):(2*n.param))
  sd.loc <- length(params)
  mu.obs <- Obs.Covars%*%params[betas]
  mu.cen <- Cen.Covars%*%params[betas]
  tau.obs <- 1/(1+exp(-1*(Obs.Covars%*%params[gammas])))
  tau.cen <- 1/(1+exp(-1*(Cen.Covars%*%params[gammas])))
  logl.obs <- sum(log(tau.obs*(exp((-1*(log(obs)-
                                          mu.obs)^2)/(2*(params[sd.loc]^2)))/(obs*sqrt(2*pi)*params[sd.loc]))))
  logl.cen <- sum(log((1-tau.cen)+tau.cen*pnorm((log(d)-mu.cen)/params[sd.loc])))
  logl <- logl.obs+logl.cen
  return(-logl)
}


mixtureLogL.m <- function(params, obs, cen, obs.covars, cen.covars){
  d <- min(obs)
  Obs.Covars <- cbind(rep(1,length(obs)), obs.covars)
  Cen.Covars <- cbind(rep(1,length(cen)), cen.covars)
  n.param <- ncol(Obs.Covars)
  betas <- 1
  gammas <- c(2:(1+n.param))
  sd.loc <- length(params)
  mu.obs <- params[betas]
  mu.cen <- params[betas]
  tau.obs <- 1/(1+exp(-1*(Obs.Covars%*%params[gammas])))
  tau.cen <- 1/(1+exp(-1*(Cen.Covars%*%params[gammas])))
  #print(paste("muobs=",mu.obs))
  #print(paste("mucen=",mu.cen))
  #print(paste("pobs=",tau.obs))
  #print(paste("pcen=",tau.cen))
  #print(paste("sd=",params[sd.loc]))
  logl.obs <- sum(log(tau.obs*(exp((-1*(log(obs)-
                                          mu.obs)^2)/(2*(params[sd.loc]^2)))/(obs*sqrt(2*pi)*params[sd.loc]))))
  logl.cen <- sum(log((1-tau.cen)+tau.cen*pnorm((log(d)-mu.cen)/params[sd.loc])))
  logl <- logl.obs+logl.cen
  return(-logl)
}


mixtureLogL.p <- function(params, obs, cen, obs.covars, cen.covars){
  d <- min(obs)
  Obs.Covars <- cbind(rep(1,length(obs)), obs.covars)
  Cen.Covars <- cbind(rep(1,length(cen)), cen.covars)
  n.param <- ncol(Obs.Covars)
  betas <- 1:n.param
  gammas <- 1+n.param
  sd.loc <- length(params)
  mu.obs <- Obs.Covars%*%params[betas]
  mu.cen <- Cen.Covars%*%params[betas]
  tau.obs <- 1/(1+exp(-1*(params[gammas])))
  tau.cen <- 1/(1+exp(-1*(params[gammas])))
  logl.obs <- sum(log(tau.obs*(exp((-1*(log(obs)-
                                          mu.obs)^2)/(2*(params[sd.loc]^2)))/(obs*sqrt(2*pi)*params[sd.loc]))))
  logl.cen <- sum(log((1-tau.cen)+tau.cen*pnorm((log(d)-mu.cen)/params[sd.loc])))
  logl <- logl.obs+logl.cen
  return(-logl)
}


try_mixture = function(x, covar){
  trym = try(Mixture(x, covar),silent=TRUE)
  if ('try-error' %in% class(trym)){
    # error when data has no NA, return NA
    return(NA)
  } else {
  return(trym$P.value)
  }
}
