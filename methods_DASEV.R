#Get prior distribution for sigma using MLE with group means seperated
Getprior <- function(params.prior, lod, covars, yvec) {
  #print(params.prior)
  n.Params <- ncol(covars)
  mu.all <- covars%*%params.prior[1:n.Params]
  sd0 <- exp(params.prior[n.Params+1])
  mu.nonzero <- mu.all[yvec>-Inf]
  nonzero <- yvec[yvec>-Inf]
  A <- pnorm(lod, mu.nonzero, sd0) # control the CDF at detection limit is less than 1
  if (sum(A>0.9999999999)>0){
    mle= -100000000000000
  }else{
    mle <- -sum(log((1 - A)*sqrt(2*pi)*sd0)+(nonzero-mu.nonzero)^2/(2*sd0^2))
  }
  return(-mle)
}


loglikestep1 <- function(Params, lsd, selectmodel, lod, s0, d0, nbeta, cov.matrix, test_cov, yvec){
  n.betas <- 1:nbeta #mean paras
  n.gammas <- (nbeta+1):length(Params)#zero paras
  betas<-Params[n.betas]
  gammas <-Params[n.gammas]
  
  if(selectmodel=="F"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix%*%gammas)))
    }
  }
  if(selectmodel=="N"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix[, -test_cov]%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas,length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix[, -test_cov]%*%gammas)))
    }
  }
  if(selectmodel=="M"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix[, -test_cov]%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix%*%gammas)))
    }
  }
  if(selectmodel=="P"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix[, -test_cov]%*%gammas)))
    }
  }
  mu.nonzero <- mu.all[yvec>-Inf]
  p.nonzero <- p.all[yvec>-Inf]
  mu.zero <- mu.all[yvec==-Inf]
  p.zero <- p.all[yvec==-Inf]
  nonzero <- yvec[yvec>-Inf] # non-zero obs
  zero <- yvec[yvec==-Inf] #zero obs
  sd <- exp(lsd)
  if (length(mu.zero)>0){
    A <- pnorm(lod, mu.zero, sd) # control the CDF at detection limit is less than 1 #A[A=="NaN"]<-1
    #  print(sum(A>0.5)>0)
    if (sum(A>0.99999999999999)>0){
      logl= -100000000000000
    }else{
      loglzero <- sum(log(p.zero + (1 - p.zero)*A))
      loglnonzero <- sum(log((1 - p.nonzero))) - sum((nonzero - mu.nonzero)^2/(2*sd^2)) - log(sqrt(2*pi)*sd)*length(nonzero)
      logl <- loglzero + loglnonzero #+ loglprior
    }
  }else{
    loglnonzero <- sum(log((1 - p.nonzero))) - sum((nonzero - mu.nonzero)^2/(2*sd^2)) - log(sqrt(2*pi)*sd)*length(nonzero)
    logl <- loglnonzero #+ loglprior
  }
  return <- (-logl)
}


loglikestep2 <- function(Params, betagamma, selectmodel, lod, s0, d0,
                         nbeta, cov.matrix, test_cov, yvec){
  n.betas <- 1:nbeta #mean paras
  n.gammas <- (nbeta+1):length(betagamma)#zero paras
  betas<-betagamma[n.betas]
  gammas <-betagamma[n.gammas]
  
  if(selectmodel=="F"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix%*%gammas)))
    }
  }
  if(selectmodel=="N"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix[, -test_cov]%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix[, -test_cov]%*%gammas)))
    }
  }
  if(selectmodel=="M"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix[, -test_cov]%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix%*%gammas)))
    }
  }
  if(selectmodel=="P"){
    if(length(betas)==1){
      mu.all <- rep(betas, length(yvec))
    }else{
      mu.all <- cov.matrix%*%betas
    }
    if(length(gammas)==1){
      p.all <- 1/(1+exp((rep(gammas, length(yvec)))))
    }else{
      p.all <- 1/(1+exp((cov.matrix[, -test_cov]%*%gammas)))
    }
  }
  
  mu.nonzero <- mu.all[yvec>-Inf]
  p.nonzero <- p.all[yvec>-Inf]
  mu.zero <- mu.all[yvec==-Inf]
  p.zero <- p.all[yvec==-Inf]
  nonzero <- yvec[yvec>-Inf] # non-zero obs
  zero <- yvec[yvec==-Inf] #zero obs
  sd <- exp(Params)
  
  if (length(mu.zero)>0){
    A <- pnorm(lod, mu.zero, sd) # control the CDF at detection limit is less than 1
    if (sum(A>0.999999999999)>0){
      logl= -100000000000000
    }else{
      loglzero <- sum(log(p.zero + (1 - p.zero)*A))
      loglnonzero <- sum(log((1 - p.nonzero))) -sum((nonzero - mu.nonzero)^2/(2*sd^2)) - log(sqrt(2*pi)*sd)*length(nonzero)
      loglprior <- (d0/2)*log(d0*s0^2/2) - (d0/2 +1)*log(sd^2) - log(gamma(d0/2)) - d0*s0^2/(2*sd^2)
      logl <- loglzero + loglnonzero + loglprior
    }
  }else{
    loglnonzero <- sum(log((1 - p.nonzero))) -sum((nonzero - mu.nonzero)^2/(2*sd^2)) - log(sqrt(2*pi)*sd)*length(nonzero)
    loglprior <- (d0/2)*log(d0*s0^2/2) - (d0/2 +1)*log(sd^2) - log(gamma(d0/2)) - d0*s0^2/(2*sd^2)
    logl <- loglnonzero + loglprior
  }
  
  return <- (-logl)
}


DASEV <- function(indata, cov.matrix, test_cov, test_cov_conti=NULL, min.non0n=3, requiredn=10,
                  requiredn2=30, DL_method= "Fixed Difference", DL_value=0.1,
                  maxit_MLE=100, maxit=10000, test_model=c("Both", "Mean", "Pzero")){
  if (!(DL_method == "Fixed Difference"|DL_method == "Fixed Rate"|DL_method ==
        "Fixed Value")){
    stop("Please enter the correct method to calculate the detect limit.")
  }
  
  if (is.na(DL_value)){
    stop("Please enter the correct value to calculate the detect limit.")
  }
  
  if (DL_method == "Fixed Rate" & DL_value < 1) {
    message("You are using the fixed rate method to calculate the detect limit,
            your input DL_value is less than 1, if this is not correct, please
            stop the program and check your input.")
  }
  
  if(min.non0n<3){
    stop("The minimum number of nonzero observation should be 3.")
  }
  
  #get data ready based on feature inclusion criteria
  n.ParamsF <- ncol(as.matrix(cov.matrix))
  n.ParamsN <- ncol(as.matrix(cov.matrix[, -test_cov]))
  ldata <- as.matrix(indata)
  idint<- rowSums(ldata> 0) >= min.non0n
  ldata <- ldata[idint,]
  
  check_cov <- length(test_cov) - length(test_cov_conti)
  if (length(test_cov_conti)>0){
    test_cov_dummy <- test_cov[-test_cov_conti]
  }else{
    test_cov_dummy <- test_cov
  }
  if (check_cov > 0){
    if(length(test_cov_dummy) > 1){
      test.cov <- cov.matrix[, test_cov_dummy]
      fslc <- function(test.cov){
        f.slc <- function(ldata){
          ind<- (!all(ldata[test.cov==0]!=0)
                 &!all(ldata[test.cov==0]==0)
                 &!all(ldata[test.cov==1]!=0)
                 &!all(ldata[test.cov==1]==0))
        }
        indall <-apply(ldata, 1, f.slc)
      }
      allind <- apply(test.cov, 2, fslc)
      allindsub <- rowSums(allind==TRUE)==ncol(allind)
    }else {
      f.slc <- function(ldata){
        ind<- (!all(ldata[test.cov==0]!=0)
               &!all(ldata[test.cov==0]==0)
               &!all(ldata[test.cov==1]!=0)
               &!all(ldata[test.cov==1]==0))
      }
      test.cov <- cov.matrix[, test_cov]
      allindsub <-apply(ldata, 1, f.slc)
    }
    ldata <- ldata[allindsub,]
  }

  ldata <- log(ldata)
  feature.names <- rownames(ldata)
  
  varprior <- c()
  for (i in 1:nrow(ldata)) {
    yvec <- ldata[i,]
    n <- (yvec > -Inf)
    covar.nonzero <- cov.matrix[n,]
    estbeta <- lm(yvec[n]~0+covar.nonzero)
    betasF <- coefficients(estbeta)
    betasF[is.na(betasF)] <- 0
    lsd <- log(sd(yvec[n]))
    if (DL_method == "Fixed Difference") {
      lod<-min(yvec[n])-DL_value
    } else if (DL_method == "Fixed Rate") {
      lod<-min(yvec[n])/DL_value
    } else if (DL_method == "Fixed Value") {
      lod <- DL_value
    } else {
      stop("Please enter the correct method to calculate the detect limit.")
    }
    mleest <- optim(par= c(betasF,lsd), fn=Getprior, yvec=yvec, lod=lod, covars=cov.matrix,
                    method="BFGS", control = list(maxit = maxit))
    
    if(mleest$convergence!=0)
      cat("\nGet prior distribution Convergence problems for feature ", i)
    varprior <- c(varprior, (exp(mleest$par[length(mleest$par)]))^2)
  }
  
  idint1<- rowSums(ldata> 0) >= requiredn
  varprior.sub <- varprior[idint1]
  count <- rowSums(ldata> -Inf)
  if(sum(idint1) < requiredn2){
    sort.varprior <- varprior[order(-count)]
    varprior.sub <- sort.varprior[1:requiredn2]
  }
  p.mean1 <- mean(varprior.sub)
  p.var1 <- var(varprior.sub)
  d0 <- 2*p.mean1^2/p.var1 + 4
  s0 <- sqrt(p.mean1*(d0-2)/d0)
  
  
  
  
  ################################################################################
  # Getting post statistics
  Para.Full <- c()
  MLE.Full <- c()
  MLE.Null <- c()
  MLE.Mean <- c()
  MLE.Pzero <- c()
  
  Pvalue <- c() #pvalue for Full model vs null model
  Pmean <-c() #pvalue for Full model vs same mean
  Pzero <-c() #pvalue for Full model vs same zero proportion
  
  DL <- c()
  for (i in 1:nrow(ldata)) {
    #print(i)
    yvec <- ldata[i,]
    n <- (yvec>-Inf)
    covar.nonzero <- cov.matrix[n,]
    if (sum(covar.nonzero[, 1] == covar.nonzero[, 2]) == nrow(covar.nonzero)){
      covar.nonzero[, 1]=1-covar.nonzero[, 1]
    }
    estbeta <- lm(yvec[n]~0+covar.nonzero)
    betasF <- coefficients(estbeta)
    betasF[is.na(betasF)] <- 0
    
    if (DL_method == "Fixed Difference") {
      lod<-min(yvec[n])-DL_value
    } else if (DL_method == "Fixed Rate") {
      lod<-min(yvec[n])/DL_value
    } else if (DL_method == "Fixed Value") {
      lod <- DL_value
    } else {
      stop("Please enter the correct method to calculate the detect limit.")
    }
    DL <- c(DL, lod)
    lsd0 <- log(sqrt(varprior[i]))# initial value for variance parameters
    mu <- mean(yvec[n])
    p <- 1 - length(yvec[n])/length(yvec)/(1-pnorm(lod, mu, exp(lsd0)))#true zero proportion starting value
    if (p>0){
      gamma0 <- log((1-p)/p)
    }else{
      gamma0 <- 1e+200
    }
    gammasF <-  c(gamma0, rep(0, ncol(cov.matrix)-1))
    
    covar.nonzero <- cov.matrix[n,]
    estbetaN <- lm(yvec[n]~0+covar.nonzero[, -test_cov])
    betasN <- coefficients(estbetaN)
    betasN[is.na(betasN)] <- 0
    gammasN <- gammasF[-test_cov]
    
    xF <- c(betasF, gammasF)#Full model params
    xN <- c(betasN, gammasN)#Null model params
    xM <- c(betasN, gammasF)#Params for model with different zero percentages but same means
    xP <- c(betasF, gammasN)#Params for model with different means but same zero percentages
    
    
    #Full Model
    nbeta <- length(betasF)
    Diff <- 1
    time <- 0
    lsd <- lsd0
    xFP <- xF
    while (Diff > 0.01 & time < maxit_MLE){
      time <- time +1
      outF <- optim(xFP, fn=loglikestep1, method = "BFGS", lsd=lsd, selectmodel="F",
                    lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                    cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
      xFP <- outF$par
      if(outF$convergence!=0)
        cat("\nStep 1 Get post statistics full model Convergence problems for feature ",i)
      out <- optim(lsd, fn=loglikestep2, selectmodel="F", method = "BFGS", betagamma=xFP,
                   lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                   cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
      Diff <- abs(out$par-lsd)
      lsd <- out$par
      if(out$convergence!=0)
        cat("\nStep 2 Get post statistics full model Convergence problems for feature ",i)
    }
    
    Para.Full <- rbind(Para.Full, c(outF$par, exp(out$par)))
    MLE.Full <- c(MLE.Full, outF$value)
    
    #Null model
    if ("Both" %in% test_model){
      nbeta <- length(betasN)
      Diff <- 1
      time <- 0
      lsd <- lsd0
      xNP <- xN
      while (Diff > 0.01 & time < maxit_MLE){
        time <- time +1
        outN <- optim(xNP, fn=loglikestep1, method = "BFGS", lsd=lsd, selectmodel="N",
                      lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                      cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        xNP <- outN$par
        if(outN$convergence!=0)
          cat("\nStep 1 Get post statistics null model Convergence problems for feature ",i)
        out <- optim(lsd, fn=loglikestep2, selectmodel="N", method = "BFGS", betagamma=xNP,
                     lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                     cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        Diff <- abs(out$par-lsd)
        lsd <- out$par
        if(out$convergence!=0)
          cat("\nStep 2 Get post statistics null model Convergence problems for feature ",i)
        
      }
      MLE.Null<-c(MLE.Null, outN$value)
    }else{
      MLE.Null <- c()
    }
    
    
    #Test mean
    if ("Mean" %in% test_model){
      nbeta <- length(betasN)
      Diff <- 1
      time <- 0
      lsd <- lsd0
      xMP <- xM
      while (Diff > 0.01 & time < maxit_MLE){
        time <- time +1
        outM <- optim(xMP, fn=loglikestep1, method = "BFGS", lsd=lsd, selectmodel="M",
                      lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                      cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        xMP <- outM$par
        if(outM$convergence!=0)
          cat("\nStep 1 Get post statistics mean model Convergence problems for feature ",i)
        out <- optim(lsd, fn=loglikestep2, selectmodel="M", method = "BFGS", betagamma=xMP,
                     lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                     cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        Diff <- abs(out$par-lsd)
        lsd <- out$par
        if(out$convergence!=0)
          cat("\nStep 2 Get post statistics mean model Convergence problems for feature ",i)
        
      }
      MLE.Mean <-c(MLE.Mean, outM$value)
    }else{
      MLE.Mean <- c()
    }
    
    
    #Test zero proportion
    if ("Pzero" %in% test_model){
      nbeta <- length(betasF)
      Diff <- 1
      time <- 0
      lsd <- lsd0
      xPP <- xP
      while (Diff > 0.01 & time < maxit_MLE){
        time <- time +1
        outP <- optim(xPP, fn=loglikestep1, method = "BFGS", lsd=lsd, selectmodel="P",
                      lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                      cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        xPP <- outP$par
        if(outP$convergence!=0)
          cat("\nStep 1 Get post statistics mean model Convergence problems for feature ",i)
        out <- optim(lsd, fn=loglikestep2, selectmodel="P", method = "BFGS", betagamma=xPP,
                     lod=lod, d0=d0, s0=s0, nbeta=nbeta, test_cov=test_cov,
                     cov.matrix=cov.matrix, yvec=yvec, control = list(maxit=maxit))
        Diff <- abs(out$par-lsd)
        lsd <- out$par
        if(out$convergence!=0)
          cat("\nStep 2 Get post statistics mean model Convergence problems for feature ",i)
        
      }
      MLE.Pzero <-c(MLE.Pzero, outP$value)
    }else{
      MLE.Pzero <- c()
    }
  }
  n.Params = ncol(cov.matrix)
  for (k in 1:n.Params){
    colnames(Para.Full)[k] <- paste("beta.", k, sep="")
  }
  for (k in (n.Params+1):(2*n.Params)){
    colnames(Para.Full)[k] <- paste("gamma.", k-n.Params, sep="")
  }
  colnames(Para.Full)[2*n.Params+1] <- "sd"
  
  if ("Both" %in% test_model){
    Pvalue <- c(Pvalue, 1-pchisq(2*(MLE.Null-MLE.Full), 2))
  }else{
    Pvalue <- c()
  }#Null vs Full
  if ("Mean" %in% test_model){
    Pmean <- c(Pmean, 1-pchisq(2*(MLE.Mean-MLE.Full), 1))
  }else{
    Pmean <- c()
  }#Mean vs Full
  if ("Pzero" %in% test_model){
    Pzero <- c(Pzero, 1-pchisq(2*(MLE.Pzero-MLE.Full), 1))
  }else{
    Pzero <- c()
  }#Pzero vs Full
  
  result <- list(feature_names =feature.names,
                 pvalue_both=Pvalue,
                 pvalue_mean=Pmean,
                 pvalue_zero=Pzero,
                 DL=DL,
                 estimates=Para.Full)
  
  return(result)
}
