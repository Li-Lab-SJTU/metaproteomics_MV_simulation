t_test<-function(x){
  p_value<-try(t.test(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)])$p.value,silent = T)
  if ('try-error' %in% class(p_value)){
    return (NA)
  }
  else{
    return (p_value)
  }
}

wilcoxon_test<-function(x){
  p_value<-try(wilcox.test(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)],exact=F)$p.value,silent=T)
  if ('try-error' %in% class(p_value)){
    return (NA)
  }
  else{
    return (p_value)
  }
}

modt<-function(dataset){
  ncontrol<-ncol(dataset)/2
  ncase<-ncontrol
  design=model.matrix(~-1+factor(c(rep(1,ncontrol),rep(2,ncase))))
  colnames(design)=c('control','case')
  fit=lmFit(dataset,design)
  contrast.matrix=makeContrasts(case-control,levels=design)
  fit1=contrasts.fit(fit,contrast.matrix)
  fit2=eBayes(fit1)
  return(fit2$p.value)
}

twopart<-function(data,test){
  Index1 <- seq(1,length(data)/2)
  Group1 <- data[Index1]
  Group2 <- data[setdiff(seq(1,length(data)),Index1)]
  n1 <- length(Group1)
  n2 <- length(Group2)
  obs <- c(n1, n2)
  success <- c(sum(!is.na(Group1)), sum(!is.na(Group2)))
  pointmass <- obs-success
  uniq1 <- length(unique(Group1[which(!is.na(Group1))]))
  uniq2 <- length(unique(Group2[which(!is.na(Group2))]))
  T2<-0
  B2<-0
  if(uniq1<2 | uniq2<2){
    T2<-0
    if(sum(pointmass)==0 | sum(pointmass)==n1+n2){
      B2<-0
    }else{
      B2<-suppressWarnings(prop.test(pointmass, obs)$statistic)
    }
  }
  else if(sum(pointmass)==0){
    B2<-0
    if(test=="t.test"){
      T2 <- t.test(Group1,Group2)$statistic^2
    }
    if(test=="wilcoxon"){
      W <- wilcox.test(Group1,Group2, exact=FALSE)$statistic
      mu <- (n1*n2)/2
      sigma <- sqrt((n1*n2*(n1+n2+1))/12)
      T2 <- ((abs(W-mu)-0.5)/sigma)^2
    }
  }
  else{
    B2 <- suppressWarnings(prop.test(pointmass, obs)$statistic)
    Group1<-Group1[which(!is.na(Group1))]
    Group2<-Group2[which(!is.na(Group2))]
    if(test=="t.test"){
      T2 <- t.test(Group1,Group2)$statistic^2
    }
    if(test=="wilcoxon"){
      W <- wilcox.test(Group1,Group2, exact=FALSE)$statistic
      n1<-length(Group1)
      n2<-length(Group2)
      mu <- (n1*n2)/2
      sigma <- sqrt((n1*n2*(n1+n2+1))/12)
      T2 <- ((abs(W-mu)-0.5)/sigma)^2
    }
  }
  X2 <- B2+T2
  if ((T2==0)|(B2==0)) {
    X2pv <- 1-pchisq(X2,1)
  }
  else if(T2==0 & B2==0){
    X2pv<-NA
  }
  else {
    X2pv <- 1-pchisq(X2,2)
  }
  
  return(X2pv)
}

################################################################################
# AFT function is wrapper function for formatting and fitting accelerated
# failure time model with survreg function in survival package
################################################################################
try_AFT <- function(data, covar){
  # Make event indicator
  # Z=0 censored, Z=1 event observed
  Z <- ifelse(data > 0,1,0)
  # Put in minimum value for censored values
  d <- min(data[data > 0])
  data[data==0] <- d
  # Create survival object
  surv <- Surv(data,Z,type="left")
  censor <- try(survreg(surv~covar, dist="lognormal"),silent=TRUE)
  if (!'try-error' %in% class(censor)){
    ans <- summary(censor)
    return(ans$chi)
  } else {
    return(NA)
  }
}
