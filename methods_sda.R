############################################################SDA#######################################
library(S4Vectors)
library(SummarizedExperiment)
library(trust)

SDA <- function(sumExp, VOI = NULL, ...){

    if (!is(sumExp, "SummarizedExperiment"))
        stop("Input must be an object of SummarizedExperiment class!")

    rawfeature <- assay(sumExp); rawcoldata <- colData(sumExp)
    newdata <- data_clean(rawfeature)
    newfeature <- newdata$feature
    newcoldata <- as(rawcoldata,"data.frame")
    feat.names <- newdata$feat.names

    rawresult <- lapply(seq_len(dim(newfeature)[1]), function (i) SDA.unit(
        featurevec=newfeature[i,], phenodata = newcoldata, VOI = VOI))
    results <- Reduce('comb', rawresult)
    qv_1part <- apply(results$X1pvalue, 2, qvalue::qvalue,pi0 = 1,...)
    qv_2part <- qvalue::qvalue(results$X2pvalue,pi0 = 1,...)
    nparams = dim(results$pointest)[2]/2
    df.results <- list(gamma = as.matrix(results$pointest[,1:nparams]),
                        beta = as.matrix(results$pointest[,-(1:nparams)]),
                        pv_gamma = as.matrix(results$X1pvalue[,1]),
                        pv_beta = as.matrix(results$X1pvalue[,2]),
                        qv_gamma = as.matrix(qv_1part[[1]]$qvalues),
                        qv_beta = as.matrix(qv_1part[[2]]$qvalues),
                        pv_2part = results$X2pvalue,
                        qv_2part = qv_2part$qvalues,
                        feat.names = feat.names)
    return(df.results)
}

aft_model <- function(data_feature, phenodata, bw = 1){
    y <- data_feature
    n <- length(y)
    dmatrix = model.matrix(~., data = as.data.frame(phenodata))[,-1]
    dmatrix = as.matrix(dmatrix)
    binit = rep(0, ncol(dmatrix))

    if (bw == 1){an <- 1.144*sd(log(y)-as.vector(dmatrix%*%as.matrix(binit)))*
        n^(-1/5)}
    if (bw == 2){an <- sd(log(y)-as.vector(dmatrix%*%as.matrix(binit)))*n^(-1/5)}
    if (bw == 3){an <- sd(log(y)-as.vector(dmatrix%*%as.matrix(binit)))*n^(-1/7)}
    if (bw == 4){an <- sd(log(y)-as.vector(dmatrix%*%as.matrix(binit)))*n^(-1/9)}
    if (bw == 5){an <- 4^(1/3)*
      min(sd(log(y)-as.vector(dmatrix%*%as.matrix(binit))),
            IQR(log(y)-as.vector(dmatrix%*%as.matrix(binit)))/1.34)*n^(-1/5)}
    if (bw == 6){an <- (8*sqrt(2)/3)^(1/5)*
        min(sd(log(y)-as.vector(dmatrix%*%as.matrix(binit))),
            IQR(log(y)-as.vector(dmatrix%*%as.matrix(binit)))/1.34)*n^(-1/5)}
    if (bw == 7){an <- 4^(1/3)*
        min(sd(log(y)-as.vector(dmatrix%*%as.matrix(binit))),
            IQR(log(y)-as.vector(dmatrix%*%as.matrix(binit)))/1.34)*n^(-1/3)}

    kern <- dnorm
    kern.1st <- function(x){-x*dnorm(x)}
    kern.2nd <- function(x){(x^2-1)*dnorm(x)}
    kern.cdf <- pnorm

    e_diff <- function(beta.iter,dmatrix,y){
        e.diff<- outer(-log(y)+as.vector(t(dmatrix%*%as.matrix(beta.iter))),
                       log(y)-as.vector(t(dmatrix%*%as.matrix(beta.iter))), '+')
        return(e.diff)
    }

    loglikf_1 <- function(beta.iter, dmatrix, y){
        e.diff<- e_diff(beta.iter,dmatrix,y)
        loglik_value<- -1*sum(log(y))+sum(log(rowSums(kern(e.diff/an))/n/an))
        return(loglik_value)
    }

    fbeta_1 <- function(beta.iter, dmatrix, y){
        e.diff<- e_diff(beta.iter,dmatrix,y)
        fbeta<- rep(NA, length(beta.iter))
        for (p in 1:length(beta.iter)) {
            x_vec <- dmatrix[,p]
            x.diff<- outer(x_vec,-x_vec,'+')
            fbeta[p]<- 1/an*sum(rowSums(kern.1st(e.diff/an)*x.diff)/
                                    rowSums(kern(e.diff/an)))
        }
        return(fbeta)
    }

    fbeta_dev_1 <- function(beta.iter,dmatrix,y){
        e.diff <- e_diff(beta.iter,dmatrix,y)
        hessian.m <- matrix(NA, nrow = length(beta.iter),
                            ncol = length(beta.iter))
        for (i in 1:length(beta.iter)) {
            x_i <- dmatrix[,i]
            x_i_diff<- outer(x_i,-x_i,'+')
            for (j in 1:i) {
                x_j <- dmatrix[,j]
                x_j_diff <- outer(x_j,-x_j,'+')
                hessian.m[j,i] <- hessian.m[i,j]<- 1/an/an*
                    sum((rowSums(kern.2nd(e.diff/an)*x_i_diff*x_j_diff)*
                             rowSums(kern(e.diff/an))-
                             rowSums(kern.1st(e.diff/an)*x_i_diff)*
                             rowSums(kern.1st(e.diff/an)*x_j_diff))/
                            (rowSums(kern(e.diff/an)))^2)
            }
        }
        return(hessian.m)
    }

    objfun_1 <- function(beta.iter){
        stopifnot(is.numeric(beta.iter))
        f <- loglikf_1(beta.iter,dmatrix,y)
        g <- fbeta_1(beta.iter,dmatrix,y)
        B <- as.matrix(fbeta_dev_1(beta.iter,dmatrix,y))
        list(value =f, gradient = g, hessian = B)
    }

  trust.results <- try(trust(objfun_1, binit, 1, 5, minimize = FALSE),
                       silent = TRUE)
  if (class(trust.results) == 'try-error') {
    return(list(pointest = rep(NA, ncol(dmatrix)),
                seest = rep(NA, ncol(dmatrix)),
                null.deviance = NA, residual.deviance = NA))
  }
  else {
    beta.est <- trust.results$argument
    se.est <- try(sqrt(diag(-1*solve(trust.results$hessian))), silent = TRUE)
    if (class(se.est) == 'try-error'){se.est = rep(NA, ncol(dmatrix));
    null.deviance = 0; residual.deviance = 0}
    else {
      se.est <- se.est
      null.deviance <- -2*loglikf_1(rep(0, ncol(dmatrix)), dmatrix, y)
      residual.deviance <- -2*loglikf_1(beta.est, dmatrix, y)
    }
    return(list(pointest = beta.est,
                seest = se.est,
                null.deviance = null.deviance,
                residual.deviance = residual.deviance))
  }
}


SDA_1cov <- function(featurevec, phenodata, VOI = NULL, bw = 1){
    
    dmatrix <- model.matrix(~., data = phenodata)
    dmatrixNoInt <- as.matrix(dmatrix[,-1])
    if(is.null(VOI)){
        data0 <- data.frame(featurevec, phenodata)
        ind_coef <- grepl(colnames(data0)[2], colnames(dmatrix))
    } else {
        data0 <- data.frame(featurevec, phenodata[,VOI])
        ind_coef <- grepl(VOI, colnames(dmatrix))
    }
    colnames(data0)[2] <- 'grouping'
    data_binary <- data0
    data_binary$featurevec[data_binary$featurevec>0] <- 1
    data_AFT <- data0[data0$featurevec>0,]

    non0_cnt <- aggregate((data_binary$featurevec==1),
                          by = list(data_binary$grouping), FUN = sum)$x
    group_size <- table(data0$grouping)

    if(any(non0_cnt == group_size)){
        coef_logit <- rep(NA, sum(ind_coef)); se_logit <- rep(NA, sum(ind_coef))
        groups = as.data.frame(phenodata[data0$featurevec>0,])
        aft_summary <- aft_model(data_AFT$featurevec, groups, bw = bw)
        coef_aft <- aft_summary$pointest[(ind_coef[-1])];

        logitNullDev <- 0; logitDev <- 0
        aftNullDev <- aft_summary$null.deviance
        aftDev <- aft_summary$residual.deviance
    }
    else if(any(non0_cnt<2) ||
            all(data_AFT$featurevec[1]==data_AFT$featurevec)){
        coef_aft <- rep(NA, sum(ind_coef)); se_aft <- rep(NA, sum(ind_coef))
        logit_reg <- glm(data_binary$featurevec ~ 1 + dmatrixNoInt,
                        family = 'binomial', control = list(maxit = 50))
        logit_summary <- summary(logit_reg)
        coef_logit <- logit_summary$coefficients[ind_coef,1]

        logitNullDev <- logit_summary$null.deviance
        logitDev <- logit_summary$deviance
        aftNullDev <- 0; aftDev <- 0
    }
    else{
        logit_reg <- glm(data_binary$featurevec ~ 1 + dmatrixNoInt,
                        family = 'binomial', control = list(maxit = 50))
        logit_summary <- summary(logit_reg)
        coef_logit <- logit_summary$coefficients[ind_coef,1]
        groups = as.data.frame(phenodata[data0$featurevec>0,])
        aft_summary <- aft_model(data_AFT$featurevec, groups, bw = bw)
        coef_aft <- aft_summary$pointest[(ind_coef[-1])]
        logitNullDev <- logit_summary$null.deviance
        logitDev <- logit_summary$deviance
        aftNullDev <- aft_summary$null.deviance
        aftDev <- aft_summary$residual.deviance
    }

    return(list(logitNullDev = logitNullDev, logitDev = logitDev,
                aftNullDev = aftNullDev, aftDev = aftDev,
                pointest = c(coef_logit, coef_aft)))

}

SDA_cont <- function(featurevec, phenodata, VOI = NULL, bw = 1){
    dmatrix <- model.matrix(~., data = phenodata)
    dmatrixNoInt <- as.matrix(dmatrix[,-1])
    data0 <- data.frame(featurevec, phenodata[,VOI])
    ind_coef <- grepl(VOI, colnames(dmatrix))

    colnames(data0)[2] <- 'grouping'
    data_binary <- data0
    data_binary$featurevec[data_binary$featurevec>0] <- 1
    data_AFT <- data0[data0$featurevec>0,]
    logit_reg <- glm(data_binary$featurevec ~ 1 + dmatrixNoInt,
                     family = 'binomial', control = list(maxit = 50))
    logit_summary <- summary(logit_reg)
    coef_logit <- logit_summary$coefficients[ind_coef,1]
    groups = as.data.frame(phenodata[data0$featurevec>0,])
    aft_summary <- aft_model(data_AFT$featurevec, groups, bw = bw)
    coef_aft <- aft_summary$pointest[(ind_coef[-1])]
    logitNullDev <- logit_summary$null.deviance
    logitDev <- logit_summary$deviance
    aftNullDev <- aft_summary$null.deviance
    aftDev <- aft_summary$residual.deviance
    return(list(logitNullDev = logitNullDev, logitDev = logitDev,
                aftNullDev = aftNullDev, aftDev = aftDev,
                pointest = c(coef_logit, coef_aft)))

}
SDADev <- function(featurevec, phenodata, bw = 1){
    dmatrix <- model.matrix(~., data = phenodata)
    dmatrixNoInt <- as.matrix(dmatrix[,-1])
    data0 <- data.frame(featurevec, phenodata)

    colnames(data0)[2] <- 'grouping'
    data_binary <- data0
    data_binary$featurevec[data_binary$featurevec>0] <- 1
    data_AFT <- data0[data0$featurevec>0,]
    logit_reg <- glm(data_binary$featurevec ~ 1 + dmatrixNoInt,
                     family = 'binomial', control = list(maxit = 50))
    logit_summary <- summary(logit_reg)
    coef_logit <- NA
    groups = as.data.frame(phenodata[data0$featurevec>0,])
    aft_summary <- aft_model(data_AFT$featurevec, groups, bw = bw)
    coef_aft <- NA
    logitNullDev <- NA
    logitDev <- logit_summary$deviance
    aftNullDev <- NA
    aftDev <- aft_summary$residual.deviance
    return(list(logitNullDev = logitNullDev, logitDev = logitDev,
                aftNullDev = aftNullDev, aftDev = aftDev,
                pointest = c(coef_logit, coef_aft)))

}
getX2pv <- function(logitNullDev, logitDev, aftNullDev, aftDev){

    diff.dev.logit <- logitNullDev - logitDev
    diff.dev.aft <- aftNullDev - aftDev
    if (diff.dev.logit == 0) {p_logit <- NA}
    else {p_logit <- 1 - pchisq(diff.dev.logit, 1)}
    if (diff.dev.aft == 0) {p_aft <- NA}
    else {p_aft <- 1 - pchisq(diff.dev.aft, 1)}

    diff.total <- diff.dev.logit + diff.dev.aft

    if ((diff.dev.logit == 0)|(diff.dev.aft == 0)) {
        X2pv <- 1 - pchisq(diff.total, 1)
    } else {
        X2pv <- 1 - pchisq(diff.total, 2)
    }

    return(list(X1pvalue = c(p_logit, p_aft),
                X2pvalue = X2pv)
    )
}
#---------------------- deal with more than one covariate ----------------------
SDA.unit <- function(featurevec, phenodata, VOI = NULL, bw = 1){

    if(dim(phenodata)[2] == 1 || is.null(VOI)) {
        results1 = SDA_1cov(featurevec = featurevec, phenodata = phenodata,
                                bw = bw)
        res_pv = getX2pv(results1$logitNullDev, results1$logitDev,
                          results1$aftNullDev, results1$aftDev)

    } else if (is.factor(phenodata[,VOI])){
        results1 = SDA_1cov(featurevec = featurevec, phenodata = phenodata,
                            VOI = VOI, bw = bw)
        results2 = SDADev(featurevec = featurevec,
                          phenodata = phenodata[,-which(colnames(phenodata)==VOI)],
                            bw = bw)
        res_pv = getX2pv(results2$logitDev, results1$logitDev,
                         results2$aftDev, results1$aftDev)

    } else {
        results1 = SDA_cont(featurevec = featurevec, phenodata = phenodata,
                            VOI = VOI, bw = bw)
        results2 = SDADev(featurevec = featurevec,
                          phenodata = phenodata[,-which(colnames(phenodata)==VOI)],
                          bw = bw)
        res_pv = getX2pv(results2$logitDev, results1$logitDev,
                         results2$aftDev, results1$aftDev)
    }
    return(list(pointest = as.numeric(results1$pointest),
                X1pvalue = res_pv$X1pvalue, X2pvalue = res_pv$X2pvalue))
}


######## read data (two ways)   ########
#-------- get data from seperate matrix ---------
createSEFromMatrix <- function(feature, colData) {

    subjectName = rownames(colData)
    cName = colnames(colData)

    if (is(colData, "data.frame")) {
        colData <- as(colData, "DataFrame")
    }
    if (is(colData, "matrix")) {
        colData <- as(colData,'DataFrame')
    }

    feature <- as.matrix(feature)
    if (ncol(feature)!=nrow(colData)) {
        stop("Feature data and column data do not match!")
    }


    result <- SummarizedExperiment(assays = SimpleList(counts=feature),
                                   colData = colData)
    return(result)
}

# ------- import data from csv files ----------
createSEFromCSV <- function(featurePath, colDataPath, rownames1 = 1,
                            rownames2 = 1, header1 = TRUE, header2 = TRUE){

    feature <- read.csv(featurePath, row.names = rownames1, header = header1,
                        check.names = FALSE)
    colData <- read.csv(colDataPath, row.names = rownames2, header = header2,
                        check.names = FALSE)
    
    result <- createSEFromMatrix(feature = feature, colData = colData)

    return(result)

}

comb <- function(list1, list2){
    res <- list()

    res$pointest <- rbind(list1$pointest, list2$pointest)
    res$X1pvalue <- rbind(list1$X1pvalue, list2$X1pvalue)
    res$X2pvalue <- rbind(list1$X2pvalue, list2$X2pvalue)

    return(res)
}

data_clean <- function(rawfeature){

    #ind <- which(rowSums(rawfeature>0)>=10)
    #feature <- rawfeature[ind, ]
    feature <- rawfeature
    feat.names <- rownames(feature)

    return(list(feature = feature,
                feat.names = feat.names)
            )
}

sda <-function(dataset){# 
  featureinfo<-as.matrix(dataset)
  featureinfo[which(is.na(featureinfo))]<-0
  rownames(featureinfo)<-as.character(seq(1,nrow(featureinfo)))
  colnames(featureinfo)<-as.character(seq(1,ncol(featureinfo)))
  groupinfo<-append(rep(0,ncol(featureinfo)/2),rep(1,ncol(featureinfo)/2))
  groupinfo<-data.frame(grouping=groupinfo)
  rownames(groupinfo) <- colnames(featureinfo)
  exampleSE <- createSEFromMatrix(feature = featureinfo, colData = groupinfo)
  results <- SDA(exampleSE)
  return(results)
}