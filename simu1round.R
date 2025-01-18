simu1round = function(a){
  rownames(a)<-as.character(seq(1,nrow(a)))
  real_predict<-as.numeric(a[,"label"])
  data.na = a[,1:ncol(a)-1]
  rm(a)
  
  data.na.log<-log2(data.na)
  rownames(data.na.log)<-as.character(seq(1,nrow(data.na.log)))
  
  data.zr=data.na
  data.zr[is.na(data.zr)] = 0
  rownames(data.zr)<-as.character(seq(1,nrow(data.zr)))
  
  group.info = as.matrix(c(rep(0,samplesize/2),rep(1,samplesize/2)))         

  ## No imputation
  modt_time = proc.time()
  modt_result_raw = modt(data.na.log)
  modt_result<-p.adjust(modt_result_raw,"BH")
  modt_time = proc.time() - modt_time

  t_time = proc.time()
  t_result_raw = apply(data.na.log,MARGIN = 1,t_test)
  t_result<-p.adjust(t_result_raw,"BH")
  t_time = proc.time() - t_time

  wilco_time = proc.time()
  wilco_result_raw = apply(data.na,MARGIN = 1,wilcoxon_test)
  wilco_result<-p.adjust(wilco_result_raw,"BH")
  wilco_time = proc.time() - wilco_time

  ## Two part
  twot_time = proc.time()
  twot_result_raw = apply(data.na.log,MARGIN = 1,twopart,test="t.test")
  twot_result<-p.adjust(twot_result_raw,"BH")
  twot_time = proc.time() - twot_time

  twowilco_time = proc.time()
  twowilco_result_raw = apply(data.na,MARGIN = 1,twopart,test="wilcoxon")
  twowilco_result<-p.adjust(twowilco_result_raw,"BH")
  twowilco_time = proc.time() - twowilco_time

  ## sda
  sda_time = proc.time()
  sda_out <-sda(data.na)
  sda_result_raw=rep(NA,nrow(data.na))
  sda_result=rep(NA,nrow(data.na))
  if (!('try-error' %in% class(sda_out))){
    sda_result_raw[as.numeric(sda_out$feat.names)]=sda_out$pv_2part
    sda_result[as.numeric(sda_out$feat.names)]=sda_out$qv_2part
  }
  sda_time = proc.time() - sda_time
  print('sda ok')

  ## AFT
  AFT_time = proc.time()
  AFT.chi = apply(data.zr, MARGIN = 1, function(x) try_AFT(x, covar=group.info))
  AFT_result_raw = pchisq(AFT.chi,1,lower=F)
  AFT_result = p.adjust(AFT_result_raw,'BH')
  AFT_time = proc.time() - AFT_time
  print('aft ok')

  ## Mixture
  mixture_time = proc.time()
  try_mixture = function(x){
    trym = try(Mixture(x, covar=group.info),silent=T)
    if ('try-error' %in% class(trym)){
      # error when data has no NA, return NA
      return(NA)
    } else {
    return(trym$P.value)
    }
  }

  mixture_result_raw = apply(data.zr, MARGIN = 1, try_mixture)
  mixture_result = p.adjust(mixture_result_raw,'BH')
  mixture_time = proc.time() - mixture_time
  print('mixture ok')

  ## DASEV
  DASEV_time = proc.time()
  DASEV.simu <- DASEV(indata=data.zr, cov.matrix=cbind(rep(1, samplesize),group.info), 
                      test_cov=2, min.non0n=3, requiredn=10, requiredn2=30, 
                      DL_method= "Fixed Difference", DL_value=0.1, maxit_MLE=100, 
                      maxit=10000, test_model=c("Both"))
  DASEV_result_raw = rep(NA,nrow(data.zr))
  DASEV_result = rep(NA,nrow(data.zr))

  DASEV_result_raw[as.numeric(DASEV.simu$feature_names)] = DASEV.simu$pvalue_both
  DASEV_result[as.numeric(DASEV.simu$feature_names)] = DASEV.simu$pvalue_both
  DASEV_result = p.adjust(DASEV_result,'BH')
  DASEV_time = proc.time() - DASEV_time
  print('DASEV ok')

  ## Imputation
  colsamp=NAcol(data.na)
  BPCA=TRUE
  if (!colsamp[[1]]){
    ### No NA col ###
    bPCA_time = proc.time()
    pc=try(pca(data.na,method='bpca',nPcs=2,verbose=F),silent = T)
    if ('try-error' %in% class(pc)){
      # write.table(as.data.frame(pc[1]), paste(cd, "/simu_results/Error_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
      pcsp=strsplit(pc[1],'computationally')[[1]]
      if (length(pcsp)==2){
        ## mod_bPCA is a modified version of bPCA imputation method implemented in the package 'pcaMethods'. ##
        source(paste(cd,'/../mod_bpca/pca.R',sep=''))
        source(paste(cd,'/../mod_bpca/bpca.R',sep=''))
        source(paste(cd,'/../mod_bpca/BPCA_initmodel.R',sep=''))
        source(paste(cd,'/../mod_bpca/BPCA_dostep.R',sep=''))
        source(paste(cd,'/../mod_bpca/checkData.R',sep=''))
        source(paste(cd,'/../mod_bpca/prep.R',sep=''))
        source(paste(cd,'/../mod_bpca/repmat.R',sep=''))
        source(paste(cd,'/../mod_bpca/errorHierarchic.R',sep=''))
        source(paste(cd,'/../mod_bpca/derrorHierarchic.R',sep=''))
        source(paste(cd,'/../mod_bpca/AllClasses.R',sep=''))
        source(paste(cd,'/../mod_bpca/AllGenerics.R',sep=''))
        source(paste(cd,'/../mod_bpca/methods-pcaRes.R',sep=''))
        my_pc=try(mod_pca(data.na,method='mod_bpca',nPcs=2,verbose=F),silent=T)
        if ('try-error' %in% class(my_pc)){
          # write.table(as.data.frame(my_pc[1]), paste(cd, "/simu_results/Error_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, "_2.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
          # write.table(data.na, paste(cd, "/simu_results/DataError_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
          BPCA=FALSE
        } else{
          data.bPCA=completeObs(my_pc)
        }
        library(pcaMethods)
      } else{
        # write.table(data.na, paste(cd, "/simu_results/DataError_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
        BPCA=FALSE
      }
    } else{
      data.bPCA=completeObs(pc)
    }
    bPCA_time = proc.time() - bPCA_time

    RF_time = proc.time()
    data.RF = missForest(data.na, ntree = 100, maxiter = 10)[['ximp']]
    RF_time = proc.time() - RF_time

    QR_time = proc.time()
    data.QR = 2**(impute.QRILC(data.na.log)[[1]])
    QR_time = proc.time() - QR_time

    KNN_time = proc.time()
    data.KNN = impute.knn(as.matrix(data.na), k = 10, rowmax = 1, colmax = 1)[['data']]
    KNN_time = proc.time() - KNN_time

    SampMin_time = proc.time()
    data.SampMin=SampMin(data.na)
    SampMin_time = proc.time() - SampMin_time
  } else{
    ### NA col ###
    # colsamp[[3]]: data.na but remove all NA col
    bPCA_time = proc.time()
    data.Tsm=SampMin(colsamp[[3]])
    # global min  % before %  imputation
    ris <- integer(ncol(data.Tsm)+length(colsamp[[2]]))
    if (length(ris)!=samplesize) print('Error in bPCA imputation!')
    ris[colsamp[[2]]] <- ncol(data.Tsm)+1L
    ris[-colsamp[[2]]] <- seq_len(ncol(data.Tsm))
    
    pc=try(pca(colsamp[[3]],method='bpca',nPcs=2,verbose=F),silent = T)
    if ('try-error' %in% class(pc)){
      # write.table(as.data.frame(pc[1]), paste(cd, "/simu_results/Error_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
      pcsp=strsplit(pc[1],'computationally')[[1]]
      if (length(pcsp)==2){
        source(paste(cd,'/../mod_bpca/pca.R',sep=''))
        source(paste(cd,'/../mod_bpca/bpca.R',sep=''))
        source(paste(cd,'/../mod_bpca/BPCA_initmodel.R',sep=''))
        source(paste(cd,'/../mod_bpca/BPCA_dostep.R',sep=''))
        source(paste(cd,'/../mod_bpca/checkData.R',sep=''))
        source(paste(cd,'/../mod_bpca/prep.R',sep=''))
        source(paste(cd,'/../mod_bpca/repmat.R',sep=''))
        source(paste(cd,'/../mod_bpca/errorHierarchic.R',sep=''))
        source(paste(cd,'/../mod_bpca/derrorHierarchic.R',sep=''))
        source(paste(cd,'/../mod_bpca/AllClasses.R',sep=''))
        source(paste(cd,'/../mod_bpca/AllGenerics.R',sep=''))
        source(paste(cd,'/../mod_bpca/methods-pcaRes.R',sep=''))
        my_pc=try(mod_pca(colsamp[[3]],method='mod_bpca',nPcs=2,verbose=F),silent=T)
        if ('try-error' %in% class(my_pc)){
          # write.table(as.data.frame(my_pc[1]), paste(cd, "/simu_results/Error_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, "_2.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
          # write.table(data.na, paste(cd, "/simu_results/DataError_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
          BPCA=FALSE
        } else{
          data.temp=completeObs(my_pc)
          data.bPCA=cbind(data.temp,min(data.Tsm))[,ris]
        }
        library(pcaMethods)
      } else{
        # write.table(data.na, paste(cd, "/simu_results/DataError_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
        BPCA=FALSE
      }
    } else{
      data.temp=completeObs(pc)
      data.bPCA=cbind(data.temp,min(data.Tsm))[,ris]
    }
    bPCA_time = proc.time() - bPCA_time
    
    RF_time = proc.time()
    data.temp = missForest(colsamp[[3]], ntree = 100, maxiter = 10)[['ximp']]
    data.RF = cbind(data.temp,min(data.Tsm))[,ris]
    RF_time = proc.time() - RF_time

    QR_time = proc.time()
    data.temp = 2**(impute.QRILC(log2(colsamp[[3]]))[[1]])
    data.QR = cbind(data.temp,min(data.Tsm))[,ris]
    QR_time = proc.time() - QR_time

    KNN_time = proc.time()
    data.temp = impute.knn(as.matrix(colsamp[[3]]), k = 10, rowmax = 1, colmax = 1)[['data']]
    data.KNN = cbind(data.temp,min(data.Tsm))[,ris]
    KNN_time = proc.time() - KNN_time

    SampMin_time = proc.time()
    data.SampMin=cbind(data.Tsm,min(data.Tsm))[,ris]
    SampMin_time = proc.time() - SampMin_time
  }
  ## imputation done
  print('imputation ok')
  
  ## log transform
  if (BPCA==TRUE){
    data.bPCA.log = log2(data.bPCA)
  }
  data.SampMin.log = log2(data.SampMin)
  data.QR.log = log2(data.QR)
  data.RF.log = log2(data.RF)
  data.KNN.log = log2(data.KNN)
  
  ## imputation test
  if (BPCA==TRUE){
    modt_bPCA_result_raw<-modt(data.bPCA.log)
    modt_bPCA_result<-p.adjust(modt_bPCA_result_raw,"BH")

    t_bPCA_time = proc.time()
    t_bPCA_result_raw<-apply(data.bPCA.log,MARGIN = 1,t_test)
    t_bPCA_result<-p.adjust(t_bPCA_result_raw,"BH")
    t_bPCA_time = proc.time() - t_bPCA_time + bPCA_time

    wilco_bPCA_time = proc.time()
    wilco_bPCA_result_raw<-apply(data.bPCA,MARGIN = 1,wilcoxon_test)
    wilco_bPCA_result<-p.adjust(wilco_bPCA_result_raw,"BH")
    wilco_bPCA_time = proc.time() - wilco_bPCA_time + bPCA_time

  } else {

    modt_bPCA_result_raw<-rep(NA,nrow(data.na))
    modt_bPCA_result<-rep(NA,nrow(data.na))

    t_bPCA_time = proc.time()
    t_bPCA_result_raw<-rep(NA,nrow(data.na))
    t_bPCA_result<-rep(NA,nrow(data.na))
    t_bPCA_time = proc.time() - t_bPCA_time + bPCA_time

    wilco_bPCA_time = proc.time()
    wilco_bPCA_result_raw<-rep(NA,nrow(data.na))
    wilco_bPCA_result<-rep(NA,nrow(data.na))
    wilco_bPCA_time = proc.time() - wilco_bPCA_time + bPCA_time
  }

  modt_SampMin_result_raw<-modt(data.SampMin.log)
  modt_SampMin_result<-p.adjust(modt_SampMin_result_raw,"BH")
  
  modt_QR_result_raw = modt(data.QR.log)
  modt_QR_result<-p.adjust(modt_QR_result_raw,"BH")

  modt_RF_result_raw = modt(data.RF.log)
  modt_RF_result<-p.adjust(modt_RF_result_raw,"BH")

  modt_KNN_result_raw = modt(data.KNN.log)
  modt_KNN_result<-p.adjust(modt_KNN_result_raw,"BH")
  
  t_SampMin_time = proc.time()
  t_SampMin_result_raw = apply(data.SampMin.log,MARGIN = 1,t_test)
  t_SampMin_result<-p.adjust(t_SampMin_result_raw,"BH")
  t_SampMin_time = proc.time() - t_SampMin_time + SampMin_time

  t_QR_time = proc.time()
  t_QR_result_raw = apply(data.QR.log,MARGIN = 1,t_test)
  t_QR_result<-p.adjust(t_QR_result_raw,"BH")
  t_QR_time = proc.time() - t_QR_time + QR_time

  t_RF_time = proc.time()
  t_RF_result_raw = apply(data.RF.log,MARGIN = 1,t_test)
  t_RF_result<-p.adjust(t_RF_result_raw,"BH")
  t_RF_time = proc.time() - t_RF_time + RF_time

  t_KNN_time = proc.time()
  t_KNN_result_raw = apply(data.KNN.log,MARGIN = 1,t_test)
  t_KNN_result<-p.adjust(t_KNN_result_raw,"BH")
  t_KNN_time = proc.time() - t_KNN_time + KNN_time

  wilco_SampMin_time = proc.time()
  wilco_SampMin_result_raw = apply(data.SampMin,MARGIN = 1,wilcoxon_test)
  wilco_SampMin_result<-p.adjust(wilco_SampMin_result_raw,"BH")
  wilco_SampMin_time = proc.time() - wilco_SampMin_time + SampMin_time

  wilco_QR_time = proc.time()
  wilco_QR_result_raw = apply(data.QR,MARGIN = 1,wilcoxon_test)
  wilco_QR_result<-p.adjust(wilco_QR_result_raw,"BH")
  wilco_QR_time = proc.time() - wilco_QR_time + QR_time

  wilco_RF_time = proc.time()
  wilco_RF_result_raw = apply(data.RF,MARGIN = 1,wilcoxon_test)
  wilco_RF_result<-p.adjust(wilco_RF_result_raw,"BH")
  wilco_RF_time = proc.time() - wilco_RF_time + RF_time

  wilco_KNN_time = proc.time()
  wilco_KNN_result_raw = apply(data.KNN,MARGIN = 1,wilcoxon_test)
  wilco_KNN_result<-p.adjust(wilco_KNN_result_raw,"BH")
  wilco_KNN_time = proc.time() - wilco_KNN_time + KNN_time

  output.data <- data.frame(modt_result,
                            sda_result,
                            t_result,
                            wilco_result,
                            twot_result,
                            twowilco_result,
                            modt_bPCA_result,
                            modt_SampMin_result,
                            modt_QR_result,
                            modt_RF_result,
                            modt_KNN_result,
                            AFT_result,
                            mixture_result,
                            DASEV_result,
                            t_bPCA_result,
                            t_SampMin_result,
                            t_QR_result,
                            t_RF_result,
                            t_KNN_result,
                            wilco_bPCA_result,
                            wilco_SampMin_result,
                            wilco_QR_result,
                            wilco_RF_result,
                            wilco_KNN_result,
                            real_predict
                            )
    colnames(output.data) <- c("ModT",
                                'SDA', 
                                'T-test', 
                                'Wilcoxon', 
                                'twoT', 
                                'twoWilcox', 
                                'Modt_bPCA',
                                'Modt_SampMin',
                                'Modt_QR',
                                'Modt_RF',
                                'Modt_KNN',
                                'AFT',
                                'Mixture',
                                'DASEV',
                                'T_bPCA',
                                'T_SampMin',
                                'T_QR',
                                'T_RF',
                                'T_KNN',
                                'Wilcox_bPCA',
                                'Wilcox_SampMin',
                                'Wilcox_QR',
                                'Wilcox_RF',
                                'Wilcox_KNN',
                                "TrueDE")

    output.data.raw <- data.frame(modt_result_raw,
                            sda_result_raw,
                            t_result_raw,
                            wilco_result_raw,
                            twot_result_raw,
                            twowilco_result_raw,
                            modt_bPCA_result_raw,
                            modt_SampMin_result_raw,
                            modt_QR_result_raw,
                            modt_RF_result_raw,
                            modt_KNN_result_raw,
                            AFT_result_raw,
                            mixture_result_raw,
                            DASEV_result_raw,
                            t_bPCA_result_raw,
                            t_SampMin_result_raw,
                            t_QR_result_raw,
                            t_RF_result_raw,
                            t_KNN_result_raw,
                            wilco_bPCA_result_raw,
                            wilco_SampMin_result_raw,
                            wilco_QR_result_raw,
                            wilco_RF_result_raw,
                            wilco_KNN_result_raw,
                            real_predict
                            )
    colnames(output.data.raw) <- c('ModT_raw',
                            'SDA_raw', 
                            'T-test_raw', 
                            'Wilcoxon_raw', 
                            'twoT_raw', 
                            'twoWilcox_raw', 
                            'Modt_bPCA_raw',
                            'Modt_SampMin_raw',
                            'Modt_QR_raw',
                            'Modt_RF_raw',
                            'Modt_KNN_raw',
                            'AFT_raw',
                            'Mixture_raw',
                            'DASEV_raw',
                            'T_bPCA_raw',
                            'T_SampMin_raw',
                            'T_QR_raw',
                            'T_RF_raw',
                            'T_KNN_raw',
                            'Wilcox_bPCA_raw',
                            'Wilcox_SampMin_raw',
                            'Wilcox_QR_raw',
                            'Wilcox_RF_raw',
                            'Wilcox_KNN_raw',
                            'TrueDE')
      
      time.data <- data.frame(modt_time,
                            sda_time,
                            t_time,
                            wilco_time,
                            twot_time,
                            twowilco_time,
                            AFT_time,
                            mixture_time,
                            DASEV_time,
                            t_bPCA_time,
                            t_SampMin_time,
                            t_QR_time,
                            t_RF_time,
                            t_KNN_time,
                            wilco_bPCA_time,
                            wilco_SampMin_time,
                            wilco_QR_time,
                            wilco_RF_time,
                            wilco_KNN_time
                            )
    colnames(time.data) <- c("ModT",
                              'SDA', 
                              'T-test', 
                              'Wilcoxon', 
                              'twoT', 
                              'twoWilcox', 
                              'AFT',
                              'Mixture',
                              'DASEV',
                              'T_bPCA',
                              'T_SampMin',
                              'T_QR',
                              'T_RF',
                              'T_KNN',
                              'Wilcox_bPCA',
                              'Wilcox_SampMin',
                              'Wilcox_QR',
                              'Wilcox_RF',
                              'Wilcox_KNN')
    if (BPCA==TRUE){
    return(list(output.data,
      output.data.raw,
      data.bPCA,
      data.SampMin,
      data.QR,
      data.RF,
      data.KNN,
      time.data))
    } else {
    return(list(output.data,
      output.data.raw,
      NA,
      data.SampMin,
      data.QR,
      data.RF,
      data.KNN,
      time.data))
    }
}