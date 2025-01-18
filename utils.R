generate_row<-function(x,mu_0,sigma_0,shape,scale,fc,sample_size,diff_ratio){
  con_mu<-rnorm(1,mean = mu_0,sd = sigma_0)
  con_sigma<-rinvgamma(1,shape = shape,scale=scale)
  exe_mu<-con_mu+log(fc)
  probs<-runif(1)
  if(probs<=diff_ratio)
  {
    con_tmp_intensity<-rlnorm(sample_size/2,meanlog = con_mu,sdlog=con_sigma)
    exe_tmp_intensity<-rlnorm(sample_size/2,meanlog = exe_mu,sdlog=con_sigma)
    temp_data<-append(con_tmp_intensity,exe_tmp_intensity)
    temp_data<-append(temp_data,TRUE)
    return(temp_data)
  }
  else{
    temp_data<-rlnorm(sample_size,meanlog = con_mu,sdlog=con_sigma)
    temp_data<-append(temp_data,FALSE)
    return(temp_data)
  }
}

generate_matrix_simu<-function(mu_0,sigma_0,shape,scale,fc,samplesize,featuresize,mr_all,mr_random,mr_nrandom,diff_ratio){
  intensity_data<-matrix(NA,nrow = featuresize,ncol= samplesize)
  intensity_data<-as.data.frame(t(apply(intensity_data,1,generate_row,mu_0=mu_0,sigma_0=sigma_0,shape=shape,scale=scale,fc=fc,sample_size=samplesize,diff_ratio=diff_ratio)))
  complete_data = intensity_data
  labels = intensity_data[dim(intensity_data)[2]]
  intensity_matrix = as.matrix(intensity_data[1:samplesize])
  log_intensity_matrix = log(intensity_matrix)

  # threshold matrix
  qtnumber<-stats::quantile(log_intensity_matrix,na.rm = T,mr_all)

  # compare to determine MNAR
  threshold<-rnorm(length(log_intensity_matrix),mean = qtnumber,sd = 0.01)

  id_1<-which(log_intensity_matrix<threshold)
  MNARnumber<-round(featuresize*samplesize*mr_all*mr_nrandom)
  
  if(length(id_1)<=MNARnumber){
    print('Threshold MNAR pool less than sampling number!')
    id_2=id_1
  }
  else{
    id_2<-sample(id_1,MNARnumber)
  }

  intensity_matrix[id_2]<-NA

  # determine MCAR
  random_position<-sample(which(!is.na(intensity_matrix)),featuresize*samplesize*mr_all*mr_random)
  intensity_matrix[random_position]<-NA

  # final data with missing values
  intensity_data<-as.data.frame(intensity_matrix)
  colnames(intensity_data) = seq(1,samplesize)
  intensity_data[,"label"]<-labels

  # remove features of all NA in each group
  filter_flag = which(
    (apply(intensity_data,1,function(x) sum(!is.na(x[1:samplesize/2])))!=0)|(apply(intensity_data,1,function(x) sum(!is.na(x[samplesize/2+1:samplesize])))!=0)
  )
  intensity_data<-intensity_data[filter_flag,]
  complete_data<-complete_data[filter_flag,]
  
  print(dim(intensity_data))
  print(sum(is.na(intensity_data))/dim(intensity_data)[1]/(dim(intensity_data)[2]-1))
  
  return(list(complete_data, intensity_data))
}

NAcol=function(dataset){
  n=nrow(dataset)
  y=colSums(is.na(dataset))
  if (sum(y==n)==0){
    return(FALSE)
  } else{
    NAsamp=which(y==n)
    return(list(TRUE,NAsamp,dataset[,-NAsamp]))
  }
}

SampMin=function(x){
  for (i in 1:dim(x)[2])
    x[,i][is.na(x[,i])]=min(x[,i],na.rm=T)
  return(x)
}