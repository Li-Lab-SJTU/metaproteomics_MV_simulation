rm(list = ls())
library(readxl)
library(actuar)
library(impute)
library(limma)
library(log4r)
library(missForest)
library(imputeLCMD)
library(survival)
library(fitdistrplus)
library(reshape2)
library(stringr)

cd = $your_working_dir
setwd(cd)

## import functions for generating data and statistical test
source("../utils.R")
source("../methods.R")
source("../methods_sda.R")
source("../methods_DASEV.R")
source("../methods_mixture.R")
source("../simu1round.R")

## logging file
create.logger("./simu_log.txt", level=2)
appender = file_appender("./simu_log.txt")

## load data
protein_data<-read.table(file = "../datasets/DIA-PASEF/Meta_DIA_14D-Pre_3Conditions_proteome.txt",sep="\t",header=T)
protein_data_need<-data.frame(protein_data)
protein_data<-NA
rownames(protein_data_need)<-protein_data_need[,1]
protein_data_need = protein_data_need[,-1]
rownum<-nrow(protein_data_need)
colnum<-ncol(protein_data_need)
protein_data_need_group1 = which(str_detect(colnames(protein_data_need),"14D"))


# normalization
protein_data_need_samplesum = colSums(protein_data_need,na.rm = T)
protein_data_need_samplesum_median = median(protein_data_need_samplesum)
norm_ratio = protein_data_need_samplesum_median / protein_data_need_samplesum
protein_data_need = sweep(protein_data_need, 2, norm_ratio, FUN='*')

# filter NA>90%
meansd_data<-protein_data_need[which(apply(protein_data_need,1,function(x) sum(is.na(x)))<=0.9*dim(protein_data_need)[2]),]
missing_ratio_filtered = sum(is.na.data.frame(meansd_data))/(dim(meansd_data)[1]*dim(meansd_data)[2])
meansd_data<-log(meansd_data+1)

### parameters for dataset distribution
data_mean<-rowMeans(meansd_data,na.rm = T)
data_sd<-apply(meansd_data, 1, stats::sd,na.rm=T)

fit.norm = fitdist(data_mean,'norm','mge')
mu_0=fit.norm[['estimate']][1]
sigma_0=fit.norm[['estimate']][2]

fit.invgamma = fitdist(data_sd,"invgamma","mge")
shape=fit.invgamma[["estimate"]][1]
scale=fit.invgamma[["estimate"]][2]


args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop("Usage: Rscript simu.R samplesize featuresize fc mr_all mr_nrandom round")
}

samplesize <- as.integer(args[1])
featuresize <- as.integer(args[2])
fc <- as.numeric(args[3])
mr_all <- as.numeric(args[4])
mr_nrandom <- as.numeric(args[5])
i <- as.integer(args[6]) # round


diff_ratio<-0.5
simu_start_time = proc.time()

simu_result_dir = paste0(cd, "/simu_results/Simulation_", samplesize,
                                        "_",featuresize,
                                        "_",fc, 
                                        "_",mr_all, 
                                        "_",mr_nrandom)
if (!dir.exists(simu_result_dir)) {dir.create(simu_result_dir, recursive = TRUE)}                                        
simu_oks = list.files(simu_result_dir)

print(paste("Simulation_", samplesize,
                    "_", featuresize,
                    "_", fc, 
                    "_", mr_all, 
                    "_", mr_nrandom))

if (!paste0("round_", i, ".txt") %in% simu_oks){
  print(paste0("round_", i, " start"))
  repeat_start_time = proc.time()
  # random seeds
  set.seed(samplesize*featuresize*fc*(exp(mr_all)**exp(mr_nrandom))*i)

  data_return<-generate_matrix_simu(mu_0=mu_0,sigma_0=sigma_0,shape=shape,scale=scale,
                          fc=fc,samplesize=samplesize,featuresize=featuresize,
                          mr_all=mr_all,mr_nrandom=mr_nrandom,mr_random=1-mr_nrandom,
                          diff_ratio=diff_ratio)
  complete_data = data_return[[1]]
  # write complete data
  complete_data_dir = paste0(cd, "/simu_complete_data/Simulation_", samplesize,
                                        "_",featuresize,
                                        "_", fc, 
                                        "_", mr_all, 
                                        "_",mr_nrandom)
  if (!dir.exists(complete_data_dir)) dir.create(complete_data_dir, recursive = TRUE)
  write.table(complete_data, paste0(complete_data_dir, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)

  a = data_return[[2]]
  # write data
  simu_data_dir = paste0(cd, "/simu_data/Simulation_", samplesize,
                                        "_",featuresize,
                                        "_", fc, 
                                        "_", mr_all, 
                                        "_",mr_nrandom)
  if (!dir.exists(simu_data_dir)) dir.create(simu_data_dir, recursive = TRUE)
  write.table(a, paste0(simu_data_dir, "/data.na_round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)

  outputs = simu1round(a)

  # write result
  write.table(outputs[[1]], paste0(simu_result_dir, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)

  # write result raw
  simu_result_dir_raw = paste0(cd, "/simu_results/Simulation_rawpv_", samplesize,
                                                    "_",featuresize,
                                                    "_", fc,
                                                    "_", mr_all,
                                                    "_",mr_nrandom)

  if (!dir.exists(simu_result_dir_raw)) {dir.create(simu_result_dir_raw, recursive = TRUE)}
  write.table(outputs[[2]], paste0(simu_result_dir_raw, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)

  # write data BPCA
  simu_data_dir_bpca = paste0(cd, "/simu_data_bpca/Simulation_", samplesize,
                                        "_",featuresize,
                                        "_", fc,
                                        "_", mr_all,
                                        "_",mr_nrandom)
  if (!dir.exists(simu_data_dir_bpca)) dir.create(simu_data_dir_bpca, recursive = TRUE)
  write.table(outputs[[3]], paste0(simu_data_dir_bpca, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)

  # write data sampmin
  simu_data_dir_sampmin = paste0(cd, "/simu_data_sampmin/Simulation_", samplesize,
                                        "_",featuresize,
                                        "_", fc,
                                        "_", mr_all,
                                        "_",mr_nrandom)
  if (!dir.exists(simu_data_dir_sampmin)) dir.create(simu_data_dir_sampmin, recursive = TRUE)
  write.table(outputs[[4]], paste0(simu_data_dir_sampmin, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)

  # write data QR
  simu_data_dir_qr = paste0(cd, "/simu_data_qr/Simulation_", samplesize,
                                        "_",featuresize,
                                        "_", fc,
                                        "_", mr_all,
                                        "_",mr_nrandom)
  if (!dir.exists(simu_data_dir_qr)) dir.create(simu_data_dir_qr, recursive = TRUE)
  write.table(outputs[[5]], paste0(simu_data_dir_qr, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)

  # write data rf
  simu_data_dir_rf = paste0(cd, "/simu_data_rf/Simulation_", samplesize,
                                        "_",featuresize,
                                        "_", fc,
                                        "_", mr_all,
                                        "_",mr_nrandom)
  if (!dir.exists(simu_data_dir_rf)) dir.create(simu_data_dir_rf, recursive = TRUE)
  write.table(outputs[[6]], paste0(simu_data_dir_rf, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)

  # write data knn
  simu_data_dir_knn = paste0(cd, "/simu_data_knn/Simulation_", samplesize,
                                        "_",featuresize,
                                        "_", fc,
                                        "_", mr_all,
                                        "_",mr_nrandom)
  if (!dir.exists(simu_data_dir_knn)) dir.create(simu_data_dir_knn, recursive = TRUE)
  write.table(outputs[[7]], paste0(simu_data_dir_knn, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)
  
  # method time
  time_data_dir = paste0(cd, "/simu_time/Simulation_", samplesize,
                         "_",featuresize,
                         "_", fc, 
                         "_", mr_all, 
                         "_",mr_nrandom)
  if (!dir.exists(time_data_dir)) dir.create(time_data_dir, recursive = TRUE)
  write.table(outputs[[8]], paste0(time_data_dir, "/round_", i, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)
  
# total time
now_time = proc.time()
all_elapsed = now_time - simu_start_time
repeat_elapsed = now_time - repeat_start_time
appender(level="info", paste("Total elapsed", all_elapsed[["elapsed"]] / 60, "min."))
appender(level="info", paste("Simulation_", samplesize,
                                        "_", featuresize,
                                        "_", fc,
                                        "_", mr_all,
                                        "_", mr_nrandom,
                                        "_repeat:", i,
                                        "elapsd", repeat_elapsed[["elapsed"]], "s."))
print(paste0("round_", i, " end"))
} else {
  appender(level="info", paste("Simulation_", samplesize,
                                          "_", featuresize,
                                          "_", fc,
                                          "_", mr_all,
                                          "_", mr_nrandom,
                                          "repeat",i,
                                          "already existed"))
  print(paste0("round_", i, " already existed"))
}