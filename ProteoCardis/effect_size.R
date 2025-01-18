rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(snowfall)
library(dplyr)
library(reshape2)
setwd($your_working_dir)
cd = $your_working_dir

fcs<-c(2)
samplesizes<-c(100)
featuresizes<-c(10000)
mr_alls<-c(0.75)
mr_nrandoms<-c(0.2,0.4,0.6,0.8)


# run the python script 'effect_size.py' to compute effect size in each simulation at first

# plot delta effect size
for (fc in fcs){
  for (samplesize in samplesizes){
    for (featuresize in featuresizes){
      for (mr_all in mr_alls){
        for (mr_nrandom in mr_nrandoms){
            data_path = paste0(cd, "/effect_size/Simulation_", samplesize,
                                  "_",featuresize,
                                  "_",fc,
                                  "_",mr_all,
                                  "_",mr_nrandom,
                                  ".txt")
            data = read.csv(data_path,sep = '\t')
            data = melt(data,id.vars=c("label"),variable.name = "method", value.name = "delta_es")
            data[data$label==0,'label'] = 'Negative'
            data[data$label==1,'label'] = 'Positive'
            data$method = as.character(data$method)
            data[data$method=='delta_na','method'] = 'Missing'
            data[data$method=='delta_sampmin','method'] = 'SMin'
            data[data$method=='delta_qr','method'] = 'QR'
            data[data$method=='delta_rf','method'] = 'RF'
            data[data$method=='delta_knn','method'] = 'KNN'
            data[data$method=='delta_bpca','method'] = 'bPCA'
            data$method = factor(data$method,levels=c('Missing','SMin','QR','KNN','bPCA','RF'))
            data$label = as.factor(data$label)
            data[data$method=='']
            p1 = ggplot(data,aes(method,delta_es,fill=label))+
                  geom_boxplot()+
                  xlab(NULL)+
                  ylab("Delta Effect Size")+
                  guides(fill=guide_legend(title="Label",
                                           keyheight = 1.25,
                                           keywidth = 0.75,
                                           label.theme = element_text(face="bold",size=18)))+
                  facet_grid(~label, scale = "free_x", space = "free") +
                  scale_fill_manual(values=c('Negative'="#9493ff", 'Positive'="#fdbfcc"))+
                  scale_y_continuous(limits = c(-4,4))+
                  theme_bw() +
                  theme(axis.title = element_text(face="bold",size=16),
                        axis.text.y = element_text(face="bold",size=16),
                        axis.text.x = element_text(face="bold",size=16,angle = 30,hjust = 1, vjust = 1),
                        legend.title = element_text(face="bold",size=16),
                        axis.ticks.x = element_blank(),
                        strip.text.x.top = element_text(size=14),
                        panel.grid.major = element_line(color = "gray70", size = .5),
                        panel.grid.minor = element_line(color = "gray70", size = .25),
                        legend.position = "none")
            ggsave(paste0(cd,"/plot/effect_size",samplesize,"_",fc,"_",mr_all,"_",mr_nrandom,".png"),
                   p1,
                   device="png",
                   dpi=1200,
                   height = 6,
                   width = 6)
        }
      }
    }
  }
}


# check exact value
for (fc in fcs){
  for (samplesize in samplesizes){
    for (featuresize in featuresizes){
      for (mr_all in mr_alls){
        for (mr_nrandom in c(0.8)){
            data_path = paste0(cd, "/effect_size/Simulation_", samplesize,
                                  "_",featuresize,
                                  "_",fc,
                                  "_",mr_all,
                                  "_",mr_nrandom,
                                  "_filteredAllNAInGroup.txt")
            data = read.csv(data_path,sep = '\t')
            data = melt(data,id.vars=c("label"),variable.name = "method", value.name = "delta_es")
            data[data$label==0,'label'] = 'Negative'
            data[data$label==1,'label'] = 'Positive'
            data$method = as.character(data$method)
            data[data$method=='delta_na','method'] = 'Missing'
            data[data$method=='delta_sampmin','method'] = 'SMin'
            data[data$method=='delta_qr','method'] = 'QR'
            data[data$method=='delta_rf','method'] = 'RF'
            data[data$method=='delta_knn','method'] = 'KNN'
            data[data$method=='delta_bpca','method'] = 'bPCA'
            data$method = factor(data$method,levels=c('Missing','SMin','QR','KNN','bPCA','RF'))
            data$label = as.factor(data$label)
            data[data$method=='']
            data_median = aggregate(data$delta_es, by=list(data$method, data$label), FUN=function(x) median(x,na.rm=T))
        }
      }
    }
  }
}
