rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(snowfall)
library(dplyr)
library(tidyr)


setwd($your_working_dir)
# result = read.csv('./simu_metrics/simu_metrics_100_10000_0.75_simulation1.txt',sep='\t')
result = read.csv('./simu_metrics/simu_metrics_simulation2.txt',sep='\t')

result = result %>% filter(!model %in% c("Modt_SampMin","Modt_bPCA","Modt_QR","Modt_RF","Modt_KNN"))

result[,"MNAR"]<-factor(result[,"MNAR"])
result[result$model %in% c("twoT","twoWilcox","ModT","SDA","AFT","Mixture","DASEV"),'Model_Type'] = 'Imputation-free'
result[result$model %in% c("T-test","T_SampMin","T_bPCA","T_QR","T_RF","T_KNN"),'Model_Type'] = 'Parametric'
result[result$model %in% c("Wilcoxon","Wilcox_SampMin","Wilcox_bPCA","Wilcox_QR","Wilcox_RF","Wilcox_KNN"),'Model_Type'] = 'Non-parametric'

result[result$model=='T-test','model'] = 'T'
result[result$model=='T_SampMin','model'] = 'T_SMin'
result[result$model=='Wilcoxon','model'] = 'Wilcox'
result[result$model=='Wilcox_SampMin','model'] = 'Wilcox_SMin'
result[result$model=='Wilcox_bPCA','model'] = 'Wilcox_bPCA'
result[result$model=='Wilcox_QR','model'] = 'Wilcox_QR'
result[result$model=='Wilcox_RF','model'] = 'Wilcox_RF'
result[result$model=='Wilcox_KNN','model'] = 'Wilcox_KNN'

result[,"model"]<-factor(result[,"model"],levels=c("ModT","twoT","twoWilcox","T",
                                                   "Wilcox","AFT","SDA","DASEV","Mixture",
                                                   "T_SMin","T_QR","T_KNN","T_bPCA","T_RF",
                                                   "Wilcox_SMin","Wilcox_QR","Wilcox_KNN","Wilcox_bPCA","Wilcox_RF"))

p1 = ggplot(result,aes(model,pAUROC,fill=MNAR))+
  geom_boxplot()+
  xlab(NULL)+
  ylab("pAUROC")+
  guides(fill=guide_legend(title="MNAR",
                           keyheight = 1.25,
                           keywidth = 0.75,
                           label.theme = element_text(face="bold",size=18)))+
  facet_grid(~Model_Type, scale = "free_x", space = "free") +
  scale_fill_manual(values=c('0.2'="#2A9471", '0.4'="#619CFF", '0.6'="#F64926",
                             '0.8'="#BE491C"))+
  theme_bw() +
  theme(axis.title = element_text(face="bold",size=16),
        axis.text.y = element_text(face="bold",size=16),
        axis.text.x = element_text(face="bold",size=16,angle = 30,hjust = 1, vjust = 1),
        legend.title = element_text(face="bold",size=16),
        axis.ticks.x = element_blank(),
        strip.text.x.top = element_text(size=14),
        panel.grid.major = element_line(color = "gray70", size = .5),
        panel.grid.minor = element_line(color = "gray70", size = .25))

p2 = ggplot(result,aes(model,FPR,fill=MNAR))+
  geom_boxplot()+
  xlab(NULL)+
  ylab("FPR")+
  guides(fill=guide_legend(title="MNAR",
                           keyheight = 1.25,
                           keywidth = 0.75,
                           label.theme = element_text(face="bold",size=18)))+
  facet_grid(~Model_Type, scale = "free_x", space = "free") +
  scale_fill_manual(values=c('0.2'="#2A9471", '0.4'="#619CFF", '0.6'="#F64926",
                             '0.8'="#BE491C"))+
  theme_bw() +
  theme(axis.title = element_text(face="bold",size=16),
        axis.text.y = element_text(face="bold",size=16),
        axis.text.x = element_text(face="bold",size=16,angle = 30,hjust = 1, vjust = 1),
        legend.title = element_text(face="bold",size=16),
        axis.ticks.x = element_blank(),
        strip.text.x.top = element_text(size=14),
        panel.grid.major = element_line(color = "gray70", size = .5),
        panel.grid.minor = element_line(color = "gray70", size = .25))

ggsave("./plot/100_10000_0.75_pAUROC.png",
       p1,
       device="png",
       dpi=1200,
       height=6,
       width = 12)
ggsave("./plot/100_10000_0.75_fpr.png",
       p2,
       device="png",
       dpi=1200,
       height=6,
       width = 12)

# find exact diff values for ProteoCardis
result_median = aggregate(result$pAUROC, by=list(result$Nsample,result$fc,result$zr,result$model, result$MNAR, result$Model_Type), FUN=function(x) median(x,na.rm=T))
colnames(result_median) = c("samplesize","fc","mr","method","mnar", "Model_Type","pAUROC")
result_median_diff = pivot_wider(result_median, id_cols = c(samplesize,fc,mr,method,Model_Type), names_from = c(mnar), values_from = c(pAUROC))
result_median_diff$diff_0.2_0.8 = result_median_diff$`0.2` - result_median_diff$`0.8`
result_median_diff$diff_0.2_0.8_ratio = result_median_diff$diff_0.2_0.8 / result_median_diff$`0.2`
mean(result_median_diff$diff_0.2_0.8_ratio[result_median_diff$diff_0.2_0.8_ratio < -0.06])

result_median_crossMNAR = aggregate(result$pAUROC, by=list(result$Nsample,result$fc,result$zr,result$model, result$Model_Type), FUN=function(x) mean(x,na.rm=T))
result_sd_crossMNAR = aggregate(result$pAUROC, by=list(result$Nsample,result$fc,result$zr,result$model, result$Model_Type), FUN=function(x) sd(x,na.rm=TRUE))

FPR_median_crossMNAR = aggregate(result$FPR, by=list(result$Nsample,result$fc,result$zr,result$model, result$Model_Type), FUN=function(x) mean(x,na.rm=T))

result_median_0.3_0.2 = result_median %>% filter(mnar == 0.2)
result_median_0.3_0.4 = result_median %>% filter(mnar == 0.4)
result_median_0.3_0.6 = result_median %>% filter(mnar == 0.6)
result_median_0.3_0.8 = result_median %>% filter(mnar == 0.8)

# plot average pAUROC for different samples
result = result %>% filter(!is.na(pAUROC))
result_ave_across_fc = aggregate(result$pAUROC, by=list(result$Nsample,result$zr,result$model), FUN=mean)
colnames(result_ave_across_fc) = c("samplesize","mr","model","pAUROC")
result_ave_across_fc_max = result_ave_across_fc %>%
  group_by(samplesize,mr) %>%
  filter(pAUROC == max(pAUROC)) %>%
  ungroup()
result_ave_across_fc_max$mr = as.factor(result_ave_across_fc_max$mr)
result_ave_across_fc_max$samplesize = as.factor(result_ave_across_fc_max$samplesize)

g<-ggplot(result_ave_across_fc_max, aes(x = samplesize, y = mr, fill = model)) +
  geom_tile(color="white",lwd=2) +
  xlab("Sample Size")+
  ylab("Missing Ratio")+
  theme_bw()+
  geom_text(aes(label = sprintf("%.3f", pAUROC)), vjust = 0.5, size = 5)+  # 在每个瓦片上标注pAUROC值
  guides(fill=guide_legend(title="Method",
                           label.theme = element_text(size=10)))+
  theme(axis.title = element_text(face="bold",size=12),
        axis.text = element_text(size=10),
        plot.title = element_text(size=14,face="bold",hjust=0.5))+
  scale_fill_manual(values=c(Wilcox = "#ED342F", T = "#19759f",
                             ModT = "#4fb587" , AFT = "#93be5b",DASEV ="#fff1cf",
                             twoT = "#a0add0",
                             T_RF = "#8583a9" , T_KNN = "#D2878c",T_bPCA="#fed865",
                             Wilcox_RF = "#9493ff",Wilcox_bPCA="#f58b6a"))

ggsave(paste(c("./plot/recommend_samplesize_missingratio.png"),collapse = ""),g,device="png",dpi=1200,height=4,width = 8)
