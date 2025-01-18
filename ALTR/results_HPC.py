import numpy as np
from sklearn.metrics import roc_curve, precision_recall_curve, auc, roc_auc_score
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests



## The two functions are modified from the ones of: https://github.com/hmsch/proteomics-simulations ##
def roc_prc_scores(is_changed, p_vals):
    """
    Calculates AUROC, AUPRC, and standardized pAUROC scores
    :param is_changed: actual labels (vector of ones and zeros, n x 1)
    :param p_vals:     predicted labels (list of vectors of non-negative pvalues, smaller more significant, k x n)
    :return:           roc_auc, prc_aucs, pauc (all lists of floats)
    """
    roc_auc = []
    prc_auc = []
    pauc = []

    #for pval in p_vals.values.transpose():
    for pval in p_vals.transpose():
        # avoid trouble with -inf in roc_curve
        # replace 0 with an arbitrary small number
        pval = [p if p != 0 else 1e-100 for p in pval]

        predicted = - np.log(pval)

        # e.g. when values are missing because a test could not be applied
        if np.all(np.isnan(predicted)):
            # No valid p-vals
            pauc.append(np.nan)
            roc_auc.append(np.nan)
            prc_auc.append(np.nan)
            continue

        # Rank measurements with nan as lowest priority
        predicted[np.isnan(predicted)] = 0
        
        fpr, tpr, _ = roc_curve(is_changed, predicted)

        pauc_std=roc_auc_score(is_changed, predicted, max_fpr=FDR)

        pauc.append(pauc_std)
        roc_auc.append(auc(fpr, tpr))
        prec, rec, _ = precision_recall_curve(is_changed, predicted)
        prc_auc.append(auc(rec, prec))

    return roc_auc, prc_auc, pauc


def power_analysis(is_changed, pvals, alpha=0.05):
    """
    Calculate Precision and Recall under a fixed pvalue threshold
    :param is_changed: actual labels (vector of ones and zeros, n x 1)
    :param p_vals:     predicted labels (list of vectors of non-negative pvalues, smaller more significant, n x k)
    :param alpha:      pvalue threshold, 0.05 by default
    """
    fps = []
    tps = []
    Precision = []
    Recall = []
    FPR = []
    for pval in pvals.transpose():
        #  e.g. when values are missing because the a test could not be applied
        if np.all(np.isnan(pval)):
            fp, tp, prec, rec, fpr = np.nan, np.nan, np.nan, np.nan, np.nan

        else:
            # take NaN as negatives: (np.nan <= alpha) return False
            # False positives
            fp = np.sum(np.logical_not(is_changed) & (pval <= alpha))
            tp = np.sum((is_changed == 1) & (pval <= alpha))
            p = np.sum(is_changed)
            n = np.sum(np.logical_not(is_changed))
            if p==len(is_changed):
                prec = np.nan
            elif fp+tp==0:
                prec = np.nan
            else:
                prec = tp/(tp+fp)
                
            rec = tp/p
            fpr = fp/n
        fps.append(fp)
        tps.append(tp)
        Precision.append(prec)
        Recall.append(rec)
        FPR.append(fpr)
    return fps, tps, Precision, Recall, FPR


#main
Nmodel=24
method_name=('ModT',
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
            'Wilcox_KNN')
result_all=pd.DataFrame()


CD=$your_working_dir


Nround=500
sample_size=[16] 
feature_sizes=[10000]
fclist=[2]
zero_ratio=[0.3]
mnar_ratio=[0.2,0.4,0.6,0.8]
FDR=0.05
path_result=CD+'/simu_metrics/simu_metrics_16_10000_0.3_simulation1.txt'


for nsample in sample_size:
    for fc in fclist:
        for zr in zero_ratio:
            for mnar in mnar_ratio:
                for feature_size in feature_sizes:
                    # 24 different models
                    result_dir= CD+f'/simu_results/Simulation_rawpv_{nsample}_{feature_size}_{fc}_{zr}_{mnar}/'
                    dataset_dir= CD+f'/simu_data/Simulation_{nsample}_{feature_size}_{fc}_{zr}_{mnar}/'
                    print(result_dir)
                    for r in range(Nround):
                        r=r+1
                        result_path = result_dir + f'round_{r}.txt'
                        try:
                            data_raw=np.genfromtxt(result_path,delimiter='\t',skip_header=1)
                        except:
                            print(f'no {result_path}')
                            continue
                        dataset_path = dataset_dir + f'data.na_round_{r}.txt'
                        dataset = np.genfromtxt(dataset_path,delimiter='\t',skip_header=1)

                        # filter out result when each group is all NA
                        mask_control = (np.sum(~np.isnan(dataset[:,:int(nsample/2)]),axis=1)!=0)
                        mask_case = (np.sum(~np.isnan(dataset[:,int(nsample/2):nsample]),axis=1)!=0)
                        filter_flag = np.logical_and(mask_control,mask_case)
                        dataset = dataset[filter_flag,]
                        data_raw = data_raw[filter_flag,]

                        ## set the NA raw pvalue to the max of the sample to make BH correction fair
                        for i in range(data_raw.shape[1]-1):
                            max_val = np.nanmax(data_raw[:, i])
                            data_raw[np.isnan(data_raw[:, i]), i] = max_val

                        data = np.ones_like(data_raw)
                        for i in range(data_raw.shape[1]-1):
                            rejected, pvals_corrected, _, _ = multipletests(data_raw[:, i], alpha=0.05, method='fdr_bh')
                            data[:, i] = pvals_corrected
                        # label
                        data[:,data_raw.shape[1]-1] = data_raw[:,data_raw.shape[1]-1]

                        actual=data[:,data.shape[1]-1]
                        p_vals=data[:,0:(data.shape[1]-1)]
                        if np.all(actual==1):
                            auroc, auprc, pauroc=[np.nan]*Nmodel,[np.nan]*Nmodel,[np.nan]*Nmodel
                        else:
                            auroc, auprc, pauroc=roc_prc_scores(actual,p_vals)
                        FP, TP, precision, recall, fpr=power_analysis(actual, pvals=p_vals)
                        result_auc=np.array([auroc,pauroc,auprc,precision,recall,fpr])
                        result_scenario=np.array([nsample,fc,zr,mnar,r])
                        result_comb=np.hstack((np.tile(result_scenario,(Nmodel,1)),result_auc.T))                   
                        result_temp=pd.DataFrame(result_comb)
                        result_temp.insert(5,'model',method_name)
                        result_all=pd.concat([result_all,result_temp])
                    

result_all.columns=['Nsample','fc','zr','MNAR','round','model','AUROC','pAUROC','AUPRC','Precision','Recall','FPR']
output_dict={
    'Nsample':'int16',
    'fc':'float16',
    'zr':'float16',
    'MNAR':'float16',
    'round':'int16',
    'model':'object',
    'AUROC':'float64',
    'pAUROC':'float64',
    'AUPRC':'float64',
    'Precision':'float64',
    'Recall':'float64',
    'FPR':'float64'}
result_all=result_all.astype(output_dict)

result_all.to_csv(path_result,sep='\t',index=False)
