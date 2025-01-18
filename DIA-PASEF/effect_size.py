import pandas as pd
import os
import numpy as np

cd = $your_working_dir

repeat_time=500
fcs=[2]
samplesizes=[24]
featuresizes=[10000]
mr_alls=[0.3]
mr_nrandoms=[0.2,0.4,0.6,0.8]

def calc_ef_size(df):
    mean_control = df.iloc[:,0:int(df.shape[-1]/2)].mean(axis=1,skipna=True)
    mean_case = df.iloc[:,int(df.shape[-1]/2):].mean(axis=1,skipna=True)
    sd = df.std(axis=1,skipna=True)
    ef = (mean_case - mean_control) / sd
    return pd.DataFrame({'mean_control':mean_control,
                            'mean_case':mean_case,
                            'sd':sd,
                            'ef':ef})

for fc in fcs:
  for samplesize in samplesizes:
    for featuresize in featuresizes:
      for mr_all in mr_alls:
        for mr_nrandom in mr_nrandoms:
          delta_na = pd.Series()
          delta_bpca = pd.Series()
          delta_sampmin = pd.Series()
          delta_qr = pd.Series()
          delta_rf = pd.Series()
          delta_knn = pd.Series()
          labels = pd.Series()
          for i in range(repeat_time):
            i = i+1
            data_path = f"{cd}/simu_complete_data/Simulation_{samplesize}_{featuresize}_{fc}_{mr_all}_{mr_nrandom}/round_{i}.txt"
            data_na_path = f"{cd}/simu_data/Simulation_{samplesize}_{featuresize}_{fc}_{mr_all}_{mr_nrandom}/data.na_round_{i}.txt"
            data_bpca_path = f"{cd}/simu_data_bpca/Simulation_{samplesize}_{featuresize}_{fc}_{mr_all}_{mr_nrandom}/round_{i}.txt"
            data_sampmin_path = f"{cd}/simu_data_sampmin/Simulation_{samplesize}_{featuresize}_{fc}_{mr_all}_{mr_nrandom}/round_{i}.txt"
            data_qr_path = f"{cd}/simu_data_qr/Simulation_{samplesize}_{featuresize}_{fc}_{mr_all}_{mr_nrandom}/round_{i}.txt"
            data_rf_path = f"{cd}/simu_data_rf/Simulation_{samplesize}_{featuresize}_{fc}_{mr_all}_{mr_nrandom}/round_{i}.txt"
            data_knn_path = f"{cd}/simu_data_knn/Simulation_{samplesize}_{featuresize}_{fc}_{mr_all}_{mr_nrandom}/round_{i}.txt"
            
            data = pd.read_table(data_path,sep='\t',header=0)
            label = data.iloc[:,-1]
            data = data.iloc[:,0:data.shape[-1]-1]
            
            data_na = pd.read_table(data_na_path,sep='\t',header=0)
            data_na = data_na.iloc[:,0:data_na.shape[-1]-1]
            data_bpca = pd.read_table(data_bpca_path,sep='\t',header=0)
            data_sampmin = pd.read_table(data_sampmin_path,sep='\t',header=0)
            data_qr = pd.read_table(data_qr_path,sep='\t',header=0)
            data_rf = pd.read_table(data_rf_path,sep='\t',header=0)
            data_knn = pd.read_table(data_knn_path,sep='\t',header=0)
            if not data.shape==data_na.shape==data_bpca.shape==data_sampmin.shape==data_qr.shape==data_rf.shape==data_knn.shape:
              print(data.shape,data_na.shape,data_bpca.shape,data_sampmin.shape,data_qr.shape,data_rf.shape,data_knn.shape)
              print(f"round_{i}")
              exit(1)

            # filter out result when each group is all NA
            mask_control = (~pd.isna(data_na.iloc[:,:int(samplesize/2)])).sum(axis=1)!=0
            mask_case = (~pd.isna(data_na.iloc[:,int(samplesize/2):samplesize])).sum(axis=1)!=0
            filter_flag = np.logical_and(mask_control,mask_case)
            print(f"round {i} filtered {sum(~filter_flag)}")
            data = data.loc[filter_flag,]
            data_na = data_na.loc[filter_flag,]
            data_sampmin = data_sampmin.loc[filter_flag,]
            data_qr = data_qr.loc[filter_flag,]
            data_rf = data_rf.loc[filter_flag,]
            data_knn = data_knn.loc[filter_flag,]
            label = label.loc[filter_flag,]
            
            efdata = calc_ef_size(data)
            efdata_na = calc_ef_size(data_na)
            efdata_sampmin = calc_ef_size(data_sampmin)
            efdata_qr = calc_ef_size(data_qr)
            efdata_rf = calc_ef_size(data_rf)
            efdata_knn = calc_ef_size(data_knn)
            
            delta_na = pd.concat([delta_na, efdata_na.ef - efdata.ef], axis=0)
            delta_sampmin = pd.concat([delta_sampmin, efdata_sampmin.ef - efdata.ef], axis=0)
            delta_qr = pd.concat([delta_qr, efdata_qr.ef - efdata.ef], axis=0)
            delta_rf = pd.concat([delta_rf, efdata_rf.ef - efdata.ef], axis=0)
            delta_knn = pd.concat([delta_knn, efdata_knn.ef - efdata.ef], axis=0)

            if (sum(data_bpca.shape) != 2):
              data_bpca = data_bpca.loc[filter_flag,]
              efdata_bpca = calc_ef_size(data_bpca)
              delta_bpca = pd.concat([delta_bpca, efdata_bpca.ef - efdata.ef], axis=0)
            else:
              delta_bpca = pd.concat([delta_bpca, pd.Series([np.nan] * efdata.shape[0])], axis=0)
              print(delta_bpca)

            labels = pd.concat([labels,label], axis=0)

          delta_na = delta_na.to_list()
          delta_sampmin = delta_sampmin.to_list()
          delta_qr = delta_qr.to_list()
          delta_rf = delta_rf.to_list()
          delta_knn = delta_knn.to_list()
          delta_bpca = delta_bpca.to_list()
          labels = labels.to_list()
          print(len(delta_na),len(delta_sampmin),len(delta_qr),len(delta_rf),len(delta_knn),len(delta_bpca),len(labels))
          out_data = pd.DataFrame({'delta_na':delta_na,
                                      'delta_sampmin':delta_sampmin,
                                      'delta_qr':delta_qr,
                                      'delta_rf':delta_rf,
                                      'delta_knn':delta_knn,
                                      'delta_bpca':delta_bpca,
                                      'label':labels})
          print('Writing.')
          out_data.to_csv(f"{cd}/effect_size/Simulation_{samplesize}_{featuresize}_{fc}_{mr_all}_{mr_nrandom}.txt",
                          sep='\t',
                          index=False)