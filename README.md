# Metaproteomics-MV-simulations
Evaluation of differential abundance analysis strategies for metaproteomics data with missing values.

We evaluated five commonly-used imputation methods and six imputation-free methods in two simulations of both DDA (data dependent acquisition) and DIA (data independent acquisition) metaproteomics datasets, including 588 data scenarios of various levels of sample size, fold change, missing ratio and MNAR ratio. Each scenario was simulated for 500 rounds in simulation 1, which simulated the DDA datasets (ProteoCardis and ALTR), and DIA dataset (DIA-PASEF). In simulation 2, each scenario was simulated for 100 times.

The code for simulation of dataset ProteoCardis, ALTR, and DIA-PASEF are **in the corresponding directory**.

## Requirements
### R (3.6+)
* invgamma
* Bioconductor
* limma
* SDAMS
* pcaMethods
* trust
* ggplot2
* gridExtra
* data.table
* fitdistrplus

### Python (3.7+)
* NumPy (1.19.2+)
* Pandas (1.2.1+)
* Scikit-Learn (0.23.2+)
* SciPy (1.5.2+)

## Simulation
The codes of simulation study are in **simulation.R**. The parameters of different simulations are given in the **simu_args_simulation_(1/2).txt**. You can conduct the simulation by running the **run_simulation.sh**.

The whole simulation work is very computationally intensive in R, so a simple implementation can be tested by the deleting the rows in **simu_args_simulation_(1/2).txt** file we presented.

Note:
1. The method of SDA was modified by cancelling the data-cleaning step of removing features with low numbers of non-zeros.
2. bPCA imputation was modified by lower the _tol_ parameter in _solve_ function to prevent the error of computationally singular. 

## Results & Plotting
The codes of summarizing the results (adjusted _p_-values) and evaluating the different methods for the 500-round / 100-round simulation are in **results_HPC.py**.

The codes of plotting pAUROC and FDR figures are in **plot.R**.

Calculation of the delta effect size between complete data, missing data, and imputed data were conducted by **effect_size.py**. The codes of plotting the result are in **effect_size.R**.


Only the results data for generating the figures is in this repository due to the limit of file size.
The results data of all methods under all scenarios had huge file size and can not be posted on github. It is available upon request. 
If you have problems, please contact jing.li@sjtu.edu.cn
