import pandas as pd
import scipy.stats as stats
filepath="/data/zhangyuntian/project/scSplice/test/PSI_pheno/CD14+_Monocyte_scSplice_potential_ES"
label=""

with open(filepath,'r') as f:
    for line in f:
        l=line.strip().split(' ')
        for i in range(1,len(l)):


