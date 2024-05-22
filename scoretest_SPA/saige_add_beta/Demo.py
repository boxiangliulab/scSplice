# -*- coding: UTF-8 -*-

"""pheno1 geno64;  pheno 5 geno 895"""
import sys
from Step1_FittGLMM import Step1_FitGLMM
from Step2_SPAtest import Step2_SPAtest
import time
import numpy as np
import warnings
import csv

warnings.filterwarnings("ignore")
Geno_name = '21'
Cell_name = 'CD4+_T_naive'
pheno_index = 79 # 哪一条phenotype的指标，跑实验的话需要循环起来 80


# geno_index = [4,12,14,15,20,22,23,26,27,29,30,37,38,48,50,54,57,64,66,68,71,75,79,81,83,86,90,92,94,96,97,99] # 哪一条genotype来做spa test
# geno_index = range(101,500)
geno_index = [1]
step1_start_time = time.time()

GLMM_Result = Step1_FitGLMM(#prefix = sys.path[0] + "/ZYTExperiment_data/",
                            # prefix = sys.path[0] + "/Experiment_data/",
                            prefix = sys.path[0] + "/simulation_data/",
                            Geno = Geno_name,
                            Cell = Cell_name,
                            index = pheno_index,
                            M = 20,
                            # family_indicator = 'Binomial',
                            family_indicator = 'Beta',
                            Tau_dic = {"tauInit": [1, 1], "tau": [1, 1], "fixtau": [1, 0]},
                            maxiterPCG=100, tolPCG=1e-5, maxiter=100, tol=1e-2, traceCVcutoff=1e-2,
                            verbose=True)

np.savetxt(sys.path[0] + "/fitted_values.csv", GLMM_Result['tau'], delimiter="," )
step1_end_time = time.time()





with open('pvalues.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['geno_index', 'p_value', 'p_adjusted'])

    for i in geno_index:

        pval_noadjust, pval = Step2_SPAtest(GLMM_Result, pheno_index, i)
        print('p-value for index {}: {}'.format(i, pval))
        print('p_adjusted for index {}: {}'.format(i, pval_noadjust))
        writer.writerow([i, pval, pval_noadjust])

step2_end_time = time.time()
print("step1:",step1_end_time-step1_start_time)
print("step2:",step2_end_time-step1_end_time)





#
#
#
#
# import sys
# from Step1_FittGLMM import Step1_FitGLMM
# from Step2_SPAtest import Step2_SPAtest
# import time
# import numpy as np
# import warnings
# import csv
#
# warnings.filterwarnings("ignore")
#
# # Initialization
# Geno_name = '21'
# Cell_name = 'CD4+_T_naive'
# num_pheno = 79  # Total number of phenotypes
# num_geno = 1000  # Total number of genotypes
#
# # Data structures
# p_adjusted_matrix = np.zeros((num_pheno, num_geno))
#
# # Iterate over all phenotypes
# for pheno_index in range(num_pheno):
#     # print("Step1_FitGLMM index", pheno_index)
#     # Step 1: Fit GLMM
#     step1_start_time = time.time()
#     GLMM_Result = Step1_FitGLMM(
#         prefix = sys.path[0] + "/simulation_data/",
#         Geno = Geno_name,
#         Cell = Cell_name,
#         index = pheno_index,
#         M = 20,
#         family_indicator = 'Beta',
#         Tau_dic = {"tauInit": [1, 1], "tau": [1, 1], "fixtau": [1, 0]},
#         maxiterPCG = 100, tolPCG = 1e-5, maxiter = 100, tol = 1e-2, traceCVcutoff = 1e-2,
#         verbose = True
#     )
#     # step1_end_time = time.time()
#     # print(f"Step 1 for phenotype {pheno_index} completed in {step1_end_time - step1_start_time} seconds.")
#
#     # Iterate over all genotypes
#     for geno_index in range(num_geno):
#         pval_noadjust, pval = Step2_SPAtest(GLMM_Result, pheno_index, geno_index+1)
#         p_adjusted_matrix[pheno_index, geno_index] = pval_noadjust
#
#
# # Save results to CSV
# with open('pvalues.csv', 'w', newline='') as file:
#     writer = csv.writer(file)
#     # Write header
#     writer.writerow(          ['Pheno_Index'] + list( range(1,num_geno+1) )              )
#     # Write data
#     for pheno_index in range(num_pheno):
#         writer.writerow( [pheno_index+1] + list(p_adjusted_matrix[pheno_index]))
#
# print("Data processing complete.")
