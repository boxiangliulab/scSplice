# -*- coding: UTF-8 -*-

import sys
from Step1 import *
from Step2 import *
import time
import numpy as np
import warnings
import csv

warnings.filterwarnings("ignore")
clu_index = 0
path_pca = "C:\study\poisson\data\simu_genoPC_5_28.txt"
path_expression= "C:\study\poisson\data\metacell_splice_vs_total.txt"
geno_PATH = ""

GLMM_Result = GLMM(clu_index = clu_index,
                   path_pca = path_pca,
                   path_expression = path_expression,
                   family = "Poisson")
print(GLMM_Result)

ScoreTest_Result = ScoreTest(geno_PATH, GLMM_Result)
print(ScoreTest_Result)

