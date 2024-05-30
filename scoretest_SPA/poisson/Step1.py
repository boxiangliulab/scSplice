from Step1_tools import *
from Prepare_data import *

def GLMM(clu_index, path_pca, path_expression,family):
    X, y = read_cov_and_pheno(clu_index, path_pca, path_expression)
    null_model_result = fit_null_model(X, y, family)

    return null_model_result