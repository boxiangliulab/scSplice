from Prepare_data import  *

import numpy as np
import pandas as pd
import statsmodels.api as sm
import numpy as np
from scipy.stats import chi2, norm,binom, poisson


def score_test_null_model(X, V,alpha):
    first_column = X[:, 0:1]
    alpha_column = np.full_like(first_column, alpha)
    X = np.hstack((alpha_column, X))

    XV = (X * V).T
    XVX_inv = np.linalg.inv(X.T @ (X * V))
    XXVX_inv = X @ XVX_inv

    return XV,XVX_inv,XXVX_inv

def fit_null_model(X,y, family):
    if family == "Poisson":
        poisson_model = sm.GLM(y, sm.add_constant(X), family=sm.families.Poisson())
        poisson_results = poisson_model.fit()
        mu = poisson_results.fittedvalues
        eta = np.exp(mu)
        print(poisson_results.summary())
        alpha = poisson_results.params[0]
        mu_eta = mu
        res = y-mu
        V = mu
        V = V.reshape(-1, 1)
        XV, XVX_inv, XXVX_inv  = score_test_null_model(X, V,alpha)

        null_model_result = {
            "x":X,
            "y": y,
            "mu":mu,
            "eta":eta,
            "mu_eta":mu_eta,
            "res":res,
            "XV":XV,
            "XVX_inv":XVX_inv,
            "XXVX_inv":XXVX_inv
        }
    elif family == "Beta":
        pass


    return null_model_result







# clu_index = 0
# path_pca = "C:\study\poisson\data\simu_genoPC_5_28.txt"
# path_expression= "C:\study\poisson\data\metacell_splice_vs_total.txt"
# X,y = read_cov_and_pheno(clu_index,path_pca,path_expression)
# family = "Poisson"
# null_model_result = fit_null_model(X,y, family)
#
# print(null_model_result)


