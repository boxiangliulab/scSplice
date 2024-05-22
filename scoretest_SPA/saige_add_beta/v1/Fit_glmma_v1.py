import numpy as np
import pandas as pd

from SAIGE_geno import geno
from SAIGE.SAIGE_FitGLMM_tools import *
from SAIGE_DataPreparation_tools import *


def fit_glmma(prefix, 
              Tau_dic, 
              verbose=True, 
              isDiagofKinSetAsOne=True, 
              maxiterPCG=100, tolPCG=1e-6, maxiter=100, tol=1e-4, traceCVcutoff=1e-4, 
              LOCO=False):
    '''
    Function to fit GLMM

    Args:
        subPheno (character): data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers);
        Tau_dic (dictionary): The parameter, vector for initial values for the variance component parameter estimates;
        Verbose (logi): Whether to print out the intermediate information;
        isDiagofKinSetAsOne (logi): Whether set the diag elements to 1;
        maxiterPCG (int): Max iteration times for PCG;
        tolPCG (float): The tolerance of whether break down the PCG;
        maxiter (int): Max iteration times for optimize Tau;
        tol (float): The tolerance of whether break down the optimize Tau loop;
        traceCVcutoff (float): Threshold for the coefficient of variation for trace estimation;
    
    Returns:
        A dictionary contains all glmma optimize result.

    ''''

        # 初始化需要调用的类
    X, offset, y = read_in_covariate(prefix + '.csv')
    if verbose:
        print("Phenotype and Covariate reading is done")

    Geno = geno()
    t_begin = pd.to_datetime("today")

    if verbose:
        print(t_begin)

    subSampleInGeno = subPheno['IndexGeno']

    if verbose:
        print("Start reading genotype plink file here")

    # setgeno function - Uncomment and implement this function if needed
    # 读取genofile中的数据
    # re1 = time.time({Geno.setgeno(genofile, subSampleInGeno, isDiagofKinSetAsOne)})



    '''Fit null GLMM model'''
    # subpheno_path = "/home/svu/e0950333/SAIGE/sample_data/subPheno.csv"
    data_path = "/home/svu/e0950333/SAIGE/sample_data/data_new.csv"
    fit_null_model_result = fit_nullGLMM_model(data_path, family)
    y = fit0['y'] # 对应fit0$y
    n = fit0['n'] # n是有多少个样本量
    X = fit0['X'] # X = model.matrix(fit0)
    offset = fit0['offset']
    if offset is None:
        offset = np.zeros(n)

    family = fit0['family']
    eta = fit0['eta']
    mu = fit0['mu']
    mu_eta = fit0['mu_eta']

    Y = eta - offset + (y - mu) / mu_eta
    alpha0 = fit0['coef']
    eta0 = eta

    if family in ["poisson", "binomial"]:
        tau[0] = 1
        fixtau[0] = 1
    
    

    '''Optimzie ML for the fix-coefficients alpha and random effects b'''
    Geno =  read_and_create_GRM(prefix + '.vcf', convert_data)
    ValueVec = Geno.data
    row = Geno.row
    col = Geno.col
    LocationMat = np.concatenate((Geno.row, Geno.col), axis = 1)

    fit_tools = fitGLMM_tools(Geno, ValueVec, LocationMat, dimNum)
    re_coef = fit_tools.get_coef(y, X, tau, family, alpha0, eta0, offset, verbose=verbose, maxiterPCG=maxiterPCG, tolPCG=tolPCG, maxiter=maxiter)
    
    '''AI_REML for variance components (phi, tau)'''
    nrun = 10
    re = fit_tools.getAIScore(re_coef['Y'], X, re_coef['W'], tau, re_coef['Sigma_iY'], re_coef['Sigma_iX'], re_coef['cov'],
                               nrun, maxiterPCG, tolPCG=tolPCG, traceCVcutoff=traceCVcutoff)
    tau[1] = max(0, tau0[1] + tau0[1]**2 * (re['YPAPY'] - re['Trace']) / n)

    if verbose:
        print("Variance component estimates:", tau)

    for i in range(maxiter):

        alpha0 = re_coef['alpha']
        tau0 = tau
        print("tau0_v1:", tau0)

        eta0 = eta

        # use Get_Coef before getAIScore
        t_begin_get_coef = pd.to_datetime("today")
        # 获取稳定系数
        re_coef = fit_tools.get_coef(y, X, tau, family, alpha0, eta0, offset, verbose=verbose, maxiterPCG=maxiterPCG,
                           tolPCG=tolPCG, maxiter=maxiter)
        t_end_get_coef = pd.to_datetime("today")
        print("t_end_get_coef - t_begin_get_coef")
        print(t_end_get_coef - t_begin_get_coef)

        fit = fit_tools.fit_glmmaiRPCG(re_coef['Y'], X, re_coef['W'], tau, re_coef['Sigma_iY'], re_coef['Sigma_iX'],
                             re_coef['cov'], nrun, maxiterPCG, tolPCG, tol=tol, traceCVcutoff=traceCVcutoff)

        t_end_fit_glmmaiRPCG = pd.to_datetime("today")
        print("t_end_fit_glmmaiRPCG - t_end_get_coef")
        print(t_end_fit_glmmaiRPCG - t_end_get_coef)

        tau = [float(i) for i in fit['tau']]
        tau = np.array(tau)
        tau0 = np.array(tau0)
        cov = re_coef['cov']
        alpha = re_coef['alpha']
        eta = re_coef['eta']
        Y = re_coef['Y']
        mu = re_coef['mu']

        print(np.abs(tau - tau0) / (np.abs(tau) + np.abs(tau0) + tol))
        print("tau:", tau)
        print("tau0:", tau0)

        if verbose:
            print("\nIteration", i, ":", tau)

        if tau[1] == 0:
            break

        # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
        if np.max(np.abs(tau - tau0) / (np.abs(tau) + np.abs(tau0) + tol)) < tol:
            break

        if np.max(tau) > tol**(-2):
            Warning("Large variance estimate observed in the iterations, model not converged...", call=False)
            i = maxiter
            break

    if verbose:
        print("Final", tau, ":")

    # added these steps after tau is estimated 04-14-2018
    re_coef = fit_tools.get_coef(y, X, tau, family, alpha, eta, offset, verbose=verbose, maxiterPCG=maxiterPCG,
                       tolPCG=tolPCG, maxiter=maxiter)
    cov = re_coef['cov']
    alpha = re_coef['alpha']
    eta = re_coef['eta']
    Y = re_coef['Y']
    mu = re_coef['mu']

    converged = i < maxiter
    res = y - mu

    # if(isCovariateTransform & hasCovariate){   out.transform
    
    # if out.transform is not None:
    #     coef_alpha = covariate_transform_back(alpha, out.transform['Param.transform'])
    # else:
    #     coef_alpha = alpha
    coef_alpha = alpha
    mu2 = mu * (1 - mu)
    # obj_noK = score_test_null_model(mu, mu2, y, X)

    glmm_result = {
        'theta': tau,
        'coefficients': coef_alpha,
        'linear.predictors': eta,
        'fitted.values': mu,
        'Y': Y,
        'residuals': res,
        'cov': cov,
        'converged': converged,
        # 'sampleID': subPheno['IID'],
        # 'obj.noK': obj_noK,
        'y': y,
        'X': X,
        'traitType': "binary"
    }

    # if LOCO:
    #     set_diagof_stdgeno_LOCO()
    #     glmm_result['LOCO'] = LOCO
    #     glmm_result['LOCOResult'] = []

    #     for j in range(1, 23):
    #         start_index = chromosome_start_index_vec[j]
    #         end_index = chromosome_end_index_vec[j]

    #         if not pd.isna(start_index) and not pd.isna(end_index):
    #             print("leave chromosome", j, "out")
    #             set_start_end_index(start_index, end_index, j - 1)

    #             t_begin_get_coef_loco = pd.to_datetime("today")
    #             re_coef_loco = get_coef_loco(y, X, tau, family, alpha, eta, offset, verbose=verbose,
    #                                          maxiterPCG=maxiterPCG, tolPCG=tolPCG, maxiter=maxiter)
    #             t_end_get_coef_loco = pd.to_datetime("today")
    #             print("t_end_get_coef_loco - t_begin_get_coef_loco")
    #             print(t_end_get_coef_loco - t_begin_get_coef_loco)

    #             cov = re_coef_loco['cov']
    #             alpha = re_coef_loco['alpha']
    #             eta = re_coef_loco['eta']
    #             Y = re_coef_loco['Y']
    #             mu = re_coef_loco['mu']

    #             mu2 = mu * (1 - mu)
    #             res = y - mu

    #             if out.transform is not None:
    #                 coef_alpha = covariate_transform_back(alpha, out.transform['Param.transform'])
    #             else:
    #                 coef_alpha = alpha

    #             obj_noK = score_test_null_model(mu, mu2, y, X)

    #             glmm_result['LOCOResult'].append({
    #                 'isLOCO': True,
    #                 'coefficients': coef_alpha,
    #                 'linear.predictors': eta,
    #                 'fitted.values': mu,
    #                 'Y': Y,
    #                 'residuals': res,
    #                 'cov': cov,
    #                 'obj.noK': obj_noK
    #             })

    #         else:
    #             glmm_result['LOCOResult'].append({'isLOCO': False})

    t_end_null = pd.to_datetime("today")
    print("t_end_null - t_begin, fitting the NULL model without LOCO took")
    print(t_end_null - t_begin)

    return glmm_result
