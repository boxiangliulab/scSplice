import numpy as np
import pandas as pd
import time

# from SAIGE_geno import geno
from SAIGE_FitGLMM_tools import *
from SAIGE_DataPreparation_tools import *
import numpy as np
from scipy.sparse import csr_matrix

def Step1_FitGLMM(prefix, Geno, Cell, index, M, family_indicator,
              Tau_dic,
              maxiterPCG=100, tolPCG=1e-6, maxiter=100, tol=1e-4, traceCVcutoff=1e-4, 
              verbose=True):
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

    '''
    
    t_begin = time.time()


    '''Read in Variables'''
    read_in_tools = Read_in_tools(prefix, Cell, Geno)
    X, y= read_in_tools.Covariate_and_Phenotype(index,family_indicator)



    N = y.shape[0]

    if verbose:
        print("Phenotype and Covariate reading is done")

    '''Optimzie ML for the fix-coefficients alpha and random effects b'''

    """试了一下GRM都设成1，会算出来tau=nan"""
    # row_indices = np.repeat(np.arange(N), N)
    # col_indices = np.tile(np.arange(N), N)
    # data = np.ones(N * N)
    # Geno = csr_matrix((data, (row_indices, col_indices)), shape=(N, N))

    """GRM非对角线元素正态分布"""
    # indices = np.triu_indices(N, 1)
    # off_diagonal_values = np.random.normal(loc=0.5, scale=0.1, size=len(indices[0]))
    # full_matrix = np.zeros((N, N))
    # np.fill_diagonal(full_matrix, 1)
    # full_matrix[indices] = off_diagonal_values
    # full_matrix[(indices[1], indices[0])] = off_diagonal_values  # 镜像到下三角
    # Geno = csr_matrix(full_matrix)



    """读取GRM (saige)--------------------------------------------------"""

    sparseGRMFile = "./simulation_data/sparseGRM_relatednessCutoff_0.125_100_randomMarkersUsed.sparseGRM.mtx"
    sparseGRMSampleIDFile = "./simulation_data/sparseGRM_relatednessCutoff_0.125_100_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
    relatednessCutoff = 0.125
    geno_index = range(101, 501)
    modelID = []
    for i in geno_index:
        modelID.append("S" + str(i))

    Geno = read_in_tools.getsubGRM(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, modelID)
    """读取GRM (from zyt)--------------------------------------------------"""
    # sparseGRMFile = "./Experiment_data/CD4+_T_naive_chr21.grm"
    # sparseGRMSampleIDFile = "./Experiment_data/CD4+_T_naive_chr21.grm.id"
    # Cell_path = sys.path[0]+"/ZYTExperiment_data/Covariate/CD4+_T_naive_PC.txt"
    # X_new = pd.read_csv(Cell_path, sep='\t')
    # X = X_new.iloc[-2:]
    # X = np.array(X.T)[1:, :]
    # modelID = np.array(X_new.columns.values.tolist()[1:])
    # relatednessCutoff = 0.05
    #
    #
    #
    # num_individuals = 100  # 矩阵的维度
    #
    # with open(sparseGRMFile, 'rb') as file:
    #     # num_individuals*num_individuals是矩阵的总元素数量
    #     grm_data = np.fromfile(file, dtype=np.float64, count=num_individuals * num_individuals)
    #     grm_matrix = grm_data.reshape((num_individuals, num_individuals))
    #
    # print(grm_matrix)

    # Geno = read_in_tools.getsubGRM(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, modelID)

    """---------------------------------------------"""


    ValueVec = Geno.data

    # 从稀疏矩阵Geno中提取非零元素的行和列索引
    rows, cols = Geno.nonzero()
    # 将行和列索引重组成位置矩阵
    #LocationMat = np.concatenate((rows.reshape(-1, 1), cols.reshape(-1, 1)), axis=1)
    LocationMat = np.vstack((rows, cols))
    print("LocationMat shape:", LocationMat.shape)

    if verbose:
        print("Genotype reading is done")

    '''Fit null GLMM model'''
    fit_null_model_result = fit_nullGLMM_model(X, y, family_indicator)

    offset = fit_null_model_result['offset']
    if offset is None:
        offset = np.zeros(N)

    eta = fit_null_model_result['eta']
    mu = fit_null_model_result['mu']
    mu_eta = fit_null_model_result['mu_eta']
    Y = eta - offset + (y - mu) / mu_eta

    alpha0 = fit_null_model_result['coef']
    eta0 = eta

    tau = Tau_dic['tau']
    tauInit = Tau_dic['tauInit']
    fixtau = Tau_dic['fixtau']
    if family_indicator in ["Poisson", "Binomial"]:
        tau[0] = 1
        fixtau[0] = 1

    id = 0
    for ii in fixtau:
        if ii != 0 :
            id += 1
        else:
            tau[id] = 0.1
            id += 1

    print("initial tau is ", tau)
    tau0 = tau
    dimNum = N

    glmm_result = {
        'tau': tau,
        'coefficients': alpha0,
        'linear.predictors': eta0,
        'fitted.values': mu,
        'Y': Y,
        # 'residuals': res,
        # 'cov': cov,
        # 'converged': converged,
        'y': y,
        'X': X,
        'traitType': family_indicator,
        # 'sparseGRM': sp.coo_matrix((fit_tools.valueVec, (fit_tools.locationMat[:, 0], fit_tools.locationMat[:, 1])),
        #              shape=(fit_tools.dimNum, fit_tools.dimNum)).tocsc(),
        'sparseGRM': sp.coo_matrix((ValueVec, (LocationMat[0, :], LocationMat[1, :])),
                                   shape=(dimNum, dimNum)).tocsc(),
        'read_in_tools': read_in_tools,
    }

    return glmm_result




    #
    # '''Optimzie ML for the fix-coefficients alpha and random effects b'''
    # X=np.ones((dimNum,1))
    # fit_tools = fitGLMM_tools(ValueVec, LocationMat, dimNum)
    # re_coef = fit_tools.get_coef(y, X, tau, family_indicator, alpha0, eta0, offset, verbose=verbose, maxiterPCG=maxiterPCG, tolPCG=tolPCG, maxiter=maxiter)
    #
    #
    #
    # if family_indicator in ["Poisson", "Binomial"]:
    #     '''AI_REML for variance components (phi, tau)'''
    #     nrun = 10
    #     re = fit_tools.getAIScore(family_indicator, re_coef['Y'], X, re_coef['W'], tau, re_coef['Sigma_iY'],
    #                               re_coef['Sigma_iX'], re_coef['cov'],
    #                               nrun, maxiterPCG, tolPCG=tolPCG, traceCVcutoff=traceCVcutoff)
    #     print("AI")
    #     print(re['AI'])
    #
    #     tau[1] = max(0, tau0[1] + tau0[1]**2 * (re['YPAPY'] - re['Trace']) / N)
    #
    #     if verbose:
    #         print("Variance component estimates:", tau)
    #
    #     for i in range(maxiter):
    #
    #         alpha0 = re_coef['alpha']
    #         tau0 = tau
    #         print("tau0_v1:", tau0)
    #         eta0 = eta
    #
    #         # use Get_Coef before getAIScore
    #         t_begin_get_coef = pd.to_datetime("today")
    #         # Get stable coefficients
    #         re_coef = fit_tools.get_coef(y, X, tau, family_indicator, alpha0, eta0, offset, verbose=verbose, maxiterPCG=maxiterPCG,
    #                            tolPCG=tolPCG, maxiter=maxiter)
    #         t_end_get_coef = pd.to_datetime("today")
    #         print("t_end_get_coef - t_begin_get_coef:")
    #         print(t_end_get_coef - t_begin_get_coef)
    #
    #         fit = fit_tools.fit_glmmaiRPCG(family_indicator,re_coef['Y'], X, re_coef['W'], tau, re_coef['Sigma_iY'], re_coef['Sigma_iX'],
    #                              re_coef['cov'], nrun, maxiterPCG, tolPCG, tol=tol, traceCVcutoff=traceCVcutoff)
    #
    #         t_end_fit_glmmaiRPCG = pd.to_datetime("today")
    #         print("t_end_fit_glmmaiRPCG - t_end_get_coef")
    #         print(t_end_fit_glmmaiRPCG - t_end_get_coef)
    #
    #         tau = [float(i) for i in fit['tau']]
    #         tau = np.array(tau)
    #         tau0 = np.array(tau0)
    #         cov = re_coef['cov']
    #         alpha = re_coef['alpha']
    #         eta = re_coef['eta']
    #         Y = re_coef['Y']
    #         mu = re_coef['mu']
    #
    #         print(np.abs(tau - tau0) / (np.abs(tau) + np.abs(tau0) + tol))
    #
    #         if verbose:
    #             print("tau:", tau)
    #             print("tau0:", tau0)
    #
    #         if tau[1] == 0:
    #             break
    #
    #         # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
    #         if np.max(np.abs(tau - tau0) / (np.abs(tau) + np.abs(tau0) + tol)) < tol:
    #             break
    #
    #         if np.max(tau) > tol**(-2):
    #             Warning("Large variance estimate observed in the iterations, model not converged...", call=False)
    #             i = maxiter
    #             break
    #
    #     if verbose:
    #         print("Final tau:", tau)
    #
    #     # added these steps after tau is estimated
    #     re_coef = fit_tools.get_coef(y, X, tau, family_indicator, alpha, eta, offset, verbose=verbose, maxiterPCG=maxiterPCG,
    #                        tolPCG=tolPCG, maxiter=maxiter)
    #     cov = re_coef['cov']
    #     alpha = re_coef['alpha']
    #     eta = re_coef['eta']
    #     Y = re_coef['Y']
    #     mu = re_coef['mu']
    #
    #     converged = i < maxiter
    #     res = y - mu
    #     coef_alpha = alpha
    #     mu2 = mu * (1 - mu)
    #
    #     glmm_result = {
    #         'tau': tau,
    #         'coefficients': coef_alpha,
    #         'linear.predictors': eta,
    #         'fitted.values': mu,
    #         'Y': Y,
    #         'residuals': res,
    #         'cov': cov,
    #         'converged': converged,
    #         'y': y,
    #         'X': X,
    #         'traitType': "binary",
    #         # 'sparseGRM': sp.coo_matrix((fit_tools.valueVec, (fit_tools.locationMat[:, 0], fit_tools.locationMat[:, 1])),
    #         #              shape=(fit_tools.dimNum, fit_tools.dimNum)).tocsc(),
    #         'sparseGRM': sp.coo_matrix((fit_tools.valueVec, (fit_tools.locationMat[0, :], fit_tools.locationMat[1, :])),
    #                                    shape=(fit_tools.dimNum, fit_tools.dimNum)).tocsc(),
    #         'read_in_tools': read_in_tools,
    #     }
    #
    #     t_end_null = time.time()
    #     print("t_end_null - t_begin, fitting the NULL model without LOCO took")
    #     print(t_end_null - t_begin)
    #
    #     return glmm_result
    #
    # elif family_indicator == "Beta":
    #
    #
    #     '''AI_REML for variance components (phi, tau)'''
    #     nrun = 10
    #     re = fit_tools.getAIScore(family_indicator, re_coef['Y'], X, re_coef['W'], tau, re_coef['Sigma_iY'],
    #                               re_coef['Sigma_iX'], re_coef['cov'],
    #                               nrun, maxiterPCG, tolPCG=tolPCG, traceCVcutoff=traceCVcutoff)
    #     print("AI")
    #     print(re['AI'])
    #     tau[1] = max(0, tau0[1] + tau0[1] ** 2 * (re['YPAPY'] - re['Trace'][1]) / N)
    #     tau[0] = max(0, tau0[0] + tau0[0] ** 2 * (re['YPAPY'] - re['Trace'][0]) / N)
    #
    #     if verbose:
    #         print("Variance component estimates:", tau)
    #     for i in range(maxiter):
    #         if verbose:
    #             # print(f"\nIteration {i + 1}:\n")
    #             print("\n", "Iteration :", i + 1, "\n")
    #         alpha0 = re_coef['alpha']
    #         tau0 = tau
    #         eta0 = eta
    #         re_coef = fit_tools.get_coef(y, X, tau, family_indicator, alpha0, eta0, offset, verbose=verbose,
    #                                      maxiterPCG=maxiterPCG,
    #                                      tolPCG=tolPCG, maxiter=maxiter)
    #
    #         fit = fit_tools.fitglmmaiRPCG_q(family_indicator, re_coef['Y'], X, re_coef['W'], tau, re_coef['Sigma_iY'],
    #                                         re_coef['Sigma_iX'],
    #                                         re_coef['cov'], nrun, maxiterPCG, tolPCG, tol=tol,
    #                                         traceCVcutoff=traceCVcutoff)
    #
    #         tau = np.array(fit['tau'])
    #         cov = re_coef['cov']
    #         alpha = re_coef['alpha']
    #         eta = re_coef['eta']
    #
    #         print("cov: ", cov)
    #
    #         if verbose:
    #             print("Variance component estimates:")
    #             print(tau)
    #             print("Fixed-effect coefficients:")
    #             print(alpha)
    #
    #         Y = re_coef['Y']
    #         mu = re_coef['mu']
    #
    #         if tau[0] <= 0:
    #             raise Exception("ERROR! The first variance component parameter estimate is 0")
    #
    #         if tau[1] <= 0: break
    #
    #         if np.max(np.abs(tau - tau0) / (np.abs(tau) + np.abs(tau0) + tol)) < tol: break
    #
    #         if np.max(tau) > tol ** (-2):
    #             print("Large variance estimate observed in the iterations, model not converged...")
    #             i = maxiter
    #             break
    #
    #     if verbose:
    #         print("Final tau:", tau)
    #
    #     t_end_null = time.time()
    #     print("t_end_null - t_begin, fitting the NULL model without LOCO took")
    #     print(t_end_null - t_begin)
    #
    #     re_coef = fit_tools.get_coef(y, X, tau, family_indicator, alpha, eta, offset, verbose=verbose,
    #                                  maxiterPCG=maxiterPCG,
    #                                  tolPCG=tolPCG, maxiter=maxiter)
    #     cov = re_coef['cov']
    #     alpha = re_coef['alpha']
    #     eta = re_coef['eta']
    #     Y = re_coef['Y']
    #     mu = re_coef['mu']
    #
    #     converged = True if i < maxiter else False
    #     res = y - mu
    #     mu2 = np.repeat((1 / tau[1]), len(res))
    #     """data 没定义
    #     out_transform = Covariate_Transform (data)
    #     # out.transform <- Covariate_Transform(formula.null, data = dataMerge_sort)
    #     if out_transform is not None:
    #         coef_alpha = Covariate_Transform_Back(alpha, covariate_transform['Param_transform'])
    #     else:
    #         coef_alpha = alpha"""
    #
    #     """先这样写，跑通了再说"""
    #     coef_alpha = alpha
    #     # obj_noK = ScoreTest_NULL_Model(mu, mu2, y, X)
    #     # lmmResult = {
    #     glmm_result = {
    #         # 'theta': tau,
    #         'tau': tau,
    #         'coefficients': coef_alpha,
    #         'linear.predictors': eta,
    #         'fitted.values': mu,
    #         'Y': Y,
    #         'residuals': res,
    #         'cov': cov,
    #         'converged': converged,
    #         # 'obj.noK': obj_noK,
    #         'y': y,
    #         'X': X,
    #         'traitType': "quantitative",
    #         'read_in_tools': read_in_tools,
    #         # 'sparseGRM': sp.coo_matrix((fit_tools.valueVec, (fit_tools.locationMat[:, 0], fit_tools.locationMat[:, 1])),
    #         #                            shape=(fit_tools.dimNum, fit_tools.dimNum)).tocsc(),
    #         'sparseGRM': sp.coo_matrix((fit_tools.valueVec, (fit_tools.locationMat[0, :], fit_tools.locationMat[1, :])),
    #                                    shape=(fit_tools.dimNum, fit_tools.dimNum)).tocsc(),
    #     }
    #
    #     t_end_null = time.time()
    #     print("t_end_null - t_begin, fitting the NULL model without LOCO took")
    #     print(t_end_null - t_begin)
    #     return glmm_result

if __name__ == "__main__":
    Geno_name = '21'
    Cell_name = 'CD4+_T_naive'
    pheno_index = 1 # 哪一条phenotype的指标，跑实验的话需要循环起来
    geno_index = 2 # 哪一条genotype来做spa test

    GLMM_Result = Step1_FitGLMM(prefix = sys.path[0] + "/Experiment_data/",
                                Geno = Geno_name,
                                Cell = Cell_name,
                                index = pheno_index, 
                                M = 20, 
                                family_indicator = 'Binomial', 
                                Tau_dic = {"tauInit": [0, 0], "tau": [0, 0], "fixtau": [0, 0]}, 
                                maxiterPCG=3, tolPCG=1e-1, maxiter=3, tol=1e-1, traceCVcutoff=1e-1,
                                verbose=True)

    # 优化结束后，打印优化结果
    print('tau', GLMM_Result['tau'])
    print('Coefficients', GLMM_Result['coefficients'])
    np.savetxt(sys.path[0] + "/fitted_values.csv", GLMM_Result['fitted.values'], delimiter="," )
