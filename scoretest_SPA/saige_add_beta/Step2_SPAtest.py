import numpy as np

# from SAIGE_FitGLMM_tools import fitGLMM_tools
# from SAIGE_DataPreparation_tools import read_in_Genotype
from SAIGE_DataPreparation_tools import *
from SAIGE_SPACalculation_tools import *













def Step2_SPAtest(glmm_result, pheno_index, geno_index):

    # 0. Prepare data
    read_in_tools = glmm_result['read_in_tools']
    G = read_in_tools.Specific_Genotype(index = geno_index).reshape(-1, 1)
    X = glmm_result['X']

    y = glmm_result['y']
    alpha = glmm_result['coefficients']
    print("alpha",alpha)

    mu = glmm_result['fitted.values']
    # print("mu",mu)
    Psi = glmm_result['sparseGRM']
    V = mu * (1 - mu)
    V = V.reshape(-1, 1)
    tau = glmm_result['tau'][1]


    # 保留每个元素的第一个值
    first_column = X[:, 0:1]
    alpha_column = np.full_like(first_column, alpha)
    X = np.hstack((alpha_column, X))

    # 用于消除共线性
    # X1 = ScoreTest_wSaddleApprox_Get_X1(X)



    XV = (X * V).T
    XVX_inv = np.linalg.inv(X.T @ (X * V))
    XXVX_inv = X @ XVX_inv


    n_g = np.sum(G)
    if n_g / (2 * len(G)) > 0.5:
        G = 2 - G
        n_g = np.sum(G)

    # NAset = np.where(G == 0)[0]  # 获取 G 中为 0 的元素索引

    G1 = G - XXVX_inv @ (XV @ G)
    # q = np.sum(G1 * y)
    q = np.dot(G1.flatten(), y)




    # 1. Compute the constant ratio r
    r,c = Estimation_r(Psi, V, tau)

    # 2. Calculate the Score Statistics T
    # T_adj,g = Estimation_T_adj(G, X, V, y, mu, r)

    # 3. Use SPA caculate p-value
    # spa_model = SPA(g,mu,r,c,V,y,T_adj)
    spa_model = SPA(G1, mu, r, c, V, y, q)
    pval_noadjust,qinv,var1,m1 = spa_model.calculate_pval_noadjust()
    p_value = spa_model.calculate_F(pval_noadjust,qinv,var1,m1,2)

    return pval_noadjust,p_value
#
#
# if __name__ == "__main__":
#     import sys
#     from scipy import sparse as sp
#     from SAIGE_DataPreparation_tools import read_and_create_GRM
#
#
#     prefix = sys.path[0] + "/sample_data/nfam_100_nindep_0_step1_includeMoreRareVariants_poly"
#
#     '''prepare input data'''
#     Geno_name = '21'
#     Cell_name = 'CD4+_T_naive'
#     geno_index = 2
#
#     read_in_tools = Read_in_tools(sys.path[0] + "/Experiment_data/", Cell_name, Geno_name)
#     G = read_in_tools.Specific_Genotype(index = geno_index).reshape(-1, 1)
#     Geno =  read_in_tools.Create_GRM(M = 20)
#     ValueVec = Geno.data
#     LocationMat = np.concatenate((Geno.row.reshape(-1, 1), Geno.col.reshape(-1, 1)), axis = 1)
#     grm = sp.coo_matrix((ValueVec, (LocationMat[:, 0], LocationMat[:, 1])),
#                      shape=(1000, 1000)).tocsc()
#
#     mu = np.loadtxt(sys.path[0] + "/sample_data/spa_tem/fitxted_values.csv", delimiter=",")
#     glmm_result = {
#         'tau': np.array([1.0, 4.59353571]),
#         'coefficients': np.array([-0.25676325, -0.15797047]),
#         'mu': mu,
#         'sparseGRM': grm
#     }
#
#     p_value = Step2_SPAtest(prefix, glmm_result)
#     print('p-value is:', p_value)