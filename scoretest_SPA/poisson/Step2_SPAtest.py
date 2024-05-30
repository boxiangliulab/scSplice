import numpy as np

# from SAIGE_FitGLMM_tools import fitGLMM_tools
# from SAIGE_DataPreparation_tools import read_in_Genotype
from SAIGE_DataPreparation_tools import *
from SAIGE_SPACalculation_tools import *


# def Step2_SPAtest(glmm_result, pheno_index, geno_index):
#
#     # 0. Prepare data
#     read_in_tools = glmm_result['read_in_tools']
#     # G = read_in_tools.Specific_Genotype(index = geno_index).reshape(-1, 1)
#     N = 10000
#     maf = 0.05
#     num_ones = int(N * maf)
#     num_zeros = N - num_ones
#     G = np.array([1] * num_ones + [0] * num_zeros)
#     np.random.shuffle(G)
#     G = G.reshape(-1, 1)
#     print("G",G)
#
#
#     X = glmm_result['X']
#     y = glmm_result['y']
#     # print("Correlation between y and G:", np.corrcoef(y, G[:, 0])[0, 1])
#     alpha = glmm_result['coefficients']
#     print("alpha",alpha)
#
#     mu = glmm_result['fitted.values']
#     # print("mu",mu)
#     Psi = glmm_result['sparseGRM']
#     V = mu * (1 - mu)
#     V = V.reshape(-1, 1)
#     tau = glmm_result['tau'][1]
#
#
#     # 保留每个元素的第一个值
#     first_column = X[:, 0:1]
#     alpha_column = np.full_like(first_column, alpha)
#     X = np.hstack((alpha_column, X))
#
#     # 用于消除共线性
#     # X1 = ScoreTest_wSaddleApprox_Get_X1(X)
#
#
#
#     XV = (X * V).T
#     XVX_inv = np.linalg.inv(X.T @ (X * V))
#     XXVX_inv = X @ XVX_inv
#
#
#     n_g = np.sum(G)
#     if n_g / (2 * len(G)) > 0.5:
#         G = 2 - G
#         n_g = np.sum(G)
#
#     # NAset = np.where(G == 0)[0]  # 获取 G 中为 0 的元素索引
#
#     G1 = G - XXVX_inv @ (XV @ G)
#     # q = np.sum(G1 * y)
#     q = np.dot(G1.flatten(), y)
#
#
#
#
#     # 1. Compute the constant ratio r
#     r,c = Estimation_r(Psi, V, tau)
#
#     # 2. Calculate the Score Statistics T
#     # T_adj,g = Estimation_T_adj(G, X, V, y, mu, r)
#
#     # 3. Use SPA caculate p-value
#     # spa_model = SPA(g,mu,r,c,V,y,T_adj)
#     spa_model = SPA(G1, mu, r, c, V, y, q)
#     pval_noadjust,qinv,var1,m1 = spa_model.calculate_pval_noadjust()
#     p_value = spa_model.calculate_F(pval_noadjust,qinv,var1,m1,2)
#
#     return pval_noadjust,p_value




import numpy as np

# from SAIGE_FitGLMM_tools import fitGLMM_tools
# from SAIGE_DataPreparation_tools import read_in_Genotype
from SAIGE_DataPreparation_tools import *
from SAIGE_SPACalculation_tools import *


def Step2_SPAtest(glmm_result, pheno_index, geno_index):

    # 0. Prepare data
    read_in_tools = glmm_result['read_in_tools']
    # G = read_in_tools.Specific_Genotype(index = geno_index).reshape(-1, 1)
    N = 10000
    maf = 0.05
    num_ones = int(N * maf)
    num_zeros = N - num_ones
    G = np.array([1] * num_ones + [0] * num_zeros)
    np.random.shuffle(G)
    G = G.reshape(-1, 1)
    print("G",G)


    X = glmm_result['X']
    y = glmm_result['y']
    # print("Correlation between y and G:", np.corrcoef(y, G[:, 0])[0, 1])
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
    XVG0 = np.dot(XV, G)

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




