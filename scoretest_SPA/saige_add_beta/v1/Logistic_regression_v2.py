from SAIGE.Step1_FittGLMM import fit_glmma
import pandas as pd
import sys



# fit_glmma 函数优化模型
glmm_result = fit_glmma(prefix = sys.path[0] + "/sample_data/nfam_100_nindep_0_step1_includeMoreRareVariants_poly", 
                        Tau_dic = {"tauInit": [0, 0], "tau": [0, 0], "fixtau": [0, 0]}, 
                        verbose=True, 
                        isDiagofKinSetAsOne=True, 
                        maxiterPCG=100, tolPCG=1e-6, maxiter=100, tol=1e-4, traceCVcutoff=1e-4,
                        LOCO=False)

# 优化结束后，打印优化结果
print('theta', glmm_result['theta'])
print('Coefficients', glmm_result['coefficients'])
