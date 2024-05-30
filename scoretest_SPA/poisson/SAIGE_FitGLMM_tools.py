import sys
import time

from scipy.linalg import inv, pinv
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix, csr_matrix
import scipy.sparse as sp
from scipy.special import expit
from numpy.random import binomial
from numpy import dot
from sklearn.linear_model import LogisticRegression, PoissonRegressor
from sklearn.linear_model import LinearRegression



import numpy as np
import pandas as pd
from rpy2.robjects import r, Formula
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.linear_model import LogisticRegression

import statsmodels.api as sm
from statsmodels.formula.api import ols
import patsy
import os

class fitGLMM_tools():
    def __init__(self, valueVec, locationMat, dimNum):
        self.locationMat = locationMat
        self.valueVec = valueVec  
        self.dimNum = dimNum

        # 占位符 需确定具体含义
        self.isUseSparseSigmaforInitTau = False

        self.isUseSparseSigmaforModelFitting = True
        #self.isUseSparseSigmaforModelFitting = False

        self.isDiagofKinSetAsOne = False

        self.isUsePrecondM = True
        #self.isUsePrecondM = False
        

    def getCrossprod(self, bVec, wVec, tauVec):
        crossProdVec = np.zeros_like(bVec, dtype=np.float32)

        # Added by SLEE, 04/16/2017
        if tauVec[1] == 0:
            crossProdVec = tauVec[0] * (bVec * (1 / wVec))
            return crossProdVec

        # Assuming getCrossprodMatAndKin is a function that returns crossProd1
        crossProd1 = self.getCrossprodMatAndKin(bVec)

        # Replace the loop with vectorized operations
        crossProdVec = tauVec[0] * (bVec * (1 / wVec)) + tauVec[1] * crossProd1

        return crossProdVec



    def gen_sp_Sigma(self, wVec, tauVec):
        dtVec = (1 / wVec) * tauVec[0]
        
        valueVecNew = self.valueVec * tauVec[1]
        nnonzero = self.valueVec.size
        # for i in range(nnonzero):
        #     if self.locationMat[i, 0] == self.locationMat[i, 1]:
        #         valueVecNew[i] += dtVec[self.locationMat[i, 0]]
        #         if valueVecNew[i] < 1e-4:
        #             valueVecNew[i] = 1e-4
        #
        # # Create sparse matrix
        # result = sp.coo_matrix((valueVecNew, (self.locationMat[:, 0], self.locationMat[:, 1])), shape=(self.dimNum, self.dimNum)).tocsc()

        #改成2xn以后的
        for i in range(nnonzero):
            if self.locationMat[0, i] == self.locationMat[1, i]:
                valueVecNew[i] += dtVec[self.locationMat[0, i]]
                valueVecNew[i] = max(valueVecNew[i], 1e-4)

        result = coo_matrix((valueVecNew, (self.locationMat[0], self.locationMat[1])), shape=(self.dimNum, self.dimNum)).tocsc()

        return result

    def gen_spsolve_v4(self, wVec, tauVec, yVec):


        # Convert input to appropriate data types
        yVec2 = np.array(yVec, dtype=float)

        # Call the function to generate sparse Sigma matrix
        result = self.gen_sp_Sigma(wVec, tauVec)

        if result.shape[0] != result.shape[1] or (result.diagonal() == 0).any():
            raise ValueError("Sigma matrix is singular or not square, cannot solve")

        # Solve the sparse linear system
        x = spsolve(result, yVec2)

        # Convert the result back to the original data type
        z = x.astype(np.float32)
        return z
    
    def spsolove_for_Tadj(self, XWX, XVec):
        sXWX = sp.coo_matrix(XWX).tocsc()
        # Solve the sparse linear system
        x = spsolve(sXWX, XVec)
        # z = np.array(x, dtype=float)
        return x

    def getDiagOfSigma(self, wVec, tauVec):
        Nnomissing = self.dimNum

        #M = self.geno.getM() #没用上
        #MminMAF = self.geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM()
        MminMAF = 0.01
        diagVec = np.zeros(Nnomissing, dtype=float)

        setKinDiagtoOne = True  # GRM对角线是0还是1
        """geno.Get_Diagof_StdGeno()返回了一个指针，我好想死"""
        if not setKinDiagtoOne:
            diagVec = tauVec[1] * (self.geno.Get_Diagof_StdGeno() / MminMAF) + tauVec[0] / wVec
        else:
            diagVec = tauVec[1] + tauVec[0] / wVec

        for i in range(Nnomissing):
            if diagVec[i] < 1e-4:
                diagVec[i] = 1e-4

        return diagVec

    def getCrossprodMatAndKin(self, bVec):
        '''
        Build a sparse matrix for ValueVec (attributes in class) and then operate the dot multiply between Valuevec and bVec

        Args:
            bVec (array): Multiply matrixx, (n, ).

        Returns:
            crossProdVec (array): ValueVec * bVec.
        '''

        # self.isUseSparseSigmaforModelFitting = True
        if self.isUseSparseSigmaforInitTau or self.isUseSparseSigmaforModelFitting:
            # # print("use sparse kinship to estimate initial tau and for getCrossprodMatAndKin")
            # result = coo_matrix((self.valueVec, (self.locationMat[:, 0], self.locationMat[:, 1])), shape=(self.dimNum, self.dimNum))
            # x = result.dot(bVec)
            # crossProdVec = np.array(x).flatten().astype(float)


            rows = self.locationMat[0, :]  # 第一行是行索引
            cols = self.locationMat[1, :]  # 第二行是列索引
            # 构建稀疏矩阵
            result = coo_matrix((self.valueVec, (rows, cols)), shape=(self.dimNum, self.dimNum))
            # 矩阵向量乘法
            x = result.dot(bVec)  # 确保 b_vec 是一个列向量或一维数组
            crossProdVec = x.astype(np.float32)





        else:
            # crossProdVec = self.parallelCrossProd(bVec)
            pass
        return crossProdVec
    
    def getCrossprodMatAndKin_Tadj(self, bVec):
        '''
        Build a sparse matrix for ValueVec (attributes in class) and then operate the dot multiply between Valuevec and bVec

        Args:
            bVec (array): Multiply matrixx, (n, ).

        Returns:
            crossProdVec (array): ValueVec * bVec.
        '''

        self.isUseSparseSigmaforModelFitting = True
        if self.isUseSparseSigmaforInitTau or self.isUseSparseSigmaforModelFitting:
            # print("use sparse kinship to estimate initial tau and for getCrossprodMatAndKin")
            # print("---------------")
            # print("self.locationMat",self.locationMat)
            # print("self.valueVec",self.valueVec)

            result = coo_matrix((self.valueVec, (self.locationMat[:, 0], self.locationMat[:, 1])), shape=(self.dimNum, self.dimNum))
            # print("result",result)
            x = result.dot(bVec)
            # crossProdVec = np.array(x).flatten().astype(float)
        else:
            # crossProdVec = self.parallelCrossProd(bVec)
            pass
        return x
    
    def linkinv(self, eta0):
        if self.family == 'Binomial':
            tem = expit(eta0)
        else:
            tem = np.exp(eta0)
        return tem
    
    def mu_eta(self, eta0):
        if self.family == 'Binomial':
            tem = np.exp(eta0) / (1 + np.exp(eta0))**2
        else:
            tem = 1 / eta0
        return tem
    
    def getPCG1ofSigmaAndVector(self, wVec, tauVec, YVec, maxiterPCG, tolPCG):
        # Implement the PCG method to get Sigma for a given vector

        Nnomissing = self.dimNum
        xVec = np.zeros(Nnomissing, dtype=np.float32)


        # Check if sparse Sigma should be used for initialization of tau
        if self.isUseSparseSigmaforInitTau:
            print("use sparse kinship to estimate initial tau")
            xVec = self.gen_spsolve_v4(wVec, tauVec, YVec)
        # Check if sparse Sigma should be used for model fitting
        elif self.isUseSparseSigmaforModelFitting:
            print("use sparse kinship to fit the model")
            xVec = self.gen_spsolve_v4(wVec, tauVec, YVec)
        else:
            rVec =  YVec.copy()
            # xVec = np.zeros(Nnomissing, dtype=np.float32)
            # crossProdVec = np.zeros(Nnomissing, dtype=np.float32)
            # zVec = np.zeros(Nnomissing, dtype=np.float32)
            minvVec = np.zeros(Nnomissing, dtype=np.float32)

            # Check if precondM should be used
            if not self.isUsePrecondM: # 不使用precondition M
                print("using inverse of diagonal of Sigma for preconditioning")
                minvVec = 1 / self.getDiagOfSigma(wVec, tauVec)
                zVec = minvVec * rVec
            else: # 使用precondition M
                # zVec = spsolve(wVec, rVec)
                zVec = self.gen_spsolve_v4(wVec, tauVec, rVec)

            iter = 0
            sumr2 = np.sum(rVec ** 2)

            pVec =  zVec.copy()
            while sumr2 > tolPCG and iter < maxiterPCG:
                iter += 1

                # Perform PCG iteration
                # preA = np.dot(rVec, zVec) / np.dot(pVec, ApVec)
                # a = preA.item()
                ApVec = self.getCrossprod(pVec, wVec, tauVec)
                a = np.dot(rVec.T, zVec) / np.dot(pVec.T, ApVec)

                xVec += a * pVec
                r1Vec = rVec - a * ApVec
                # z1Vec = np.zeros(Nnomissing, dtype=np.float32)

                # Check if precondM should be used
                isUsePrecondM = True
                if not isUsePrecondM:
                    z1Vec = minvVec * r1Vec
                else:
                    # z1Vec = spsolve(wVec, r1Vec)
                    z1Vec = self.gen_spsolve_v4(wVec, tauVec, r1Vec)

                # Prebet = np.dot(z1Vec, r1Vec) / np.dot(zVec, rVec)
                # bet = Prebet.item()

                bet = np.dot(z1Vec.T, r1Vec) / np.dot(zVec.T, rVec)
                pVec = z1Vec + bet * pVec
                zVec = z1Vec
                rVec = r1Vec

                sumr2 = np.sum(rVec ** 2)

            if iter >= maxiterPCG:
                print("pcg did not converge. You may increase maxiter number.")
            print("iter from get_PCG1_of_Sigma_And_Vector", iter)

        return xVec


    def getCoefficients(self, Yvec, Xmat, wVec, tauVec, maxiterPCG, tolPCG):
        '''
        Use PCG methods to calculate several intermediate terms.

        Args:
            y (array): The response variables, [n*1, Binary];
            X (array): The input variables, Covariates and Genotypes, [n*p, Continuous];
            tau (list): List for initial values for the variance component parameter estimates;
            maxiterPCG (int): The max interation time for PCG;
            tolPCG (float): The tolerance of PCG algorithm for whether stop iterations.

        Returns:
            Sigma_iY (array): The intermediate terms, (n, );
            Sigma_iX (array): The intermediate terms, (n, );
            cov (array): The covariance matrix of inputs vairables X, (p, p);
            alpha (array): Fixed-effect coefficients, (p, );
            eta (array): The linear predictor, (n, );
        '''


        Nnomissing = self.dimNum

        # Calculate Sigma_iY using PCG method
        Sigma_iY = self.getPCG1ofSigmaAndVector(wVec, tauVec, Yvec, maxiterPCG, tolPCG)

        # Get the number of columns in Xmat
        colNumX = Xmat.shape[1]

        # Initialize Sigma_iX matrix
        Sigma_iX = np.zeros((Nnomissing, colNumX), dtype=np.float32)

        # Calculate Sigma_iX for each column of Xmat using PCG method
        for i in range(colNumX):
            XmatVecTemp = Xmat[:, i]
            Sigma_iX[:, i] = self.getPCG1ofSigmaAndVector(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG)

        # Transpose Xmat
        Xmatt = Xmat.T

        # Calculate the inverse or pseudo-inverse of (Xmatt * Sigma_iX)
        try:
            cov = inv(np.dot(Xmatt, Sigma_iX))
        except np.linalg.LinAlgError:
            cov = pinv(np.dot(Xmatt, Sigma_iX))
            print("inv_sympd failed, inverted with pinv")

        # Transpose Sigma_iX
        Sigma_iXt = Sigma_iX.T

        # Calculate Sigma_iXtY
        SigmaiXtY = np.dot(Sigma_iXt, Yvec)

        # Calculate alpha
        alpha = np.dot(cov, SigmaiXtY)

        # Calculate eta
        eta = Yvec - tauVec[0] * (Sigma_iY - np.dot(Sigma_iX, alpha)) / wVec

        # Return the results as a dictionary
        return {"Sigma_iY": Sigma_iY, "Sigma_iX": Sigma_iX, "cov": cov, "alpha": alpha, "eta": eta }



    def get_coef(self, y, X, tau, family, alpha0, eta0, offset, maxiterPCG, tolPCG, maxiter, verbose=True):
        '''
        Get Coefficients for the GLMM, logit(mu) = alpha X + beta G + b

        Args:
            y (array): The response variables, [n*1, Binary];
            X (array): The input variables, Covariates and Genotypes, [n*p, Continuous];
            tau (list): List for initial values for the variance component parameter estimates;
            alpha0 (array): alpha_0 is the null model parameter, [(p, )];
            eta0 (array): eta_0 is the linear predictor, alpha X + beta G, [(n, )];
            offset (None/array): offset for the null GLMM model, [(n, )];
            maxiterPCG (int): The max interation time for PCG;
            tolPCG (float): The tolerance of PCG algorithm for whether stop iterations.

        Returns:
            Y (array): The response variables, (n, 1);
            alpha (array): Fixed-effect coefficients, (p, );
            eta (array): The linear predictor, (n, );
            W (array): The variance components mu(1 - mu), (n, );
            cov (array): The covariance matrix of inputs vairables X, X^T*Sigma_i* X, (p, p);
            sqrtW (array): The sqrt variance components sqrt(mu(1 - mu)), (n, );
            Sigma_iY (array): The intermediate terms, (n, );
            Sigma_iX (array): The intermediate terms, (n, );
            mu (array): The probability of model, (n, );
        '''
        self.family = family
        tol_coef = 0.1 
        mu = self.linkinv(eta0) # Linear predictor 通过反向的变换函数到均值mu

        mu_eta = self.mu_eta(eta0) # 计算mu对eta求导的导数值

        Y = eta0 - offset + (y.flatten() - mu) / mu_eta # 计算tilde{Y} 就是经过修饰以后的Y

        sqrtW = mu_eta / np.sqrt(np.var(mu)) 
        W = sqrtW ** 2 

        for i in range(maxiterPCG):
            # Use PCG to update intermediate terms
            re_coef = self.getCoefficients(Y, X, W, tau, maxiterPCG=maxiterPCG, tolPCG=tolPCG)
            alpha = re_coef['alpha']
            eta = re_coef['eta'] + offset

            if verbose:
                print("Iteration times:", i)
                print("Tau:", tau)
                print("Fixed-effect coefficients:", alpha)

            mu = self.linkinv(eta)
            mu_eta = self.mu_eta(eta)

            Y = eta - offset + (y.flatten() - mu) / mu_eta
            sqrtW = mu_eta / np.sqrt(np.var(mu))
            W = sqrtW ** 2

            if np.max(np.abs((alpha - alpha0) / (np.abs(alpha) + np.abs(alpha0) + tol_coef))) < tol_coef:
                break

            alpha0 = alpha

        return {'Y': Y, 'alpha': alpha, 'eta': eta, 'W': W, 'cov': re_coef['cov'], 'sqrtW': sqrtW, 'Sigma_iY': re_coef['Sigma_iY'], 'Sigma_iX': re_coef['Sigma_iX'], 'mu': mu}


    def fit_glmmaiRPCG(self, family_indicator,Yvec, Xmat, wVec, tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, tol, traceCVcutoff):
        re = self.getAIScore(family_indicator,Yvec, Xmat, wVec, tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff)
        YPAPY = re["YPAPY"]
        Trace = re["Trace"]
        score1 = YPAPY - Trace
        AI1 = re["AI"]
        Dtau = score1 / AI1
        tau0 = np.copy(tauVec)
        tauVec[1] = tau0[1] + Dtau

        for i in range(len(tauVec)):
            if tauVec[i] < tol:
                tauVec[i] = 0

        step = 1.0
        while tauVec[1] < 0.0:
            step *= 0.5
            tauVec[1] = tau0[1] + step * Dtau

        for i in range(len(tauVec)):
            if tauVec[i] < tol:
                tauVec[i] = 0

        return {"tau": tauVec}

    def fitglmmaiRPCG_q(self, family_indicator, Yvec, Xmat, wVec, tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG,
                        tolPCG, tol,
                        traceCVcutoff):
        tauVec = np.array(tauVec)
        zeroVec = (tauVec < tol)  # 对应 arma::uvec zeroVec = (tauVec < tol);

        re = self.getAIScore(family_indicator, Yvec, Xmat, wVec, tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG,
                             tolPCG, traceCVcutoff)

        YPAPY = re["YPAPY"]
        YPA0PY = re["YPA0PY"]
        Trace = re["Trace"]

        score0 = YPA0PY - Trace[0]
        score1 = YPAPY - Trace[1]
        scoreVec = np.array([score0, score1])

        AI = re["AI"]
        # Dtau = np.linalg.solve(AI, scoreVec)  # 解线性方程组 AI是奇异矩阵 不可逆
        Dtau = np.dot(np.linalg.pinv(AI), scoreVec)
        print("Dtau===", Dtau)
        tau0 = tauVec.copy()
        tauVec = tau0 + Dtau

        # 条件更新tauVec的元素为0，类似于C++中的elem(find()).zeros()操作
        condition = zeroVec & (tauVec < tol)
        tauVec[condition] = 0.0

        step = 1.0
        while tauVec[0] < 0.0 or tauVec[1] < 0.0:
            step = step * 0.5
            tauVec = tau0 + step * Dtau

            condition = zeroVec & (tauVec < tol)
            tauVec[condition] = 0.0

        # 最后一次更新小于tol的元素为0
        tauVec[tauVec < tol] = 0.0

        return {"tau": tauVec}

    def calCV(self, xVec):
        veclen = len(xVec)
        vecMean = np.mean(xVec)
        vecSd = np.std(xVec)
        vecCV = (vecSd / vecMean) / veclen
        return vecCV

    def GetTrace(self, Sigma_iX, Xmat, wVec, tauVec, cov1, nrun, maxiterPCG, tolPCG, traceCVcutoff):
        '''
        Get Trace of which PW_i and Pphi by using Randomized trace estimator

        Args:
            Sigma_iX (array): The intermediate terms, (n, );
            Xmat (array): The input variables, Covariates and Genotypes, [n*p, Continuous];
            wVec (array): The variance components mu(1 - mu), (n, );
            tauVec (list): List for initial values for the variance component parameter estimates;
            cov1 (array): The covariance matrix of inputs vairables X, (p, p);
            nrun (int): iteraion times of ?;
            maxiterPCG (int): The max interation time for PCG;
            tolPCG (float): The tolerance of PCG algorithm for whether stop iterations;
            traceCVcutoff (float): Threshold for the coefficient of variation for trace estimation.

        Returns:
            YPAPY (array): The matrix YPAPY, (n, 1);
            Trace (float): The trace of ?;
            PY1 (array): The matrix P in paper multiply Y, ();
            AI (array): The matrix YPAPAPY AI_{\tau\tau}, ()
        '''
        np.random.seed(200)

        Nnomissing = self.dimNum
        Sigma_iXt = Sigma_iX.T

        nrunStart = 0
        nrunEnd = nrun

        traceCV = traceCVcutoff + 0.1
        tempVec = np.zeros(nrun)

        while traceCV > traceCVcutoff:
            for i in range(nrunStart, nrunEnd):
                # uVec = binomial(1, 0.5, Nnomissing) * 2 - 1
                uVec0 = np.random.normal(0, 1, Nnomissing)
                uVec = 2 * uVec0 - 1
                # 通过PCG算法计算目的矩阵 sigma^{-1}u
                Sigma_iu = self.getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG)
                # P*u
                Pu = Sigma_iu - Sigma_iX @ (cov1 @ (Sigma_iXt @ uVec)) 
                # phi*u
                Au = self.getCrossprodMatAndKin(uVec)
                # u^T*phi*P*u
                tempVec[i] = dot(Au, Pu)

            traceCV = self.calCV(tempVec)

            if traceCV > traceCVcutoff:
                nrunStart = nrunEnd
                nrunEnd = nrunEnd + 10
                #tempVec.resize(nrunEnd)
                tempVec = np.resize(tempVec, nrunEnd)
                print(f"CV for trace random estimator using {nrun} runs is {traceCV} > {traceCVcutoff}")
                print(f"try {nrunEnd} runs")

        tra = np.mean(tempVec)
        return tra

    def GetTrace_q(self, Sigma_iX, Xmat, wVec, tauVec, cov1, nrun, maxiterPCG, tolPCG, traceCVcutoff):
        Sigma_iXt = Sigma_iX.T
        Nnomissing = self.dimNum
        tempVec = np.zeros(nrun)
        tempVec0 = np.zeros(nrun)

        nrunStart = 0
        nrunEnd = nrun
        traceCV = traceCVcutoff + 0.1
        traceCV0 = traceCVcutoff + 0.1

        while traceCV > traceCVcutoff or traceCV0 > traceCVcutoff:

            for i in range(nrunStart, nrunEnd):
                uVec = binomial(1, 0.5, Nnomissing) * 2 - 1
                Sigma_iu = self.getPCG1ofSigmaAndVector(wVec, tauVec, uVec, maxiterPCG, tolPCG)
                Pu = Sigma_iu - Sigma_iX @ (cov1 @ (Sigma_iXt @ uVec))
                Au = self.getCrossprodMatAndKin(uVec)
                tempVec[i] = np.dot(Au, Pu)
                tempVec0[i] = np.dot(uVec, Pu)

            traceCV = self.calCV(tempVec)
            traceCV0 = self.calCV(tempVec0)

            if traceCV > traceCVcutoff or traceCV0 > traceCVcutoff:
                nrunStart = nrunEnd
                nrunEnd += 10
                tempVec = np.pad(tempVec, (0, 10), 'constant', constant_values=0)
                tempVec0 = np.pad(tempVec0, (0, 10), 'constant', constant_values=0)
                print(f"CV for trace random estimator using {nrunStart} runs is {traceCV}(> {traceCVcutoff})")
                print(f"try {nrunEnd} runs")

        traVec = np.array([np.mean(tempVec0), np.mean(tempVec)])
        return traVec
    def getAIScore(self, family_indicator,Yvec, Xmat, wVec, tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff):
        '''
        Get AI (Average Information) score, Latex: \frac{- F - H}{2}, where F is the fisher-scoring matrix and H is the Hessian Matrix.

        Args:
            Yvec (array): The response variables, [n*1, Binary];
            Xmat (array): The input variables, Covariates and Genotypes, [n*p, Continuous];
            wVec (array): The variance components mu(1 - mu), (n, );
            tauVec (list): List for initial values for the variance component parameter estimates;
            Sigma_iY (array): The intermediate terms, (n, );
            Sigma_iX (array): The intermediate terms, (n, );
            cov (array): The covariance matrix of inputs vairables X, (p, p);
            nrun (int): ;
            maxiterPCG (int): The max interation time for PCG;
            tolPCG (float): The tolerance of PCG algorithm for whether stop iterations;
            traceCVcutoff (float): Threshold for the coefficient of variation for trace estimation.

        Returns:
            YPAPY (array): The matrix YPAPY, (n, 1);
            Trace (float): The trace of P*phi;
            PY1 (array): The matrix P in paper multiply Y, ();
            AI (array): The matrix YPAPAPY AI_{\tau\tau}, ()
        '''
        if family_indicator == "Binomial":
            Sigma_iXt = Sigma_iX.T
            PY1 = Sigma_iY - Sigma_iX @ (cov @ (Sigma_iXt @ Yvec))
            APY = self.getCrossprodMatAndKin(PY1)
            YPAPY = np.dot(PY1, APY)

            # 计算 tr(P*phi)
            Trace = self.GetTrace(Sigma_iX, Xmat, wVec, tauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff)

            PAPY_1 = self.getPCG1ofSigmaAndVector(wVec, tauVec, APY, maxiterPCG, tolPCG)
            PAPY = PAPY_1 - Sigma_iX @ (cov @ (Sigma_iXt @ PAPY_1))
            AI = np.dot(APY, PAPY)
            return {"YPAPY": YPAPY, "Trace": Trace, "PY": PY1, "AI": AI}
        elif family_indicator == "Beta":

            Sigma_iXt = Sigma_iX.T
            Xmatt = Xmat.T
            sym_matrix = (Xmatt @ Sigma_iX + (Xmatt @ Sigma_iX).T) / 2
            try:
                # cov1 = np.linalg.inv((Xmatt @ Sigma_iX + Sigma_iX @ Xmatt) / 2)  # Ensuring symmetry
                cov1 = np.linalg.inv(sym_matrix)
            except np.linalg.LinAlgError:
                # cov1 = np.linalg.pinv((Xmatt @ Sigma_iX + Sigma_iX @ Xmatt) / 2)
                cov1 = np.linalg.pinv(sym_matrix)
                print("inv_sympd failed, inverted with pinv")
            PY1 = Sigma_iY - Sigma_iX @ (cov1 @ (Sigma_iXt @ Yvec))
            APY = self.getCrossprodMatAndKin(PY1)  # Placeholder for actual implementation
            YPAPY = np.dot(PY1, APY)
            A0PY = PY1  # Quantitative
            YPA0PY = np.dot(PY1, A0PY)  # Quantitative
            Trace = self.GetTrace_q(Sigma_iX, Xmat, wVec, tauVec, cov1, nrun, maxiterPCG, tolPCG,
                                    traceCVcutoff)  # Placeholder for actual implementation

            AI = np.zeros((2, 2))
            PA0PY_1 = self.getPCG1ofSigmaAndVector(wVec, tauVec, A0PY, maxiterPCG,
                                                   tolPCG)  # Placeholder for actual implementation
            PA0PY = PA0PY_1 - Sigma_iX @ (cov1 @ (Sigma_iXt @ PA0PY_1))

            AI[0, 0] = np.dot(A0PY, PA0PY)

            PAPY_1 = self.getPCG1ofSigmaAndVector(wVec, tauVec, APY, maxiterPCG,
                                                  tolPCG)  # Placeholder for actual implementation
            PAPY = PAPY_1 - Sigma_iX @ (cov1 @ (Sigma_iXt @ PAPY_1))
            AI[1, 1] = np.dot(APY, PAPY)
            AI[0, 1] = np.dot(A0PY, PAPY)
            AI[1, 0] = AI[0, 1]
            print("AI", AI)
            return {"YPAPY": YPAPY, "YPA0PY": YPA0PY, "Trace": Trace, "PY": PY1, "AI": AI}


def Covariate_Transform_Back(coef, Param_transform):
    # coef: 模型系数
    # Param_transform: 变换参数
    
    coef1 = coef * np.sqrt(Param_transform['N'])
    
    # 使用 np.linalg.solve 求解线性方程组
    coef_org = np.linalg.solve(Param_transform['qrr'], coef1)
    
    # 将系数命名为 Param_transform['X_name']
    coef_org_dict = dict(zip(Param_transform['X_name'], coef_org))
    
    return coef_org_dict


#
# import statsmodels.api as sm
# import numpy as np
# import patsy
#
# def covariate_transform(formula, data):
#     # Prepare the design matrix X and response variable Y using patsy
#
#     Y, X1 = patsy.dmatrices(formula, data, return_type='dataframe')
#
#     # Get column names from the design matrix
#     X_names = X1.columns
#
#     # Fit the linear regression model without an intercept
#     model = sm.OLS(Y, X1)
#     results = model.fit()
#
#     # Check for NA coefficients which may indicate multicollinearity
#     na_coef_idx = np.where(np.isnan(results.params))[0]
#     if len(na_coef_idx) > 0:
#         print(X1.head())
#         # Drop columns that are causing multicollinearity
#         X1.drop(X1.columns[na_coef_idx], axis=1, inplace=True)
#         print(X1.head())
#         excluded_vars = X_names[na_coef_idx]
#         print("Warning: multi collinearity is detected in covariates!", excluded_vars, "will be excluded in the model")
#         X_names = X_names.delete(na_coef_idx)
#
#     # Adjust the first column name if its index is not in the NA indices
#     if 1 not in na_coef_idx:
#         X_names = X_names.insert(0, 'minus1')
#     data_new =  pd.merge(Y, X1)
#     return data_new
#



def fit_nullGLMM_model(X, y, family):
    start_time = time.time()
    # Fit null model
    if family == 'Binomial':
        y[np.where(y >= 0.5)] = 1
        y[np.where(y < 0.5)] = 0
        y = y.flatten()
        print(type(y))
        print(type(X))
        data_new = pd.DataFrame({"y_binary":y, "x1":X[:, 0], "x2":X[:, 1]})
        formula = 'y_binary ~ x1 + x2'

        # data_new = pd.DataFrame({"y_binary": y, "x1": X[:,0]})
        # formula = 'y_binary ~ x1 '

        model = smf.glm(formula=formula, data=data_new, family=sm.families.Binomial())
        fit0 = model.fit()
        print("fit0.summary()")
        print(fit0.summary())
        eta = fit0.predict(linear=True)
        mu = fit0.predict()
        mu_eta = mu * (1 - mu)

        result = {'model':fit0, 'X':X, 'y':y, 'n':y.shape[0], 'coef':fit0.params[0], 'offset': None, 'family': family, 'eta':eta, 'mu':mu, 'mu_eta':mu_eta}
        print("result",result)
        # fit0 = LogisticRegression()
        # fit0.fit(X, y)
        #
        # P = fit0.predict_proba(X)
        # eta = fit0.decision_function(X)
        # mu = fit0.predict_proba(X)[:, 1]
        # mu_eta = np.exp(mu) / (1 + np.exp(mu))**2
    elif family == 'Poisson':
        fit0 = LinearRegression()
        fit0.fit(X, y)

        # P = fit0.predict_proba(X)
        mu = fit0.predict(X)
        eta = np.log(mu)
        mu_eta = 1 / mu


    elif family == 'Beta':

        epsilon = 1e-8 #simulation data里，1e-9及以下会因为0太多报错
        y = np.where(y == 0, epsilon, y)
        y = np.where(y == 1, 1 - epsilon, y)

        betareg = importr('betareg')
        df = pd.DataFrame(X, columns=[f'X{i + 1}' for i in range(X.shape[1])])
        df['y'] = y
        # print("NAN:",df.isnull().values.any())
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_dataframe = ro.conversion.py2rpy(df)

        formula = Formula('y ~ ' + ' + '.join(df.columns[:-1]))

        try:
            fit0 = betareg.betareg(formula=formula, data=r_dataframe, link="logit")
        except Exception as e:
            print(e)

        mu = np.array(r['predict'](fit0, type="response"))
        # eta = np.log(mu / (1 - mu))
        # mu_eta = mu * (1 - mu)
        eta = np.array(r['predict'](fit0, type="link"))
        mu_eta = np.exp(eta) / (1 + np.exp(eta))**2
        mean_coefficients = np.array(fit0.rx2('coefficients')[0])
        intercept_coef = mean_coefficients[0]
        precision_coefficients = np.array(fit0.rx2('coefficients')[1])
        result = {'model': fit0, 'X': X, 'y': y, 'n': y.shape[0],'coef': intercept_coef,'offset': None, 'family': family, 'eta': eta, 'mu': mu,
                   'mu_eta': mu_eta}
        end_time = time.time()
        # print(r.summary(fit0))
        print('Time taken to Beta Regression:', end_time-start_time)

    return result