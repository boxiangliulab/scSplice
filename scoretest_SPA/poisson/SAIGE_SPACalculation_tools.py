import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.optimize import root_scalar
from scipy import stats
from scipy.stats import chi2
from SAIGE_FitGLMM_tools import fitGLMM_tools

from scipy.linalg import qr, svd

# Setting the seed
np.random.seed(1234)


def ScoreTest_wSaddleApprox_Get_X1(X1):
    q1 = X1.shape[1]
    if q1 >= 2:
        # 检查第一列和第二列是否完全相同
        if np.sum(np.abs(X1[:, 0] - X1[:, 1])) == 0:
            X1 = np.delete(X1, 1, axis=1)
            q1 -= 1

    # 进行 QR 分解
    qr1 = qr(X1, mode='economic', pivoting=True)
    rank = np.sum(np.abs(qr1[1].diagonal()) > np.finfo(float).eps * 100)

    if rank < q1:
        # 奇异值分解
        u, _, _ = svd(X1, full_matrices=False)
        X1 = u[:, :rank]

    return X1

def Estimation_r(Psi, W, tau):
    '''
    Calculate the constant ratio r.

    Args:
        Psi (sp.coomatrix): Sparse Geneomic Relationship Matrix;
        W (array): mu(1-mu) is the variance matrix;

    Returns:
        r = frac{sum_{i=1}^{n}frac{lambda_i}{ bar{w} +  hat{ tau} lambda_i}}{sum_{i=1}^{n} frac{ lambda_i}{ bar{w}}}.
    '''
    Lambda = Psi.data
    w = W
    bar_w = np.mean(1/w)
    hat_tau = tau
    r = np.sum(Lambda/(bar_w + hat_tau * Lambda)) / np.sum(Lambda / bar_w)
    c = 1/np.sqrt(np.sum(Lambda/bar_w))
    return r,c


def Estimation_T_adj(G, X, W, Y, mu, r):
    '''
    Calculate the test statistics T_adj.

    Args:
        G (array): Genotype, allele count (0, 1, 2);
        X (array): covariate;
        W (array): mu(1-mu) is the variance matrix;
        Y (array): response variable, case-control;
        mu (array): fitted probability.

    Returns:
        T_adj: Adjusted score test statistic.
    '''

    dimNum = X.shape[0]
    ValueVec = np.sqrt(W)
    # print("valueVec", ValueVec)
    LocationMat = np.concatenate((np.array([i for i in range(dimNum)]).reshape(-1, 1), np.array([i for i in range(dimNum)]).reshape(-1, 1)), axis = 1)
    # print("LocationMat",LocationMat)
    fit_tools = fitGLMM_tools(ValueVec, LocationMat, dimNum)
    # print("X",X)
    sqrt_WX = fit_tools.getCrossprodMatAndKin_Tadj(X)
    # print("sqrt_WX",sqrt_WX)
    XWX = np.dot(sqrt_WX.T, sqrt_WX)
    XWX_iXT = fit_tools.spsolove_for_Tadj(XWX, X.T)
    # print('XWX', XWX)
    # print('XWX_iXT', XWX_iXT)
    g = G - X @ XWX_iXT @ np.diag(W) @ G

    tem = r * (g.T @ np.diag(W) @ g)
    denominator = np.sqrt(tem)
    nominator = g.T @ (Y - mu.reshape(-1, 1))
    T = nominator/denominator

    return T[0, 0], g







class SPA():
    def __init__(self, G, mu,r,c,W,y,q):
        self.G = G
        self.mu = mu
        self.n = mu.shape[0] 
        self.r = r
        self.c = c
        self.W = W
        self.y = y
        # self.T_adj = T_adj
        self.T_adj = q
        self.q = q

    def Korg(self,t):
          # Initialize the output array with zeros
        n=self.n
        temp = np.log(np.ones(n) - self.mu + self.mu * np.exp(self.G * t))
          # sum(temp) seems redundant unless temp is vectorized
        return np.sum(temp)

    def K1_adj(self, t, q):
        temp1 = (np.ones(self.n) - self.mu) * np.exp(-self.G * t) + self.mu
        temp2 = self.mu * self.G
        result = temp2/temp1
        result[np.isinf(result)]=np.nan
        out = np.nansum(result) - q
        return out 

    def K2(self,t):
          temp1 = ((np.ones(self.n) - self.mu) * np.exp(-self.G * t) + self.mu)**2
          temp2 = (np.ones(self.n) - self.mu) * self.mu * self.G**2 * np.exp(-self.G * t)
          result = temp2/temp1
          result[np.isinf(result)]=np.nan
          out = np.nansum(temp2 / temp1)
          return out  

    
    def getroot_K1(self,init,q):
        newG = np.array([self.G[i][0] for i in range(len(self.G))])  
        g_pos = np.sum(newG[newG>0])
        g_neg = np.sum(newG[newG<0])
        tol = 1e-4
        conv = False
        if (q >= g_pos) | (q <= g_neg):
            conv = False
            return np.inf,0,conv
        else:
            t=init
            maxiter=1000
            K1_eval=self.K1_adj(t,q)
            prevJump = np.inf
            rep=1
            while rep<=maxiter:
                K2_eval=self.K2(t)
                tnew=t-K1_eval/K2_eval
                if tnew is None:
                    conv=False
                    break
                if abs(tnew-t)<tol:
                    conv=True
                    break
                newK1 = self.K1_adj(tnew,q)
                if np.sign(K1_eval)!=np.sign(newK1):
                    if abs(tnew-t)>prevJump-tol:
                        tnew=t+np.sign(newK1-K1_eval)*prevJump/2
                        newK1=self.K1_adj(tnew,q)
                        prevJump=prevJump/2
                    else:
                        prevJump=abs(tnew-t)
                rep=rep+1
                t=tnew
                K1_eval=newK1
            return t,rep,conv

            
    def solve_t(self,x, q):
        k1 = self.Korg(x)
        k2 = self.K2(x)
        if (k1 is not None) & (k2 is not None):
            tmp1 = x * q - k1
            w = np.sign(x) * (2*tmp1)**(1/2)
            v = x * (k2)**(1/2)
            Ztest = w + 1/w * np.log(v/w)
            if Ztest > 0:
                pval = stats.norm.sf(Ztest)
            else:
                pval = -stats.norm.cdf(Ztest)
        else:
            pval = np.inf
        return pval

    def calculate_F(self,pval_noadjust,qinv,var1,m1,cutoff=2):
        pval=pval_noadjust
        if abs(self.T_adj-m1)/np.sqrt(var1) > cutoff:

            zeta1,rep1,conv1 = self.getroot_K1(0,self.T_adj)
            zeta2,rep2,conv2 = self.getroot_K1(0,qinv)
            try:
                 p1 = self.solve_t(zeta1,self.T_adj)
            except Exception as e:
                 p1 = pval_noadjust/2
                 return p1
            try:
                 p2 = self.solve_t(zeta2,qinv)
            except Exception as e:
                 p2 = pval_noadjust/2
                 return p2
            pval = abs(p1)+abs(p2)
        else:
            pval = pval_noadjust
        if (pval!=0) & (pval_noadjust/pval > 1e3):
             return self.calculate_F(pval_noadjust,qinv,var1,m1,cutoff=cutoff*2)
        else:
             return pval_noadjust
        
    def calculate_pval_noadjust(self):
        q = self.q
        newG = [self.G[i][0] for i in range(len(self.G))]    
        newG2=[i*i for i in newG]
        var1 = np.sum(self.mu * (np.ones(self.n)-self.mu) * newG2)
        m1 = np.sum(self.mu * newG)
        S=q-m1
        val=S*S/var1

        pval_noadjust = chi2.sf(val,df=1)

        qinv= -np.sign(q-m1) * abs(S) + m1

        return pval_noadjust,qinv,var1,m1

    

    

# if __name__ == "__main__":

#     spa_stdnormal = SPA()


#     # Bootstrapping
#     # boots = np.array([np.mean(np.random.choice(x, size=len(x), replace=True)) for _ in range(100)])

#     # Generating random exponential data
#     boots = np.random.normal(0, 1, 1000)
#     # boot = [np.random.choice(sample, ) for i in range(100)]

#     # Plotting the histogram
#     plt.hist(boots, bins=30, density=True, alpha=0.6, color='g')

#     # Adding the saddlepoint approximation plot
#     fhat_values = np.linspace(-3, 3, 100)  # Adjust range and number of points as needed
#     density_values = [spa_stdnormal.f_hat(i) for i in fhat_values]
#     plt.plot(fhat_values, np.array(density_values), color='red')

#     plt.xlabel('Value')
#     plt.ylabel('Probability Density')
#     plt.title('Bootstrapped Sample Means with Saddlepoint Approximation')
#     plt.savefig('spa.png', bbox_inches='tight')

