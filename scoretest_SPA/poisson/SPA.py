import numpy as np
from scipy.stats import chi2, norm,binom, poisson




def add_logp(p1, p2):
    p1 = -np.abs(p1)
    p2 = -np.abs(p2)
    maxp = max(p1, p2)
    minp = min(p1, p2)
    return maxp + np.log(1 + np.exp(minp - maxp))



def Korg(t, mu, g):
    n_t = len(t)
    out = np.zeros(n_t)
    for i in range(n_t):
        t1 = t[i]
        temp = mu * np.exp(g * t1)  # 调整泊松均值
        out[i] = np.sum(np.log(temp))  # 计算对数概率
    return out

def K1_adj(t, mu, g, q):
    n_t = len(t)

    out = np.zeros(n_t)
    for i in range(n_t):
        t1 = t[i]
        """二项分布是temp1<-(1 - mu)* exp(-g * t1) + mu"""
        temp1 = mu * np.exp(-g * t1)  # 对应泊松分布
        temp2 = mu * g
        out[i] = np.sum(temp2 / temp1) - q
    return out

def K2(t, mu, g):
    n_t = len(t)

    out = np.zeros(n_t)
    for i in range(n_t):
        t1 = t[i]
        """二项分布是temp1<-((1 - mu)* exp(-g * t1) + mu)^2"""
        temp1 = (mu * np.exp(-g * t1)) ** 2  # 对应泊松分布
        temp2 = mu * g**2 * np.exp(-g * t1)
        out[i] = np.sum(temp2 / temp1)
    return out


def get_saddle_prob(zeta, mu, g, q, log_p=False):
    k1 = Korg(zeta, mu, g)
    k2 = K2(zeta, mu, g)

    if np.isfinite(k1) and np.isfinite(k2):
        temp1 = zeta * q - k1
        w = np.sign(zeta) * np.sqrt(2 * temp1)
        v = zeta * np.sqrt(k2)
        Z_test = w + (1 / w) * np.log(v / w)
        if Z_test > 0:
            pval = norm.sf(Z_test, loc=0, scale=1)  # lower.tail=False
        else:
            pval = -norm.sf(-Z_test, loc=0, scale=1)  # lower.tail=True
        if log_p:
            pval = np.log(pval)
    else:
        pval = -np.inf if log_p else 0

    return pval


def getroot_K1(init, mu, g, q, m1=None, tol=np.finfo(float).eps ** 0.25, maxiter=1000):
    g_pos = np.sum(g[g > 0])
    g_neg = np.sum(g[g < 0])

    if q >= g_pos or q <= g_neg:
        return {'root': np.inf, 'n_iter': 0, 'Is_converge': True}
    else:
        t = init
        K1_eval = K1_adj(t, mu, g, q)
        prevJump = np.inf
        rep = 1
        while True:
            K2_eval = K2(t, mu, g)
            tnew = t - K1_eval / K2_eval
            if np.isnan(tnew):
                conv = False
                break
            if np.abs(tnew - t) < tol:
                conv = True
                break
            if rep == maxiter:
                conv = False
                break
            newK1 = K1_adj(tnew, mu, g, q)
            if np.sign(K1_eval) != np.sign(newK1):
                if np.abs(tnew - t) > prevJump - tol:
                    tnew = t + np.sign(newK1 - K1_eval) * prevJump / 2
                    newK1 = K1_adj(tnew, mu, g, q)
                    prevJump /= 2
                else:
                    prevJump = np.abs(tnew - t)
            rep += 1
            t = tnew
            K1_eval = newK1
        return {'root': t, 'n_iter': rep, 'Is_converge': conv}


def Saddle_Prob(q, mu, g, Cutoff=2, alpha = 5*10e-8, log_p=False):
    m1 = np.sum(mu * g)

    """这里算的时候，泊松的均值=方差。如果是二项分布就需要写成 mu(1-mu)* g ** 2"""
    var1 = np.sum(mu * g ** 2)

    p1 = None
    p2 = None

    Score = q - m1
    qinv = -np.sign(q - m1) * np.abs(q - m1) + m1

    pval_noadj = chi2.sf((q - m1) ** 2 / var1, df=1)

    if log_p:  #是否以对数形式返回P值
        pval_noadj = np.log(pval_noadj)

    Is_converge = True

    if Cutoff == "BE":
        """这里算的时候，泊松的均值=方差。
        如果是二项分布就需要写成 rho = np.sum((np.abs(g)**3) * mu * (1 - mu) * (mu**2 + (1 - mu)**2))"""
        rho = np.sum(  (np.abs(g) ** 3) * mu * (mu ** 2 )  )
        B = 0.56 * rho * var1 ** (-3 / 2)
        p = B + alpha / 2
        # Cutoff = 0.01 if p >= 0.496 else norm.ppf(1 - p)
        Cutoff = norm.ppf(p, loc=0, scale=1, lower=False) if p >= 0.496 else 0.01
    elif Cutoff < 10 ** -1:
        Cutoff = 10 ** -1



    pval = pval_noadj
    """只做了score test，下面spa还没写"""
    # if np.abs(q - m1) / np.sqrt(var1) < Cutoff:
    #     pval = pval_noadj
    # else:
    #     out_uni1 = getroot_K1(0, mu=mu, g=g, q=q)
    #     out_uni2 = getroot_K1(0, mu=mu, g=g, q=qinv)
    #     if out_uni1['Is_converge'] and out_uni2['Is_converge']:
    #         try:
    #             p1 = get_saddle_prob(out_uni1['root'], mu, g, q, log_p=log_p)
    #             p2 = get_saddle_prob(out_uni2['root'], mu, g, qinv, log_p=log_p)
    #             if log_p:
    #                 pval = add_logp(p1, p2)
    #             else:
    #                 pval = np.abs(p1) + np.abs(p2)
    #         except Exception as e:
    #             if log_p:
    #                 pval = pval_noadj - np.log(2)
    #             else:
    #                 pval = pval_noadj / 2
    #         Is_converge = True
    #     else:
    #         print("Error_Converge")
    #         pval = pval_noadj
    #         Is_converge = False



    return {
        "p_value": pval,
        "p_value_NA": pval_noadj,
        "Is_converge": Is_converge,
        "Score": Score
    }


def TestSPA(G, null_model_result, Cutoff=2, alpha = 5*10e-8,   log_p=False):
    y = null_model_result['y']
    mu = null_model_result['mu']
    res = null_model_result['res']

    n_g = np.sum(G)   ##G是geno

    if n_g / (2 * len(G)) > 0.5:
        G = 2 - G
        n_g = np.sum(G)



    G1 = G - null_model_result['XXVX_inv'] @ (null_model_result['XV'] @ G)

    q = np.sum(G1 * y)
    out = Saddle_Prob(q, mu=mu, g=G1, Cutoff=Cutoff, alpha=5*10e-8, log_p=log_p)
    return out













# num_simulations = 1000
# alpha = 0.05
# n = 1000
#
# # 固定随机种子
# np.random.seed(123)
#
# # 模拟 p 值
# p_values = []
# for _ in range(num_simulations):
#     sample1 = binom.rvs(n=1, p=0.02, size=n)
#     sample2 = poisson.rvs(mu=10, size=n)
#     test_result = TestSPA(sample1, sample2)
#     p_value = test_result['p_value']
#     print(p_value)
#     p_values.append(p_value)
#
# # 计算经验 I 型错误率
# type_I_error_rate = np.mean(np.array(p_values) < alpha)

# 打印 I 型错误率
# print("Type I error rate:", type_I_error_rate)