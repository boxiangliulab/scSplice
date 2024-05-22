import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from SPA_test import SPA

import pandas as pd
import sys

data = pd.read_csv(sys.path[0] + "/sample_data/poisson_sim.csv") # y:获奖数量，x:学科内容, 数学期末考试分数 <-> id,num_awards,prog,math
print(max(data['num_awards']))
data['num_awards'][data['num_awards'] >= 1] = 1
data['num_awards'][data['num_awards'] < 1] = 0
data = data.drop(columns='id')

X = np.array(data['math']).reshape(-1, 1)
G = np.array(data['prog']).reshape(-1, 1)
y = np.array(data['num_awards']).reshape(-1, 1)
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

from sklearn.linear_model import LogisticRegression
logit_model = LogisticRegression()
logit_model.fit(X, y)
P = logit_model.predict_proba(X)
mu = P[:, 1].reshape(-1, 1)

X = np.array(X)
G = np.array(G)
W = np.diag(np.diagonal(P[:, 1].reshape(-1, 1) @ P[:, 0].reshape(-1, 1).T))
# G_tilde = G - X(X^TWX)^{-1} X^TWG
G_tilde = G - X @ np.linalg.inv(X.T @ W @ X) @ X.T @ W @ G # 预计计算会有问题矩阵很大

S = np.sum(G_tilde * (y - mu))

spa = SPA(G_tilde, mu)
p_value = spa.calculate_F(S)
print(p_value)