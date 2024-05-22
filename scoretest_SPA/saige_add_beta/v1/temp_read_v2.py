import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, cauchy

import sys

# 生成一些数据用于绘图
x = np.linspace(-10, 10, 1000)
y_cauchy = cauchy.pdf(x, 0, 1)  # 柯西分布的概率密度函数
y_norm = norm.pdf(x, 0, 1)       # 标准正态分布的概率密度函数

# 绘制图像
plt.figure(figsize=(10, 6))
plt.plot(x, y_cauchy, label='Cauchy Distribution', color='blue')
plt.plot(x, y_norm, label='Standard Normal Distribution', color='red', linestyle='--')
plt.title('Comparison of Cauchy Distribution and Standard Normal Distribution')
plt.xlabel('x')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True)
plt.savefig(sys.path[0] + '/cauchy_distribution.png', bbox_inches='tight')

