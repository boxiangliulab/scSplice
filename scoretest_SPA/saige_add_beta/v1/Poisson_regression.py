from scipy.stats import poisson
import pandas as pd
import sys
import statsmodels.api as sm
import statsmodels.formula.api as smf

data = pd.read_csv(sys.path[0] + "/sample_data/poisson_sim.csv") # y:获奖数量，x:学科内容, 数学期末考试分数 <-> id,num_awards,prog,math

poisson_model = smf.glm(formula="num_awards ~ prog + math", data=data, family=sm.families.Poisson())
result_poisson = poisson_model.fit()
print(result_poisson.summary())

