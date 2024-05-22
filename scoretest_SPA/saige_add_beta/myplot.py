# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import numpy as np
# from scipy.stats import pearsonr
#
#
#
# df1 = pd.read_csv(r"C:\study\SAIGE_debuggingGLMM\pvalues.csv",usecols=['geno_index', 'p_value'])
# df2 = pd.read_csv(r"C:\study\SAIGE_debuggingGLMM\Experiment_data\genotype_100markers_marker_plink_step1withSparseGRM_Firth.txt", sep='\t', usecols=['POS', 'p.value'])
# df3 = pd.read_csv(r"C:\study\SAIGE_debuggingGLMM\geno\pearson_coefficients.csv",usecols=['Index', 'PValue'])
# #df1_filtered = df1[df1['p_value'] < 0.05]
# #df2_filtered = df2[df2['p.value'] < 0.05]
# print(df1.head(10))
# print(df2.head(10))
# print(df3.head(10))
#
# merged_df12 = pd.merge(df1, df2, right_on='POS', left_on='geno_index')
# print(merged_df12.head(10))
# merged_df12['neg_log10_pvalue'] = -np.log10(merged_df12['p_value'])
# merged_df12['neg_log10_p_dot_value'] = -np.log10(merged_df12['p.value'])
# #
# #
# # merged_df13 = pd.merge(df1, df3, right_on='Index', left_on='geno_index')
# # print(merged_df13.head(10))
# # merged_df13['neg_log10_pvalue'] = -np.log10(merged_df13['p_value'])
# # merged_df13['neg_log10_p_dot_value'] = -np.log10(merged_df13['PValue'])
#
# #
# # merged_df23 = pd.merge(df2, df3, right_on='Index', left_on='POS')
# # print(merged_df23.head(10))
# # merged_df23['neg_log10_pvalue'] = -np.log10(merged_df23['PValue'])
# # merged_df23['neg_log10_p_dot_value'] = -np.log10(merged_df23['p.value'])
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# plt.figure(figsize=(10, 6))
# plt.scatter(merged_df12['neg_log10_pvalue'], merged_df12['neg_log10_p_dot_value'], alpha=0.5)
# plt.title('Scatter Plot of -log10(pearson) vs -log10(p_adjusted)')
# plt.xlabel('-log10(pearson)')
# plt.ylabel('-log10(p_adjusted)')
# plt.grid(True)
# plt.show()
#
#
# coefficient, p_value = pearsonr(merged_df12['neg_log10_pvalue'], merged_df12['neg_log10_p_dot_value'])
# print(coefficient, p_value )
#
#








import pandas as pd

# 读取CSV文件
df = pd.read_csv(r"pvalues.csv")

# 找到所有小于0.05的值
low_values = df[df < 0.05].stack()

# 创建一个包含行、列和值的DataFrame
result = pd.DataFrame(list(low_values.index), columns=['Row', 'Column'])
result['Value'] = low_values.values

# 保存结果到新的CSV文件
result.to_csv('significant_python_pvalues.csv', index=False)



