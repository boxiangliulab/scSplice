import pandas as pd
import numpy as np

def read_cov_and_pheno(clu_index,path_pca,path_expression):
    df_pca = pd.read_csv(path_pca, sep='\s+', index_col=0)
    df_pca.reset_index(inplace=True)
    df_pca.rename(columns={'index': 'sample1'}, inplace=True)
    df_expr = pd.read_csv(path_expression, sep='\s+', index_col=0)

    samples_list = []
    expr_index = df_expr.iloc[clu_index].index

    for x in range(len(expr_index)):
        sample1 = expr_index[x].split("SEA")[0]
        sample2 = "SEA" + expr_index[x].split("SEA")[1].rstrip('.bam')
        samples_list.append((sample1, sample2))
    samples_df = pd.DataFrame(samples_list, columns=['sample1', 'sample2'])
    expr = np.array(df_expr.iloc[clu_index])
    for x in range(len(expr)):
        splice = float(expr[x].split(":")[0])
        total = float(expr[x].split(":")[1])
        samples_df.loc[x, 'splice'] = splice
        samples_df.loc[x, 'total'] = total
    data = pd.merge(samples_df, df_pca, on='sample1', how='left')
    data['sample'] = data['sample1']+data['sample2']
    data.drop(['sample1', 'sample2'], axis=1, inplace=True)
    y = data['splice'].to_numpy(dtype=int)
    X = data.drop(['splice', 'sample'], axis=1).to_numpy()
    return X,y

def read_geno(geno_PATH):
    pass