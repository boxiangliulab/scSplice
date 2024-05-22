import numpy as np
import pandas as pd
from scipy import sparse
import sys
import gzip
from scipy.sparse import coo_matrix
from scipy.io import mmread
import os

class Read_in_tools:
    def __init__(self, prefix, cell, geno):
        self.prefix = prefix
        self.cell = cell
        self.geno = geno

    def convert_data(self, x):
        '''
        Function to help convert the GT data to categorical variables.
        '''
        if len(x) > 3:
            x = x[:3]
        else:
            pass

        if x == '0|0':
            t = 0
        elif x == '0|1' or x == '1|0':
            t = 1
        elif x == '1|1':
            t = 2
        else:
            t = 3

        # if x == '0/0':
        #     t = 0
        # elif x == '0/1':
        #     t = 1
        # elif x == '1/0':
        #     t = 2
        # else:
        #     t = 3

        return t



    def Covariate_and_Phenotype(self, index,family_indicator):
        '''
        Help to read in Phenotype and the covariates.
        TODO: Need modification based on the upstream inputs
        '''
        # if family_indicator == "Binomial":
        #     data = pd.read_csv('pheno_1000samples.txt_withdosages_withBothTraitTypes.csv', sep=" ")
        #     self.data_df = data[['IID', 'x1', 'x2', 'y_binary']]
        #     self.data_df = self.data_df.rename(columns={'IID': 'id'})
        #     print(self.data_df)
        #
        #     X = np.array(self.data_df[['x1', 'x2']])
        #     y = np.array(self.data_df[['y_binary']])
        #
        #     return X , y
        # elif family_indicator == "Beta":
        #     data = pd.read_csv('pheno_1000samples.txt_withdosages_withBothTraitTypes.csv', sep=" ")
        #     self.data_df = data[['IID', 'x1', 'x2', 'y_quantitative']]
        #     self.data_df = self.data_df.rename(columns={'IID': 'id'})
        #     print(self.data_df)
        #
        #     X = np.array(self.data_df[['x1', 'x2']])
        #     y = np.array(self.data_df[['y_quantitative']])
        #     y = 1/(1+np.exp(-y))
        #
        #     return X, y

        # Cell_path = self.prefix + 'Covariate/' + self.cell + '_PC.txt'
        Cell_path = self.prefix + 'simu_covariate.txt'
        X_new = pd.read_csv(Cell_path, sep = '\t')

        # 选择 'id' 列和 'S101' 列到最后的所有列
        "x y如果id对应的话就不需要"
        X_new = X_new.loc[:, ['id'] + [f'S{i}' for i in range(101, 501)]]


        X = X_new.iloc[-2:]
        X = np.array(X.T)[1:, :]
        id = np.array(X_new.columns.values.tolist()[1:])
        data = {'id': id, 'sex': X[:, 0], 'age': X[:, 1]}

        # Phenotype_path = self.prefix + 'Phenotype/' + self.cell + self.geno + '_beta.tmp.txt.gz'
        Phenotype_path = self.prefix + "simu_scSplice_pheno.txt.gz"
        y_new = pd.read_csv(Phenotype_path, sep = '\t' , compression='gzip')
        yy = y_new.T[index]
        yy_list = list(map(float, yy.values[0].split(' ')[1:]))
        y = np.array(yy_list)

        #正则化
        y_min = y.min()
        y_max = y.max()
        y = (y - y_min) / (y_max - y_min)

        data['phenotype'] = y

        self.data_df = pd.DataFrame(data)
        # print(self.data_df)
        return X, y

    def Get_Covariate_and_Phenotype(self):
        X = np.array(self.data_df['sex','age'])
        y = np.array(self.data_df['phenotype'])
        return X, y
    
    def Specific_Genotype(self, index):
        '''
        Help to read in the VCF files.
        '''
        # sample_data = []
        Geno_path = self.prefix + "simu_geno.vcf.gz"
        # Geno_path = self.prefix + 'Genotype/filter.chr' + self.geno + '.dose.vcf.gz'
        # Geno_path = r"./Experiment_data/Genotype/genotype_100markers.vcf.gz"
        sample_dic = {}
        with gzip.open(Geno_path, 'rb') as vcf_in:
            line_count = 0

            for line in vcf_in:
                line = line.strip()
                try:
                    line = str(line, encoding = "utf-8")
                except (UnicodeDecodeError, AttributeError):
                    pass

                if line[:2] == "##":
                    pass
                else:
                    if line_count == 0:
                        sample_header = line.split('\t')[9:]
                        #sample_header = [x.split('_')[0] for x in sample_header]
                        sample_dic['id'] = sample_header
                        line_count += 1
                    elif line_count == index:
                        line_content = line.split('\t')[9:]
                        line_tem = list(map(self.convert_data, line_content))
                        sample_dic['geno' + str(index)] = line_tem
                        line_count += 1
                        break
                    else:
                        line_count += 1
                        pass

        self.data_df_geno = pd.DataFrame(sample_dic)
        self.data_full = pd.merge(self.data_df, self.data_df_geno, on='id', how = 'left')

        # print("self.data_full",self.data_full)

        # file_name = f".\geno\genotype_{index}.txt"
        # self.data_full.to_csv(file_name, sep='\t', index=False)
        return np.array(self.data_full[['geno' + str(index)]])


    def read_in_Genotype(self, M):
        '''
        Help to read in the VCF files.
        '''

        # Geno_path = self.prefix + 'Genotype/filter.chr' + self.geno + '.dose.vcf.gz'
        Geno_path = self.prefix + "simu_geno.vcf.gz"
        # Geno_path = r"./Experiment_data/Genotype/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.vcf.gz"

        sample_dic = {}
        sample_data = []
        # with open(vcf_file_path, 'r') as vcf_in:
        with gzip.open(Geno_path, 'rb') as vcf_in:

            line_count = 0
            for line in vcf_in.readlines():
                try:
                    line = str(line, encoding = "utf-8")
                except (UnicodeDecodeError, AttributeError):
                    pass

                if line[:2] == "##":
                    pass
                else:
                    if line_count == 0:
                        sample_header = line.split('\t')[9:]
                        sample_header = [x.split('_')[0] for x in sample_header]
                        sample_dic['id'] = sample_header
                        line_count += 1
                    else:
                        line_content = line.split('\t')[9:]
                        # line_content.append(line_content[-1])
                        line_tem = list(map(self.convert_data, line_content))
                        sample_data.append(line_tem)
                        sample_dic['geno' + str(line_count)] = line_tem
                        line_count += 1
                        # For sample data to reduce time
                        if line_count == M + 1:
                            break
        
        self.data_df_geno = pd.DataFrame(sample_dic)
        self.data_full = pd.merge(self.data_df, self.data_df_geno, on='id', how = 'left')
        
        return np.array(self.data_full[['geno' + str(i+1) for i in range(M)]])

    def getsubGRM(self,sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, modelID=None):
        # 检查文件存在性
        if not os.path.exists(sparseGRMFile):
            raise FileNotFoundError(f"ERROR! Sparse GRM file {sparseGRMFile} does not exist")
        if not os.path.exists(sparseGRMSampleIDFile):
            raise FileNotFoundError(f"ERROR! Sparse GRM Sample ID file {sparseGRMSampleIDFile} does not exist")

        # 读取稀疏GRM文件
        print("Extracting sparse GRM...")
        sparseGRMLarge = mmread(sparseGRMFile).tocsr()
        print("Initial non-zero elements:", sparseGRMLarge.nnz)

        # 将小于等于相关性阈值的元素设置为零并移除
        sparseGRMLarge.data[sparseGRMLarge.data <= relatednessCutoff] = 0
        sparseGRMLarge.eliminate_zeros()
        print("Non-zero elements after applying threshold:", sparseGRMLarge.nnz)

        # 读取样本ID文件
        sparseGRMSampleID = pd.read_csv(sparseGRMSampleIDFile, header=None, names=["sampleID"], dtype={'sampleID': str})
        sparseGRMSampleID['IndexGRM'] = np.arange(sparseGRMLarge.shape[0])

        # 如果提供了模型ID
        if modelID is not None:
            sampleInModel = pd.DataFrame(modelID, columns=['IID'])
            sampleInModel['IndexInModel'] = np.arange(len(sampleInModel))
            print(len(sampleInModel), "samples have been used to fit the GLMM null model.")

            # 合并模型ID和样本ID
            mergeID = pd.merge(sampleInModel, sparseGRMSampleID, left_on="IID", right_on="sampleID", how='inner')
            if len(sampleInModel) > len(mergeID):
                missing = len(sampleInModel) - len(mergeID)
                raise Exception(f"ERROR: {missing} samples used for model fitting are not in the specified GRM")

            # 提取子矩阵
            indexIDofGRM = mergeID['IndexGRM'].values
            sparseGRM = sparseGRMLarge[indexIDofGRM, :][:, indexIDofGRM]
            return sparseGRM

        return None



    def Create_GRM(self, M):
        '''
        Help to read in the VCF files and convert to genotype.
        '''
        sample_data = self.read_in_Genotype(M)
        
        GRM = np.cov(sample_data)
        GRM = np.where(abs(GRM) > 0.001, GRM, 0)
        # GRM = np.where(GRM < -0.001, GRM, 0)
        setDiagasOne = np.eye(GRM.shape[0])
        GRM += setDiagasOne
        GRM = np.where(GRM<1, GRM, 1)
        sGRM = sparse.coo_matrix(GRM)
        
        return sGRM


if __name__ == "__main__":
    import sys

    prefix = sys.path[0]
    Geno_name = '21'
    Cell_name = 'CD4+_T_naive'
    i = 1

    read_in_tools = Read_in_tools(prefix + "/Experiment_data/", Cell_name, Geno_name)
    X, y = read_in_tools.Covariate_and_Phenotype(index = i)

    # Geno =  read_in_tools.Create_GRM(M = 20)
    geno_2 = read_in_tools.Specific_Genotype(index = 2)
    print('end')