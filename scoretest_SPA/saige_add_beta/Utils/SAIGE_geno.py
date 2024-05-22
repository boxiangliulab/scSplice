import numpy as np
import pandas as pd
import sys
import math
import time

class SparseGRMUsingOneMarker:
    def __init__(self, GRMvec):
        self.GRMvec = GRMvec

    def __call__(self, geno, begin, end):
        for i in range(begin, end):
            iint = geno.indiceVec[i][0]
            jint = geno.indiceVec[i][1]

            ival = geno.m_OneSNP_Geno(iint)
            jval = geno.m_OneSNP_Geno(jint)

            self.GRMvec[i] = geno.sKinLookUpArr[ival][jval]



class geno():
    def __init__(self):
        # 设置一些 public variables
        self.m_bits_val = [0 for i in range(8)]

    def Init_OneSNP_Geno(self):
        self.m_size_of_esi = (self.Nnomissing + 3) / 4
        k = 8
        while k > 0:
            k -= 1
            self.m_bits_val[k] = 1

    def setGenotype(c, pos, geno):
        """
        Set the genotype in the given byte at the specified position.

        Parameters:
            c: bytearray
                The byte to modify.
            pos: int
                The position within the byte (0-3).
            geno: int
                The genotype value (0, 1, 2, or 3).

        Returns:
            None
        """
        c[pos] = geno

        

    def setGenoObj(self, genofile, subSampleInGeno, memoryChunk, isDiagofKinSetAsOne):
        '''
        这个函数的主要目的是从二进制 bed 文件中读取基因型数据，并进行一些预处理，包括处理缺失数据、计算频率等
        '''

        # 初始化一些基本参数
        self.setKinDiagtoOne = isDiagofKinSetAsOne
        self.ptrsubSampleInGeno = subSampleInGeno
        self.Nnomissing = len(subSampleInGeno)

        self.alleleFreqVec = [] 
        self.MACVec = []
        self.invstdvVec = []

        self.N = 0
        self.M = 0 # TODO: 要设置成M个基因marker 个数

        bedfile = genofile + ".bed"
        bimfile = genofile + ".bim"
        famfile = genofile + ".fam"
        
        # 尝试打开famfile
        try:
            with open(famfile, 'r') as test_famfile:
                print("famfile exist!")
        except FileNotFoundError:
            print("Error! famfile not open!")
            return

        indexRow = 0
        with open(famfile, 'r') as test_famfile:
            for _ in test_famfile:
                indexRow += 1

        self.N = indexRow

        # 尝试打开bimfile
        try:
            with open(bimfile, 'r') as test_bimfile:
                print("bimfile exist!")
        except FileNotFoundError:
            print("Error! bimfile not open!")
            return
        
        indexRow = 0
        with open(bimfile, 'r') as test_bimfile:
            for _ in test_bimfile:
                indexRow += 1
        self.M = indexRow

        indexRow = 0
        buffer = 0
        TotalRead = 0

        genoVecOneMarkerOld = []
        genoVecOneMarkerNew = []

        # 计算所需的字节数
        nbyteOld = math.ceil(self.N / 4)
        nbyteNew = math.ceil(self.Nnomissing / 4)
        reserve = math.ceil(self.Nnomissing / 4) * self.M + self.M * 2

        print("nbyte:", nbyteOld)
        print("nbyte:", nbyteNew)
        print("reserve:", reserve)

        # 预留空间并调整大小
        genoVecOneMarkerOld = []

        # 尝试打开bedfile
        try:
            with open(bedfile, 'r') as test_bedfile:
                print("bedfile exist!")
        except FileNotFoundError:
            print("Error! bedfile not open!")
            return
        
        print('M:', self.M, 'N:', self.N)
        numMarkersofEachArray = math.floor((memoryChunk * pow(10.0, 9.0))/(math.ceil(self.N/4)))
		
        numofGenoArray = 0
        genoVecofPointers = []

        if self.M % self.numMarkersofEachArray == 0:
            numofGenoArray = self.M // numMarkersofEachArray
            genoVecofPointers = [(numMarkersofEachArray * math.ceil(self.N / 4)) for _ in range(numofGenoArray)]
        else:
            numofGenoArray = self.M // numMarkersofEachArray + 1
            genoVecofPointers = [(numMarkersofEachArray * math.ceil(self.N / 4)) for _ in range(numofGenoArray - 1)]
            numMarkersofLastArray = self.M - (numofGenoArray - 1) * numMarkersofEachArray
            genoVecofPointers.append((numMarkersofLastArray * math.ceil(self.N / 4)))

        print("size of genoVecofPointers:", len(genoVecofPointers))
        # 捕捉内存分配错误模块尚未添加

        # Setgeno Mark1
        print("setgeno mark1")
        self.alleleFreqVec = np.zeros(self.M)
        self.invstdvVec = np.zeros(self.M)
        self.MACVec = np.zeros(self.M)

        freq = 0.0
        Std = 0.0
        invStd = 0.0
        indexNA = []
        lengthIndexNA = 0
        indexGeno = 0
        indexBit = 0
        fillinMissingGeno = 0
        b2 = 0
        a2 = 0

        ind = 0
        geno1 = 0
        bufferGeno = 0
        u = 0
        
        # setgeno mark2
        print("setgeno mark2")

        for i in range(self.M):
            # Clear and resize the genoVecOneMarkerOld list
            genoVecOneMarkerOld = []

            # Move to the appropriate position in the file and read the data into genoVecOneMarkerOld
            test_bedfile.seek(3 + nbyteOld * i)
            test_bedfile.readinto(genoVecOneMarkerOld)

            # Process the genotypes and fill in missing values
            indexNA = []

            self.Get_OneSNP_Geno_atBeginning(i, indexNA, genoVecOneMarkerOld)

            ind = 0
            geno1 = 0

            for j in range(self.Nnomissing):
                u = j % 4
                bufferGeno = m_OneSNP_Geno[j]
                if bufferGeno == 0:
                    setGenotype(geno1, u, HOM_ALT)
                elif bufferGeno == 1:
                    setGenotype(geno1, u, HET)
                elif bufferGeno == 2:
                    setGenotype(geno1, u, HOM_REF)
                else:
                    setGenotype(geno1, u, MISSING)
                    m_OneSNP_Geno[j] = 0  # 12-18-2017

                if u == 3:
                    genoVecofPointers[i // numMarkersofEachArray].append(geno1)
                    geno1 = 0

            if Nnomissing % 4 != 0:
                genoVecofPointers[i // numMarkersofEachArray].append(geno1)

            lengthIndexNA = len(indexNA)
            freq = np.sum(m_OneSNP_Geno) / (2 * (Nnomissing - lengthIndexNA))





        self.indiceVec = [] # TODO: 每个元素是一个包含两个整数的元组 
        self.kinValueVecFinal = [] # TODO: 放float的列表 

        
        return 0

    def setgeno(self, genofile, subSampleInGeno, MemoryChunk, isDiagofKinSetAsOne):
        start_time = time.time()
        self.setGenoObj(genofile, subSampleInGeno, MemoryChunk, isDiagofKinSetAsOne)        
        stop_time = time.time()
        print("time:", (stop_time - start_time) * 1000, "milliseconds")
        

    def getAlleleFreqVec(self):
        '''
        从geno这个类中获取 allele frequency, 返回形式是一个list/array 
        '''
        return self.AlleleFreqVec

    def setRelatednessCutoff(self, relatednessCutoff):
        '''
        将相关性阈值(relatedness cutoff)设置为一个浮点数
        '''
        self.RelatednessCutoff = relatednessCutoff
        return 0
    
    def print_comb(self, N):
        # 计算组合数的总数量
        x = N * (N - 1) // 2 - 1
        for k in range(x):
            # 计算组合数对应的索引 i 和 j
            i = k // N
            j = k % N
            # 如果 j 小于 i，则交换 i 和 j 的值
            if j < i:
                i = N - i - 2
                j = N - j - 1
            print("i,j:", i, ",", j)

    
    def createSparseKinParallel(self, relatednessCutoff, nblocks = 1, ncore = 1):
        '''
        Create sparse Kinship matrix (GRM) 

        Args:
            relatednessCutoff (float): The threshold to treat two samples as unrelated if IsSparseKin is TRUE. By default, 0.125;
            nblocks (int): The multiple threads number.
            ncore (int): The multiple threads number

        Returns:
            A kin list   
        
        '''

        self.setRelatednessCutoff(relatednessCutoff)
        
        self.Get_MultiMarkersBySample_StdGeno_Mat()  

        #  Record and calculate time usage
        tp0 = time.time()
        self.printComb(3)
        tp1 = time.time()
        print("tp1 - tp0:", tp1 - tp0)

        sparseKinList = self.refineKin(relatednessCutoff)
        
        Nval = self.getNnomissingOut()

        sparseKinList.iIndex = list(sparseKinList.iIndex, [i for i in range(Nval)])
        sparseKinList.jIndex = list(sparseKinList.jIndex, [i for i in range(Nval)])

        diagKin = self.get_DiagofKin()
        sparseKinList.kinValue = list(sparseKinList.kinValue, diagKin)

        tp2 = time.time()
        print("tp2 - tp1: ", tp2-tp1)

        return sparseKinList


    def get_multi_markers_by_sample_std_geno_mat(self):
        '''
        从基因数据中获取多个标记位点的标准化基因型数据，并存储在一个矩阵中
        '''
        # 获取子标记位点的数量
        m_M_Submarker = self.get_sub_marker_num()
        # 获取没有缺失值的个体数量
        Nnomissing = self.get_n_nomissing()
        # 初始化存储标准化基因型数据的矩阵
        std_geno_multi_markers_mat = np.zeros((m_M_Submarker, Nnomissing))
        
        for k in range(m_M_Submarker):
            ind = 0
            flag = 0
            # 获取第 k 个子标记位点的索引
            SNPIdx = self.sub_marker_index[k]
            # 计算 SNP 在 vector 中的索引
            indexOfVectorPointer = SNPIdx // self.num_markers_of_each_array
            # 计算 SNP 在 vector 中的位置
            SNPIdxinVec = SNPIdx % self.num_markers_of_each_array
            # 计算 SNP 在 vector 中的起始位置
            Start_idx = self.m_size_of_esi * SNPIdxinVec
            # 获取频率和逆标准差
            freq = self.allele_freq_vec[SNPIdx]
            invStd = self.inv_stdv_vec[SNPIdx]

            while flag == 0:
                for i in range(Start_idx, Start_idx + self.m_size_of_esi):
                    geno1 = self.geno_vec_of_pointers[indexOfVectorPointer][i]
                    
                    for j in range(4):
                        b = geno1 & 1
                        geno1 >>= 1
                        a = geno1 & 1
                        # 计算标准化基因型数据并存储在矩阵中
                        std_geno_multi_markers_mat[k, ind] = ((2 - (a + b)) - 2 * freq) * invStd
                        ind += 1
                        geno1 >>= 1
                        
                        if ind == Nnomissing:
                            flag = 1
                            break

        return std_geno_multi_markers_mat
    
    def init_kin_value_vec_final(self, ni):
        self.kinValueVecFinal = [0] * ni
        return None
    
    # def get_one_snp_geno(self, SNPIdx):
    #     temp = self.one_snp_geno(SNPIdx)
    #     return temp

    def set_sparse_kin_lookup_arr(self, maf_val, invsd_val):
        '''
        根据给定的 MAF 和逆标准差值设置稀疏核心矩阵查找表。
        '''
        maf_val2 = 2 * maf_val
        a0 = (0 - maf_val2) * invsd_val
        a1 = (1 - maf_val2) * invsd_val
        a2 = (2 - maf_val2) * invsd_val
        
        s_kin_lookup_arr = [[0.0]*3 for _ in range(3)]  # 初始化一个3x3的数组
        
        s_kin_lookup_arr[0][0] = a0 * a0
        s_kin_lookup_arr[0][1] = a0 * a1
        s_kin_lookup_arr[0][2] = a0 * a2
        s_kin_lookup_arr[1][0] = s_kin_lookup_arr[0][1]
        s_kin_lookup_arr[1][1] = a1 * a1
        s_kin_lookup_arr[1][2] = a1 * a2
        s_kin_lookup_arr[2][0] = s_kin_lookup_arr[0][2]
        s_kin_lookup_arr[2][1] = s_kin_lookup_arr[1][2]
        s_kin_lookup_arr[2][2] = a2 * a2

        return s_kin_lookup_arr
    
    def parallelcalsparseGRM(self, GRMvec):
        '''
        用于并行计算稀疏遗传相关矩阵(GRM)并将结果存储在 GRMvec 中
        '''
        return None

    def parallelsumTwoVec(self, x):
        n1 = len(x)
        # allocate the output matrix
        sumVec = np.zeros(n1, dtype=np.float32)

        # function call operator that work for the specified range (begin/end)
        for i in range(n1):
            sumVec[i] = x[i]

        return sumVec


    def refine_kin(self, relatedness_cutoff):
        '''
        用于计算并返回满足相关性阈值要求的 kinship 列表
        '''
        i_index_vec2 = []
        j_index_vec2 = []
        kin_value_vec2 = []

        temp = geno.m_OneSNP_StdGeno
        temp.clear()
        ni = len(geno.indiceVec)
        print("ni:", ni)

        self.init_kin_value_vec_final(ni)

        M_marker = self.M
        M_marker_maf_gr1perc = 0

        for i in range(M_marker):
            freq_v = self.alleleFreqVec[i]
            minimum_MAF_to_construct_GRM = 0.001
            if minimum_MAF_to_construct_GRM <= freq_v <= (1 - minimum_MAF_to_construct_GRM):
                M_marker_maf_gr1perc += 1
                # self.Get_OneSNP_Geno(i)
                inv_std_v = geno.invstdvVec[i]
                self.s_kin_lookup_arr = self.setSparseKinLookUpArr(freq_v, inv_std_v)
                GRM_vec = self.parallelcalsparseGRM() # 这步不知道在干嘛
                sumVec = self.parallelsumTwoVec(GRM_vec)

        for j in range(ni):
            self.kinValueVecFinal[j] /= M_marker_maf_gr1perc
            if self.kinValueVecFinal[j] >= relatedness_cutoff:
                a1 = self.indiceVec[j][0] + 1
                a2 = self.indiceVec[j][1] + 1
                i_index_vec2.append(a1)
                j_index_vec2.append(a2)
                kin_value_vec2.append(self.kinValueVecFinal[j])

        print("kinValueVec2.size():", len(kin_value_vec2))
        return {"iIndex": i_index_vec2, "jIndex": j_index_vec2, "kinValue": kin_value_vec2}
    

    def getNnomissingOut(self):
        # 这个函数是干什么的呀！
        return self.genoObj['Nnomissing']
    
    def get_DiagofKin(self):
        M = self.genoObj['M']
        Nnomissing = self.genoObj['Nnomissing']
        
        x = np.zeros(Nnomissing)
        
        if not self.genoObj['setKinDiagtoOne']:
            x = self.getnumberofMarkerswithMAFge_minMAFtoConstructGRM
        else:
            x = np.ones(Nnomissing)
            
        return x
    









