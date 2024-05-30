from SPA import *
from Step1 import *
from Prepare_data import *

def ScoreTest(geno_PATH, null_model_result):
    G = read_geno(geno_PATH)
    result = TestSPA(G, null_model_result)
    return result