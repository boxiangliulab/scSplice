import time
import pandas as pd
import numpy as np

from SAIGE_geno import geno



def createSparseGRM(plinkFile = "",
                    outputPrefix = "",
                    numRandomMarkerforSparseKin = 1000,
                    relatednessCutoff = 0.125,
                    memoryChunk = 2,
                    isDiagofKinSetAsOne = False,
                    nThreads = 1,
                    minMAFforGRM = 0.01,
                    isSetGeno=True,
                    isWritetoFiles=True):
    '''
    Construct a sparse GRM for a given data set 

    Args:
        plinkFile (chr): Path to plink file to be used for calculating the sparse GRM;
        outputPrefix (chr): Path to the output files with prefix;
        numRandomMarkerforSparseKin (int): number of randomly selected markers (MAF >= 0.01) to be used to identify related samples for sparse GRM. By default, 1000;
        relatednessCutoff (float): The threshold to treat two samples as unrelated if IsSparseKin is TRUE. By default, 0.125;
        memoryChunk (integer or float): The size (Gb) for each memory chunk. By default, 2;
        isDiagofKinSetAsOne  (logical): Whether to set the diagnal elements in GRM to be 1. By default, FALSE;
        nThreads (integer): Number of threads to be used. By default, 1;
        minMAFforGRM (list): Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01.

    Returns:
        a file ended with sampleIDs.txt that contains sample IDs for the sparse GRM and a file ended with .sparseGRM.mtx that contains the sparse GRM.

    '''

    # Set the multiple threads for computing
    if nThreads > 1:
        # wait for multi-computing
        print(nThreads, " threads are set to be used.")

    # Set sparse GRM
    if minMAFforGRM > 0:
        print("Markers in the Plink file with MAF >= ", minMAFforGRM, " will be used to construct GRM")
    else:
        print("Markers in the Plink file with MAF > ", minMAFforGRM, " will be used to construct GRM")

    # Begin to read in the plinkFile
    famFile = plinkFile + ".fam"
    fam = pd.read_csv(famFile, header=None, dtype=str)
    fam.columns = ["character" for i in range(4)] + ["numeric" for i in range(2)]
    fam.iloc[:, 4:] = fam.iloc[:, 4:].astype(float)

    sparseGRMSampleID = fam[:, 2]
    sparseGRMSampleIDFile = print(outputPrefix,"_relatednessCutoff_",relatednessCutoff,"_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")

    # If write them to files
    if isWritetoFiles:
        print("write sample IDs for the sparse GRM to ", sparseGRMSampleIDFile)
        sparseGRMSampleID.to_csv(sparseGRMSampleIDFile, sep='\t', quoting=False, header=False, index=False)

    genoSampleIndex = list(range(1, fam.shape[0]))

    # Whether need to setGeno
    if isSetGeno:
        geno_instance = geno()
        geno_instance.setgeno(plinkFile, genoSampleIndex, memoryChunk, isDiagofKinSetAsOne)
    
    freqVec = geno_instance.getAlleleFreqVec() # TODO: try cpp to finish this function

    if minMAFforGRM > 0:
        MAFindex = [i for i, freq in enumerate(freqVec) if minMAFforGRM <= freq <= 1 - minMAFforGRM]
        print(len(MAFindex), " markers have MAF >= ", minMAFforGRM)
    else:
        MAFindex = [i for i, freq in enumerate(freqVec) if 0 < freq < 1]
        print(len(MAFindex), " markers have MAF > ", minMAFforGRM)

    print(numRandomMarkerforSparseKin, " genetic markers are randomly selected to decide which samples are related")

    if len(MAFindex) < numRandomMarkerforSparseKin:
        if minMAFforGRM > 0:
            raise ValueError("ERROR! not enough genetic markers with MAF >= ", minMAFforGRM,
                            " to detect which samples are related\nTry include at least ", numRandomMarkerforSparseKin,
                            " genetic markers with MAF >= ", minMAFforGRM, " in the plink file")
        else:
            raise ValueError("ERROR! not enough genetic markers with MAF > ", minMAFforGRM,
                            " to detect which samples are related\nTry include at least ", numRandomMarkerforSparseKin,
                            " genetic markers with MAF > ", minMAFforGRM, " in the plink file")

    # Start detecting related samples for the sparse GRM
    # print("Start detecting related samples for the sparse GRM")
    # ta = time.time()
    # setSubMarkerIndex(markerIndexforSparseM - 1)
    # tb = time.time()
    # print("Time Used:", tb - ta)
            
    sparseMList = geno_instance.createSparseKinParallel(relatednessCutoff, nblocks = nThreads, ncore = nThreads)
    
    print("length(sparseMList$iIndex): ", len(sparseMList['iIndex']), "\n")
    print(sparseMList['iIndex'][:102])
    print("length(sparseMList$jIndex): ", len(sparseMList['jIndex']), "\n")
    print(sparseMList['jIndex'][:102])
    print("length(sparseMList$kinValue): ", len(sparseMList['kinValue']), "\n")
    print(sparseMList['kinValue'][:102])

    from scipy.sparse import csr_matrix
    import numpy as np

    i = np.array(sparseMList['iIndex'])
    j = np.array(sparseMList['jIndex'])
    x = np.array(sparseMList['kinValue'])

    sparseGRM = csr_matrix((x, (i, j)), shape=(max(i)+1, max(j)+1))
    print("nrow(sparseGRM): ", sparseGRM.shape[0], "\n")
    print("ncol(sparseGRM): ", sparseGRM.shape[1], "\n")
    print("Number of non-zero entries in sparseGRM: ", sparseGRM.nnz, "\n")


    return [sparseGRMSampleID, sparseGRM]