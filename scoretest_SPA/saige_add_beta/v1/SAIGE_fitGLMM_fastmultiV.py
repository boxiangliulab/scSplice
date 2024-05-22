import os
import pandas as pd
import numpy as np
# import statsmodels.api as sm



def fitNULLGLMM_multiV(plinkFile = "",
                       phenoFile = "",
                       config = {}):
    
    # Set up output files
    modelOut = f"{outputPrefix}.rda"
    skipModelFitting = config["skipModelFitting"]
    if skipModelFitting:
        if not os.path.exists(modelOut):
            raise Exception(f"skipModelFitting=TRUE but {modelOut} does not exist")
    else:
        LOCO = config["LOCO"]
        isLowMemLOCO = config["isLowMemLOCO"]
        if LOCO and isLowMemLOCO:
            modelOut = f"{outputPrefix}_noLOCO.rda"
        os.makedirs(modelOut, exist_ok=True)
    
    # You'll need to read in the plink files, phenotype files, etc. and perform
    # the appropriate analysis
    plinkFile = config["plinkFile"]
    if plinkFile != "":
        bimFile = f"{plinkFile}.bim"
        bedFile = f"{plinkFile}.bed"
        famFile = f"{plinkFile}.fam"

    # TODO: setgenoNULL()是什么
    def setgenoNULL():
        # Placeholder for the actual functionality of 'setgenoNULL' from the R code.

        pass

    setgenoNULL()

    useGRMtoFitNULL = config["useGRMtoFitNULL"]
    if not useGRMtoFitNULL:
        useSparseGRMtoFitNULL = False
        useSparseGRMforVarRatio = False
        LOCO = False
        nThreads = 1
        print("No GRM will be used to fit the NULL model and nThreads is set to 1")

    if useSparseGRMtoFitNULL and bedFile == "":
        print("Sparse GRM is used to fit the null model and plink file is not specified, so variance ratios won't be estimated")
        skipVarianceRatioEstimation = True

    numMarkersForVarRatio = config["numMarkersForVarRatio"]
    IsOverwriteVarianceRatioFile = config["IsOverwriteVarianceRatioFile"]
    if not skipVarianceRatioEstimation:
        SPAGMMATOut = f"{outputPrefix}_{numMarkersForVarRatio}markers.SAIGE.results.txt"
        # Equivalent functionality for 'Check_OutputFile_Create' should be implemented here if necessary

        if outputPrefix_varRatio == "":
            outputPrefix_varRatio = outputPrefix

        varRatioFile = f"{outputPrefix_varRatio}.varianceRatio.txt"

        if not os.path.exists(varRatioFile):
            # Create the file. In Python, merely opening a file in write mode will create it.
            with open(varRatioFile, 'w') as file:
                pass  # File is created and immediately closed
        else:
            if not IsOverwriteVarianceRatioFile:
                raise Exception(f"WARNING: The variance ratio file {varRatioFile} already exists. The new variance ratios will be output to {varRatioFile}. In order to avoid overwriting the file, please remove the {varRatioFile} or use the argument outputPrefix_varRatio to specify a different prefix to output the variance ratio(s). Otherwise, specify --IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten with new variance ratio(s)")
            else:
                print(f"The variance ratio file {varRatioFile} already exists. IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten")
    else:
        print("Variance ratio estimation will be skipped.")
        useSparseGRMforVarRatio = False

    if useSparseGRMtoFitNULL:
        # useSparseGRMforVarRatio = False  # This line is commented out in the R code, so it's omitted here.
        LOCO = False
        nThreads = 1
        if bedFile != "":
            print("sparse GRM will be used to fit the NULL model and nThreads is set to 1")
        print("Leave-one-chromosome-out is not applied")

    # Assuming nThreads, FemaleOnly, MaleOnly, outputPrefix, FemaleCode, MaleCode, and sexCol are already defined
    # TODO：并行计算 -> multiprocessing/concurrent.futures/joblib/Dask
    if nThreads > 1:
        # In Python, you might use multiprocessing or a similar module depending on your needs
        # RcppParallel::setThreadOptions is specific to R and doesn't have a direct Python equivalent
        print(f"{nThreads} threads will be used")

    FemaleOnly = config["FemaleOnly"]
    MaleOnly = config["MaleOnly"]
    FemaleCode = config["FemaleCode"]
    MaleCode = config["MaleCode"]
    sexCol = config["sexCol"]
    if FemaleOnly and MaleOnly:
        raise Exception("Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job")

    if FemaleOnly:
        outputPrefix += "_FemaleOnly"
        print(f"Female-specific model will be fitted. Samples coded as {FemaleCode} in the column {sexCol} in the phenotype file will be included")

    elif MaleOnly:
        outputPrefix += "_MaleOnly"
        print(f"Male-specific model will be fitted. Samples coded as {MaleCode} in the column {sexCol} in the phenotype file will be included")



if __name__ == "__main__":
    import json

    with open("nullGLM_config.json",'r') as load_f:
        config_nullGLMM =json.load(load_f)
        print(config_nullGLMM)
 
    
    fitNULLGLMM_multiV(plinkFile="",
                       phenoFile="", 
                       config=config_nullGLMM)

