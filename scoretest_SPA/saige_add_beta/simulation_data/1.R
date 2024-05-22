


pheno_data <- read.table("simu_scSplice_pheno.txt", header=TRUE, check.names=FALSE)
pheno <- data.frame(pheno_data[,-1])
pheno1 <- pheno[1,]
pheno1 <- as.numeric(pheno1)

pheno1_minmax <- (pheno1 - min(pheno1)) / (max(pheno1) - min(pheno1))





convert_genotype <- function(genotype) {
  splits <- strsplit(genotype, split="/|\\|")
  counts <- sapply(splits, function(x) sum(as.numeric(x)))
  return(counts)
}
library(vcfR)
vcf <- read.vcfR("simu_geno.vcf.gz")
genos <- extract.gt(vcf)
#genos[0,101:500]
geno1_matrix <- matrix(genos[64,101:500],nrow = 1)
geno1 <- matrix(apply(geno1_matrix, 1:2, convert_genotype), nrow = 1)










cov_data <- read.table("simu_covariate.txt", header=TRUE, check.names=FALSE)
cov <- data.frame(cov_data)
cov1 <- cov[1,102:501]
cov1 <- as.numeric(cov1)




# library("SPAtest")
#result <- ScoreTest_SPA(geno1, pheno1, cov1, method=c("fastSPA"), minmac=5, Cutoff=2, alpha=5*10^-8, missing.id=NA,beta.out=FALSE, beta.Cutoff=5*10^-7, log.p=FALSE)
result <- ScoreTest_SPA(as.numeric(geno1), pheno1_minmax,cov1, method=c("fastSPA"))

print(result)








# 
# 
# library(vcfR)
# library(SPAtest)
# 
# # 读取表型数据和基因型数据
# pheno_data <- read.table("simu_scSplice_pheno.txt", header=TRUE, check.names=FALSE)
# pheno <- data.frame(pheno_data[,-1])
# 
# vcf <- read.vcfR("simu_geno.vcf.gz")
# genos <- extract.gt(vcf)
# 
# # 读取协变量数据
# cov_data <- read.table("simu_covariate.txt", header=TRUE, check.names=FALSE)
# cov <- data.frame(cov_data)
# 
# # 函数用于转换基因型数据
# convert_genotype <- function(genotype) {
#   splits <- strsplit(genotype, split="/|\\|")
#   counts <- sapply(splits, function(x) sum(as.numeric(x)))
#   return(counts)
# }
# 
# # 准备存储p值的矩阵
# p_values <- matrix(nrow = 80, ncol = 1000)
# 
# # 循环处理前80个表型和前1000个基因型
# for (i in 1:80) {
#   pheno1 <- as.numeric(pheno[i,])
#   pheno1_minmax <- (pheno1 - min(pheno1)) / (max(pheno1) - min(pheno1))
#   
#   for (j in 1:1000) {
#     geno1_matrix <- matrix(genos[j,101:500], nrow = 1)
#     geno1 <- matrix(apply(geno1_matrix, 1:2, convert_genotype), nrow = 1)
#     cov1 <- as.numeric(cov[1,102:501])
#     
#     # 运行SPA得分测试
#     result <- ScoreTest_SPA(as.numeric(geno1), pheno1_minmax, cov1, method=c("fastSPA"))
#     p_values[i, j] <- result$p.value
#   }
# }
# 
# # 将结果保存为CSV文件
# write.csv(p_values, "p_values.csv")
# 
# # 寻找所有小于0.05的p值的索引
# significant_indices <- which(p_values < 0.05, arr.ind = TRUE)
# 
# # 创建一个包含这些索引的数据框
# significant_results <- data.frame(
#   Row = significant_indices[,1],
#   Column = significant_indices[,2],
#   P_Value = p_values[significant_indices]
# )
# 
# # 将结果保存为CSV文件
# write.csv(significant_results, "significant_p_values.csv")
# 


































