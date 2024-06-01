suppressMessages(library(argparser, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(foreach, quietly=TRUE))
suppressMessages(library(data.table, quietly=TRUE))
library(SPAtest)

cauchy <- function(p_values){
    transformed_sum = sum(tan(pi*(p_values - 0.5)))
    combined_p_value<-0.5+(atan(transformed_sum)/pi)
    return(combined_p_value)
}
compare <- function(individual,input_dir,input_prefix){
    ind<-fread(individual,sep=" ",header=FALSE)
    PSI<-fread(paste(input_dir,"/",input_prefix,sep=""),sep=" ")
    PSI<-as.data.frame(PSI)
    ind<-as.data.frame(ind)
    term<-ind$V2[!duplicated(ind$V2)]
    G1<-which(ind[,2]==term[1])
    G2<-which(ind[,2]==term[2])
    sig<-as.data.frame(matrix(NA,0,2))
    for(i in 1:nrow(PSI)){
        pheno<-c(as.numeric(PSI[i,ind[G1,1]]),as.numeric(PSI[i,ind[G2,1]]))
        geno<-c(rep(0,length(G1)),rep(2,length(G2)))
        pval<-t.test(as.numeric(PSI[i,ind[G1,1]]),as.numeric(PSI[i,ind[G2,1]]))$p.value
        #pval<-ScoreTest_SPA(geno,pheno,method="fastSPA")$p.value
        newsig<-as.data.frame(matrix(NA,1,2))
        newsig[1,1]<-PSI[i,1]
        newsig[1,2]<-pval
        sig<-rbind(sig,newsig)
    }
    predict_gene<-as.data.frame(matrix(NA,nrow(sig),1))
    clu_gene<-fread("/data/zhangyuntian/simulation/PSI_centric/test/test_3_31_clu_gene.txt",sep="\t",header=FALSE)
    clu_gene<-as.data.frame(clu_gene)
    for(i in 1:nrow(sig)){
        tmp<-strsplit(sig[i,1],":")[[1]]
        clu_tmp<-paste(tmp[1],tmp[2],sep=":")
        num<-which(clu_gene[,1]==clu_tmp)
        if(length(num)>0){
            predict_gene[i,1]<-clu_gene[num[1],2]
            predict_gene[i,2]<-sig[i,2]
        }
    }
    output<-as.data.frame(matrix(NA,300,3))
    ground_truth_gene<-fread("/data/zhangyuntian/simulation/code/ground_truth_gene",header=FALSE)
    ground_truth_gene<-as.data.frame(ground_truth_gene)
    total<-fread("/data/zhangyuntian/simulation/code/simu_singlecell_3_8/total_gene",header=FALSE)
    total<-as.data.frame(total)
    output$V1<-total$V1
    for(i in 1:nrow(output)){
        num<-which(ground_truth_gene[,1]==total[i,1])
        if(length(num)>0){
            output[i,2]<-1
            tmp<-which(predict_gene[,1]==output[i,1])
            if(length(tmp)>0){
                output[i,3]<-1-cauchy(as.numeric(predict_gene[tmp,2]))
                ##adj.p<-p.adjust(as.numeric(predict_gene[tmp,2]),method="bonferroni",n=length(tmp))
                ##output[i,3]<-max(1-adj.p)
                ##output[i,3]<-max(1-as.numeric(predict_gene[tmp,2]))
            }
            else{
                output[i,3]<-0
            }
        }
        else{
            output[i,2]<-0
            tmp<-which(predict_gene[,1]==output[i,1])
            if(length(tmp)>0){
                output[i,3]<-1-cauchy(as.numeric(predict_gene[tmp,2]))
                ##adj.p<-p.adjust(as.numeric(predict_gene[tmp,2]),method="bonferroni",n=length(tmp))
                ##output[i,3]<-max(1-adj.p)
                ##output[i,3]<-max(1-as.numeric(predict_gene[tmp,2]))
            }
            else{
                output[i,3]<-0
            }
        }
    }
    output
}

p <- arg_parser("scSplice: Do differential splicing analysis")
p <- add_argument(p, "--individual", help="File which choose the compare object",default=".")
p <- add_argument(p, "--input_dir", help="Input position")
p <- add_argument(p, "--input_prefix", help="File listing all PSI of each splicing event")
p <- add_argument(p, "--output_dir", help="Output position",default=".")
p <- add_argument(p, "--output_prefix", help="Output prefix")

argv <- parse_args(p)

cat("scSplice: start DS analyses\n")

m <- compare(argv$individual, argv$input_dir,argv$input_prefix)
write.table(m, file.path(argv$output_dir,"/", argv$output_prefix), sep = "\t", quote=FALSE, row.names=FALSE,col.names=FALSE)
