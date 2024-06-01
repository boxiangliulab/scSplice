###input: non-split-reads; intron retention reads
library(Rsubread)
library(Rsamtools)
library(GenomicRanges)
library(data.table)
suppressMessages(library(argparser, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(foreach, quietly=TRUE))
suppressMessages(library(data.table, quietly=TRUE))

IR_compute<-function(input_dir, input_prefix, input_perind, output_dir,output_prefix){
non_split_reads<-read.table(paste(input_dir,input_prefix,sep="/"),sep=" ",check.names=FALSE)
for(i in 1:ncol(non_split_reads)){
    tmp<-strsplit(colnames(non_split_reads)[i],"/")[[1]]
    l_tmp<-length(tmp)
    colnames(non_split_reads)[i]<-tmp[l_tmp]
}
total<-fread(input_perind,sep=" ")
total<-as.data.frame(total)

total$chr<-NA
total$start<-NA
total$end<-NA
total$clu<-NA
total$strand<-NA
for(i in 1:nrow(total)){
    total$chr[i]<-strsplit(total$V1[i],":")[[1]][1]
    total$start[i]<-strsplit(total$V1[i],":")[[1]][2]
    total$end[i]<-strsplit(total$V1[i],":")[[1]][3]
    clu<-strsplit(total$V1[i],":")[[1]][4]
    total$clu[i]<-clu
    total$strand[i]<-strsplit(clu,"_")[[1]][3]
}

tmp<-as.data.frame(matrix(NA,0,ncol(total)-6))
for(i in 1:nrow(total)){
    clu<-total$clu[i]
    num<-which(total$clu==clu)
    if(length(num)==1){
        start_site<-paste(total$chr[i],total$strand[i],total$start[i],sep=":")
        end_site<-paste(total$chr[i],total$strand[i],total$end[i],sep=":")
        num_1<-which(rownames(non_split_reads)==start_site)
        if(length(num_1)==1){
            tmp_one<-as.data.frame(matrix(NA,1,ncol(total)-6))
            rownames(tmp_one)<-start_site
            for(k in 1:(ncol(total)-6)){
                IR<-total[i,k+1]
                nsr<-non_split_reads[num_1,colnames(total)[k+1]]
                ##tmp_one[1,k]<-IR/(IR+nsr)
                tmp_one[1,k]<-paste(IR,IR+nsr,sep=":")
            }
            tmp<-rbind(tmp,tmp_one)
        }
        num_2<-which(rownames(non_split_reads)==end_site)
        if(length(num_2)==1){
            tmp_one<-as.data.frame(matrix(NA,1,ncol(total)-6))
            rownames(tmp_one)<-end_site
            for(k in 1:(ncol(total)-6)){
                IR<-total[i,k+1]
                nsr<-non_split_reads[num_2,colnames(total)[k+1]]
                tmp_one[1,k]<-paste(IR,IR+nsr,sep=":")
                ##tmp_one[1,k]<-IR/(IR+nsr)
            }
            tmp<-rbind(tmp,tmp_one)
        }
    }
}

colnames(tmp)<-colnames(total)[2:(ncol(total)-5)]
write.table(tmp,paste(output_dir,output_prefix,sep="/"),sep=" ",row.names=TRUE,col.names=TRUE,quote=FALSE)

}

p <- arg_parser("scSplice: Identify IR site")
p <- add_argument(p, "--input_dir", help="Deposit position of non split read counts",default=".")
p <- add_argument(p, "--input_prefix", help="File containing non-split reads")
p <- add_argument(p, "--input_perind", help="File listing perind_numers.constcounts.gz")
p <- add_argument(p, "--output_dir", help="Output position",default=".")
p <- add_argument(p, "--output_prefix", help="Output prefix")
argv <- parse_args(p)

cat("scSplice: start computing IR ratios\n")

m <- IR_compute(argv$input_dir, argv$input_prefix, argv$input_perind, argv$output_dir,argv$output_prefix)




