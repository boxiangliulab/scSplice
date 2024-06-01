suppressMessages(library(argparser, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(foreach, quietly=TRUE))
suppressMessages(library(data.table, quietly=TRUE))

phenotype<-function(input_dir, input_prefix_clu,input_prefix_IR, cluname){
    clutmp<-fread(cluname,sep=" ")
    clutmp<-as.data.frame(clutmp)
    clu_info<-as.data.frame(matrix(NA,nrow(clutmp),2))
    for(i in 1:nrow(clu_info)){
        tmp<-strsplit(clutmp[i,1],":")[[1]][4]
        clu_info[i,1]<-tmp
        clu_info[i,2]<-strsplit(tmp,"_")[[1]][2]
    }
    ind<-colnames(clutmp)
    total_PSI<-as.data.frame(matrix(NA,0,length(ind)))
    colnames(total_PSI)<-ind
    file_path<-paste(input_dir,"/",input_prefix_clu,sep="")
    if(file.info(file_path)$size>0){
        PSI<-fread(paste(input_dir,"/",input_prefix_clu,sep=""),sep=" ")
        PSI<-as.data.frame(PSI)
        colnames(PSI)<-ind
        for(i in 1:nrow(PSI)){
            ###delete all the events losing in more than 1-prop of all individuals
            chr<-strsplit(PSI[i,1],":")[[1]][1]
            mid<-strsplit(PSI[i,1],":")[[1]][2]
            strand<-strsplit(mid,"_")[[1]][3]
            site<-strsplit(PSI[i,1],":")[[1]][3]
            PSI[i,1]<-paste(chr,strand,site,sep=":")
        }
        PSI<-na.omit(PSI)
    }
    file_path<-paste(input_dir,"/",input_prefix_IR,sep="")
    if(file.info(file_path)$size>0){
        IR<-fread(paste(input_dir,"/",input_prefix_IR,sep=""),sep=" ")
        IR<-as.data.frame(IR)
        IR<-na.omit(IR)
    }
    ###combine IR & PSI
    total<-rbind(PSI,IR)
    total
}

p <- arg_parser("scSplice: Prepare phenotype file from PSI results")
p <- add_argument(p, "--input_dir", help="Deposit position of PSI value",default=".")
p <- add_argument(p, "--input_prefix_clu", help="File listing all unique sites in clusters")
p <- add_argument(p, "--input_prefix_IR", help="File listing all unique sites in intron retention")
p <- add_argument(p, "--output_dir", help="Output position",default=".")
p <- add_argument(p, "--output_prefix", help="Output prefix")
p <- add_argument(p, "--cluname", help="Input the perind.counts.gz file")
argv <- parse_args(p)

cat("scSplice: start read in PSI values\n")

m <- phenotype(argv$input_dir, argv$input_prefix_clu,argv$input_prefix_IR, argv$cluname)
write.table(m, file.path(argv$output_dir,"/", argv$output_prefix), sep = " ", quote=FALSE, row.names=FALSE)