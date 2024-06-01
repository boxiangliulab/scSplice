###this script used to extract nonsplit-reads for each 5' splice site
###Intron retention: Junction reads/ (Junction reads + nonsplite-reads)
library(Rsubread)
library(Rsamtools)
library(GenomicRanges)
suppressMessages(library(argparser, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(foreach, quietly=TRUE))
suppressMessages(library(data.table, quietly=TRUE))
###R: read in intron retention list
IR_identify<-function(input_dir, input_prefix, input_bam, output_dir,output_prefix,read_len){
   intron_retention=paste(input_dir,input_prefix,sep="/")
   splice_site=list()
   lines<-readLines(intron_retention)
   for(line in lines){
      chr = strsplit(line," ")[[1]][1]
      if(length(strsplit(line," ")[[1]])==2){
         intron=strsplit(line," ")[[1]][2]
         start=strsplit(intron,":")[[1]][1]
         end=strsplit(intron,":")[[1]][2]
         splice_site[length(splice_site)+1]=paste(chr,start,sep=":")
         splice_site[length(splice_site)+1]=paste(chr,end,sep=":")
      }
   }
   splice_site<-as.character(splice_site)
   dedup_splicesite<-unique(splice_site)
  ###read in bamfiles needed to be considered
  bamfile<-read.table(input_bam)
newIR_site<-as.data.frame(matrix(NA,length(dedup_splicesite),nrow(bamfile)))
rownames(newIR_site)<-dedup_splicesite
colnames(newIR_site)<-bamfile[,1]

for(i in 1:nrow(bamfile)){
   bam <- scanBam(bamfile[i,1],sep="")
   newdat<-as.data.frame(matrix(NA,length(bam[[1]]$cigar),0))
   newdat$rname<-bam[[1]]$rname
   newdat$pos<-bam[[1]]$pos
   newdat$cigar<-bam[[1]]$cigar
   non_split_reads<-newdat[which(grepl("N",newdat$cigar)!=TRUE),]
# Convert your non-split reads to a GRanges object
   read_ranges <- GRanges(seqnames = non_split_reads$rname,
                       ranges = IRanges(start = non_split_reads$pos,
                                        end = non_split_reads$pos + as.numeric(read_len)-1 ))
# Find overlaps
   for(k in 1:length(dedup_splicesite)){
      chr=strsplit(dedup_splicesite[k],":")[[1]][1]
      site=strsplit(dedup_splicesite[k],":")[[1]][3]
      splice_sites <- GRanges(chr, IRanges(start = as.numeric(site)-1, end = as.numeric(site)+1))
      overlaps <- findOverlaps(read_ranges, splice_sites)
# Count overlaps
      newIR_site[k,i]<-length(unique(queryHits(overlaps)))}
}
    newIR_site
}

p <- arg_parser("scSplice: Identify IR site")
p <- add_argument(p, "--input_dir", help="Deposit position of all sites",default=".")
p <- add_argument(p, "--input_prefix", help="File listing cluster's sites")
p <- add_argument(p, "--input_bam", help="File listing all bamfiles need to be considered")
p <- add_argument(p, "--output_dir", help="Output position",default=".")
p <- add_argument(p, "--output_prefix", help="Output prefix")
p <- add_argument(p, "--read_len", help="read length; default=151")
argv <- parse_args(p)

cat("scSplice: start identifying IR site\n")

m <- IR_identify(argv$input_dir, argv$input_prefix, argv$input_bam, argv$output_dir,argv$output_prefix,argv$read_len)
write.table(m, file.path(argv$output_dir,"/", argv$output_prefix), sep = " ", quote=FALSE)