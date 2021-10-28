argv <- commandArgs(T)
VCF_INFO<- argv[1]

library(dplyr)
library(stringr)
#read the vcf information file
vcf_info<-read.table(VCF_INFO, header=F, stringsAsFactor=F)
colnames(vcf_info)<-c("CHR","start","id","end","startQ","endQ","REF","ALT")
vcf_info[1:2,]

#calculate length
vcf_info$LEN_ref<-vcf_info$end-vcf_info$start
vcf_info$LEN_qry<-vcf_info$endQ-vcf_info$startQ
vcf_info$max_LEN<-pmax(vcf_info$LEN_qry,vcf_info$LEN_ref)

#categorize as INS or DEL relatively to reference
vcf_info$TYPE<-substr(vcf_info$ALT,2,4)
vcf_info$TYPE[vcf_info$TYPE=="CPG"]<-"INS"
vcf_info$TYPE[vcf_info$TYPE=="CPL"]<-"DEL"
head(vcf_info,20)

write.table(vcf_info[,c(1,2,11,12)],paste0(VCF_INFO,".annot"),sep="\t", row.names=F, quote=F, col.names=F)


