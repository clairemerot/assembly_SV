argv <- commandArgs(T)
VCF_INFO<- argv[1]

library(dplyr)

#read the vcf information file
vcf_info<-read.table(VCF_INFO, header=F, stringsAsFactor=F)
colnames(vcf_info)<-c("CHR","start","id","end","REF","ALT")
vcf_info[1:2,]

#count the number of base in ref and alt
vcf_info$ref_length<-nchar(vcf_info$REF)
vcf_info$alt_length<-nchar(vcf_info$ALT)
#calculate length
vcf_info$LEN<-vcf_info$alt_length-vcf_info$ref_length

#categorize as INS or DEL relatively to reference
vcf_info$TYPE<-"SNP"
vcf_info$TYPE[vcf_info$LEN<0]<-"DEL"
vcf_info$TYPE[vcf_info$LEN>0]<-"INS"
head(vcf_info,20)

write.table(vcf_info[,c(1,2,9,10)],paste0(VCF_INFO,".annot"),sep="\t", row.names=F, quote=F, col.names=F)


