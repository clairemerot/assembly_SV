argv <- commandArgs(T)
REF_FASTA<-argv[1]
QRY_FASTA<-argv[2]
INPUT<-argv[3]
OUTPUT<-argv[3]





library(dplyr)
library(data.table)

vcf<-fread(INPUT, header=F, stringsAsFactor=F)
colnames(vcf)<-c("CHR","POS","id_SV","old_ref","old_alt","QUAL","FILTER","INFO","FORMAT","IND")
head(vcf)

ref<-fread(REF_FASTA, header=F, stringsAsFactor=F)
colnames(ref)<-c("id_SV","CHR","POS","END","strand","REF")
head(ref)

qry<-fread(QRY_FASTA, header=F, stringsAsFactor=F)
colnames(qry)<-c("id_SV","contig_qry","pos_qry","end_qry","strand_qry","ALT")
head(qry)

REF_ALT<-left_join(ref,qry)
head(REF_ALT)
dim(REF_ALT)
dim(ref)
dim(qry)

vcf_annot<-left_join(vcf, REF_ALT)
dim(vcf) #vcf has less SV than the annotation, I don't know why, I guess because of filter at SURVIVOR step
dim(vcf_annot)
head(vcf_annot,1)
vcf_new<-vcf_annot[,c(1,2,3,13,18, 6:10)]

write.table(vcf_new,OUTPUT,sep="\t", row.names=F, col.names=F, quote=F)

