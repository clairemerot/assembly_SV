argv <- commandArgs(T)
VCF_INFO<- argv[1]
BED_FILE<- argv[2]
ANNOT_OUTPUT<- argv[3]

library(dplyr)

#read the vcf information file
vcf_info<-read.table(VCF_INFO, header=F, stringsAsFactor=F)
colnames(vcf_info)<-c("CHR","start","REF","ALT")
vcf_info[1:2,]

#count the number of base in ref and alt
vcf_info$ref_length<-nchar(vcf_info$REF)
vcf_info$alt_length<-nchar(vcf_info$ALT)
#calculate length
vcf_info$LEN<-vcf_info$alt_length-vcf_info$ref_length
vcf_info$max_LEN<-pmax(vcf_info$alt_length,vcf_info$ref_length)
#categorize as INS or DEL relatively to reference
vcf_info$TYPE<-"INS"
vcf_info$TYPE[vcf_info$LEN<0]<-"DEL"
#now we have the regular indels that have just 1 bp in the ref or the alt, and the complex SV which are diff paths in the graph
#not sure whether to keep them
#in any case their start pos is not the same in the bed file so we need to corect for that to extract the stop
vcf_info$SV_CLASS<-"complex"
vcf_info$SV_CLASS[which(vcf_info$ref_length==1 | vcf_info$alt_length==1)]<-"regular"



head(vcf_info[,c(1,2,5:9)])

bed_info<-read.table(BED_FILE, header=F, stringsAsFactor=F)[,c(1,2,3,4,6,8)]
head(bed_info)
colnames(bed_info)<-c("CHR","start_bed","END", "index","inv","len_bed")
bed_info$start<-bed_info$start_bed
bed_info$start[bed_info$index==4]<-(bed_info$start_bed[bed_info$index==4])+1
head(bed_info)
vcf_bed_info<-left_join(vcf_info, bed_info)
dim(vcf_info)
dim(bed_info)#y en a plus - peut-être que certains ont été éliminé par le graph pour non diploidie
dim (vcf_bed_info)

head(vcf_bed_info[,c(1,2,5:13)])

vcf_bed_info$TYPE[vcf_bed_info$inv==1]<-"INV"
vcf_bed_info$LEN[vcf_bed_info$inv==1]<-vcf_bed_info$ref_length[vcf_bed_info$inv==1]
vcf_bed_info$SV_CLASS[vcf_bed_info$inv==1]<-"regular"
head(vcf_bed_info[vcf_bed_info$inv==1,c(1,2,5:13)])
head(vcf_bed_info,1)
tail(vcf_bed_info)

vcf_bed_info_regular<-vcf_bed_info[vcf_bed_info$SV_CLASS=="regular",]

write.table(vcf_bed_info[,c(1,2,15,9,12)],ANNOT_OUTPUT,sep="\t", row.names=F, quote=F, col.names=F)

write.table(vcf_bed_info_regular[,c(1,2,15,9,12)],paste0(ANNOT_OUTPUT,"regular_only.txt"),sep="\t", row.names=F, quote=F, col.names=F)

