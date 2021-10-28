argv <- commandArgs(T)
BED_FILE<- argv[1]





library(tidyr)
#read input bed from Assemblytics
bed<-read.table(BED_FILE, header=T, stringsAsFactor=F)
head(bed)
#make id uniq
bed$ID_CHR<-paste(bed$ID, bed$reference, sep="_")
#split query coordinates
bed.split<-separate (bed, query_coordinates,c("query","qry_pos","qry_strand"),sep=":")
head(bed.split)
bed.split2<-separate (bed.split, qry_pos,c("qry_start","qry_stop"),sep="-")
head(bed.split2)

#keep info for reference sequence coordinates
ref_info<- bed.split2[,c(15,1,2,3,6)]
head(ref_info)

#keep info for alternative sequence coordinates
qry_info<- bed.split2[,c(15,10,11,12,13)]
head(qry_info)

#format with uniq id
bed_formatted<-cbind(bed.split2[,c(1,2,3,15)], bed[,5:11])
colnames(bed_formatted)<-colnames(bed)[1:11]
head(bed_formatted)

#save output
write.table(ref_info,paste0(BED_FILE,"_ref_info.tsv"),sep="\t", row.names=F, col.names=F, quote=F)
write.table(qry_info,paste0(BED_FILE,"_qry_info.tsv"),sep="\t", row.names=F, col.names=F, quote=F)
write.table(bed_formatted,paste0(BED_FILE,"_formatted.bed"),sep="\t", row.names=F, quote=F)