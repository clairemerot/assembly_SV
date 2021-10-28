#!/bin/bash
#SBATCH -J "syri_filter"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=1-00:00
#SBATCH --mem=1G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

module load bcftools

BASE=11_syri/asm5/normalVSdwarf_pseudochrs_asm5_essaisyri
OUTPUT_VCF=11_syri/syri_asm5_filtered.vcf

grep -v ^\#\# "$BASE".vcf | wc -l #5071181

## step 0 separate vcf into vcf of small SVs and SR (remov ealigned regions, translocatins,e tc
#remove alignment region
bcftools filter -i 'ALT!="<NOTAL>" & ALT!="<SYNAL>" & ALT!="<INVAL>" & ALT!="<INVDPAL>" & ALT!="<INVTRAL>" & ALT!="<DUPAL>" & ALT!="<TRANSAL>"' -o "$BASE"_noAL.vcf -Ov "$BASE".vcf
grep -v ^\#\# "$BASE"_noAL.vcf  | wc -l # 5022324

#keep only SV & SNPs
bcftools filter -i 'ALT!="<HDR>" & ALT!="<TDM>" & ALT!="<TRANS>" & ALT!="<INVTR>" & ALT!="<SYN>"' -o "$BASE"_noAL_SV.vcf -Ov "$BASE"_noAL.vcf
grep -v ^\#\# "$BASE"_noAL_SV.vcf  | wc -l # 5022324
less "$BASE"_noAL_SV.vcf

#deal with SR, CPG, CPL
bcftools filter -i 'INFO/VarType=="SR" | ALT=="<CPG>" | ALT=="<CPL>"' -o "$BASE"_noAL_SV_SR.vcf -Ov "$BASE"_noAL_SV.vcf
grep -v ^\#\# "$BASE"_noAL_SV_SR.vcf  | wc -l # 5282
tail "$BASE"_noAL_SV_SR.vcf


#deal with small SVs
bcftools filter -i 'INFO/VarType=="ShV" & ALT!="<CPG>" & ALT!="<CPL>"' -o "$BASE"_noAL_SVsmall_SNP.vcf -Ov "$BASE"_noAL_SV.vcf
grep -v ^\#\# "$BASE"_noAL_SVsmall_SNP.vcf  | wc -l # 4979114
tail "$BASE"_noAL_SVsmall_SNP.vcf


## step 1 annotate vcf of small SVs

#Export sequences to get length & type
bcftools query -f '%CHROM %POS %ID %INFO/END %REF %ALT\n' "$BASE"_noAL_SVsmall_SNP.vcf > "$BASE"_noAL_SVsmall_SNP.info
#use R to extrct SVTYPE and SVLEN
Rscript 01_scripts/Rscripts/extract_SV_info_syri.r "$BASE"_noAL_SVsmall_SNP.info

#add info on vcf about type and length
bgzip "$BASE"_noAL_SVsmall_SNP.info.annot
tabix -s1 -b2 -e2 "$BASE"_noAL_SVsmall_SNP.info.annot.gz

#prepare the header
echo -e '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' > "$BASE"_noAL_SVsmall_SNP.info.annot.hdr
echo -e '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' >> "$BASE"_noAL_SVsmall_SNP.info.annot.hdr

#-a is the annotation file (tabix and bgzip, it needs at least CHROM and POS, -h are the header lines to add, -c are the meaning of the column in the annotation file
bgzip "$BASE"_noAL_SVsmall_SNP.vcf
tabix "$BASE"_noAL_SVsmall_SNP.vcf.gz
bcftools annotate -a "$BASE"_noAL_SVsmall_SNP.info.annot.gz -h "$BASE"_noAL_SVsmall_SNP.info.annot.hdr -c CHROM,POS,INFO/SVLEN,INFO/SVTYPE "$BASE"_noAL_SVsmall_SNP.vcf.gz > "$BASE"_noAL_SVsmall_SNP.annotated.vcf
grep -v ^\#\# "$BASE"_noAL_SVsmall_SNP.annotated.vcf | head  
grep -v ^\#\# "$BASE"_noAL_SVsmall_SNP.annotated.vcf | wc -l #4979114

#filter out above 50bp
bcftools filter -i'INFO/SVLEN<=-50 | INFO/SVLEN>=50' -o "$BASE"_noAL_SVsmall_SNP.annotated_50.vcf -Ov "$BASE"_noAL_SVsmall_SNP.annotated.vcf
echo "total number of SVs >50bp"
grep -v ^\#\# "$BASE"_noAL_SVsmall_SNP.annotated_50.vcf | wc -l #63902

## step 2 annotate vcf of larger SRs
#Export sequences to get length & type
bcftools query -f '%CHROM %POS %ID %INFO/END %INFO/StartB %INFO/EndB %REF %ALT\n' "$BASE"_noAL_SV_SR.vcf > "$BASE"_noAL_SV_SR.info
#use R to extrct SVTYPE and SVLEN
Rscript 01_scripts/Rscripts/extract_SV_info_syri_large_SR.r "$BASE"_noAL_SV_SR.info

#add info on vcf about type and length
bgzip "$BASE"_noAL_SV_SR.info.annot
tabix -s1 -b2 -e2 "$BASE"_noAL_SV_SR.info.annot.gz

#prepare the header
echo -e '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' > "$BASE"_noAL_SV_SR.info.annot.hdr
echo -e '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' >> "$BASE"_noAL_SV_SR.info.annot.hdr

#-a is the annotation file (tabix and bgzip, it needs at least CHROM and POS, -h are the header lines to add, -c are the meaning of the column in the annotation file
bgzip "$BASE"_noAL_SV_SR.vcf
tabix "$BASE"_noAL_SV_SR.vcf.gz
bcftools annotate -a "$BASE"_noAL_SV_SR.info.annot.gz -h "$BASE"_noAL_SV_SR.info.annot.hdr -c CHROM,POS,INFO/SVLEN,INFO/SVTYPE "$BASE"_noAL_SV_SR.vcf.gz > "$BASE"_noAL_SV_SR.annotated.vcf
grep -v ^\#\# "$BASE"_noAL_SV_SR.annotated.vcf | head  
grep -v ^\#\# "$BASE"_noAL_SV_SR.annotated.vcf | wc -l #5282

#filter out above 50bp
bcftools filter -i'INFO/SVLEN<=-50 | INFO/SVLEN>=50' -o "$BASE"_noAL_SV_SR.annotated_50.vcf -Ov "$BASE"_noAL_SV_SR.annotated.vcf
echo "total number of SVs >50bp"
grep -v ^\#\# "$BASE"_noAL_SV_SR.annotated_50.vcf | wc -l #5244

## step 3 join the two vcfs
#keep the filtered vcf
cat "$BASE"_noAL_SVsmall_SNP.annotated_50.vcf "$BASE"_noAL_SV_SR.annotated_50.vcf > "$BASE"_joined.vcf
grep -v ^\# "$BASE"_joined.vcf | wc -l #69144
(grep ^"#" "$BASE"_noAL_SVsmall_SNP.annotated_50.vcf; grep -v ^"#" "$BASE"_joined.vcf | sort -k1,1 -k2,2n) > "$BASE"_joined.sorted.vcf
cp "$BASE"_joined.sorted.vcf $OUTPUT_VCF
bgzip -c $OUTPUT_VCF > $OUTPUT_VCF.gz
tabix $OUTPUT_VCF.gz
less $OUTPUT_VCF
