#!/bin/bash
#SBATCH -J "svim-asm"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=1-00:00
#SBATCH --mem=100G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR


NB_CPU=1
module load bcftools/1.13

#temporary folder for intermediate files
VCF_FOLDER=13_svimAsm/TMP_VCF
INPUT_VCF=13_svimAsm/variants.vcf # we will try not to modify this one - it was generated by calling snarls from the gfa graph with vg
OUTPUT_VCF=13_svimAsm/svimAsm_filtered.vcf

mkdir $VCF_FOLDER
cp $INPUT_VCF $VCF_FOLDER/raw.vcf.gz
gunzip $VCF_FOLDER/raw.vcf.gz

#check what it looks like
#grep -v ^\#\# $VCF_FOLDER/raw.vcf | head 
#tail $VCF_FOLDER/raw.vcf 
#head $VCF_FOLDER/raw.vcf 
echo "total number of SVs"
grep -v ^\#\# $VCF_FOLDER/raw.vcf | wc -l  #172200


#filter out BND
#DUP¨are considered as insertions
bcftools filter -i'INFO/SVTYPE!="BND"' -o $VCF_FOLDER/raw.noBND.vcf -Ov $VCF_FOLDER/raw.vcf 
echo "total number of SVs restricted to INS, DEL, INV"
grep -v ^\#\# $VCF_FOLDER/raw.noBND.vcf | wc -l #170582

#Export sequences for advanced filtering
bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' $VCF_FOLDER/raw.noBND.vcf > $VCF_FOLDER/SV_data_with_seq.txt

#blacklist because of N string > 10 (possible junction of contigs 
grep -P "N{10,}" $VCF_FOLDER/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' > $VCF_FOLDER/N10_blacklist.bed
echo "SVs excluded because of >10N" 
wc -l $VCF_FOLDER/N10_blacklist.bed

bgzip -c $VCF_FOLDER/N10_blacklist.bed > $VCF_FOLDER/N10_blacklist.bed.gz
tabix -s1 -b2 -e2 $VCF_FOLDER/N10_blacklist.bed.gz

#remove blacklist of variants
bcftools view -T ^$VCF_FOLDER/N10_blacklist.bed.gz $VCF_FOLDER/raw.noBND.vcf > $VCF_FOLDER/raw.noBND_Nfiltered.vcf
grep -v ^\#\# $VCF_FOLDER/raw.noBND_Nfiltered.vcf | wc -l #170568

#keep the filtered vcf
cp $VCF_FOLDER/raw.noBND_Nfiltered.vcf $OUTPUT_VCF
bgzip -c $OUTPUT_VCF > $OUTPUT_VCF.gz
tabix $OUTPUT_VCF.gz