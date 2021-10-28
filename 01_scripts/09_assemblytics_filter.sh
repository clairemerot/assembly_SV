#!/bin/bash
#SBATCH -J "assemblytics_filter"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=1-00:00
#SBATCH --mem=10G
cd $SLURM_SUBMIT_DIR

###this script will format and filter the vcf from assmeblytics
module load htslib/1.10.2
module load bcftools/1.13
module load vcftools

#temporary folder for intermediate files
VCF_FOLDER=09_assemblytics/TMP_VCF
INPUT_VCF=09_assemblytics/TMP_VCF/SVassemblytics.withseq.vcf.gz # we will try not to modify this one
OUTPUT_VCF=09_assemblytics/assemblytics_filtered.vcf

#mkdir $VCF_FOLDER
cp $INPUT_VCF $VCF_FOLDER/raw.vcf.gz
gunzip $VCF_FOLDER/raw.vcf.gz

#check what it oolk like
grep -v ^\#\# $VCF_FOLDER/raw.vcf | head 
tail $VCF_FOLDER/raw.vcf 
echo "total number of SVs"
grep -v ^\#\# $VCF_FOLDER/raw.vcf | wc -l  #97654


###some filtering

#filter vcf -i (include, -O vcf format -o
bcftools filter -i'INFO/SVLEN<=100000 && INFO/SVLEN>=-100000' -o $VCF_FOLDER/raw_sorted.noTRA_100k.vcf -Ov $VCF_FOLDER/raw.vcf
echo "total number of SVs < 100kb"
grep -v ^\#\# $VCF_FOLDER/raw_sorted.noTRA_100k.vcf | wc -l #96658



#Export sequences for advanced filtering
bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' $VCF_FOLDER/raw_sorted.noTRA_100k.vcf > $VCF_FOLDER/SV_data_with_seq.txt

#blacklist because of N string > 10 (possible junction of contigs 
grep -P "N{10,}" $VCF_FOLDER/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' > $VCF_FOLDER/N10_blacklist.bed
echo "SVs excluded because of >10N" 
wc -l $VCF_FOLDER/N10_blacklist.bed
head $VCF_FOLDER/N10_blacklist.bed

bgzip -c $VCF_FOLDER/N10_blacklist.bed > $VCF_FOLDER/N10_blacklist.bed.gz
tabix -s1 -b2 -e2 $VCF_FOLDER/N10_blacklist.bed.gz

#remove blacklist of variants
bcftools view -T ^$VCF_FOLDER/N10_blacklist.bed.gz $VCF_FOLDER/raw_sorted.noTRA_100k.vcf > $VCF_FOLDER/raw_sorted.noTRA_100k_Nfiltered.vcf
echo "SVs after filtration for N seq" 
grep -v ^\#\# $VCF_FOLDER/raw_sorted.noTRA_100k_Nfiltered.vcf | wc -l #96555

#keep the filtered vcf
(grep ^"#" $VCF_FOLDER/raw_sorted.noTRA_100k_Nfiltered.vcf; grep -v ^"#" $VCF_FOLDER/raw_sorted.noTRA_100k_Nfiltered.vcf | sort -k1,1 -k2,2n) > $VCF_FOLDER/raw_sorted.noTRA_100k_Nfiltered_sorted.vcf
cp $VCF_FOLDER/raw_sorted.noTRA_100k_Nfiltered_sorted.vcf $OUTPUT_VCF
bgzip -c $OUTPUT_VCF > $OUTPUT_VCF.gz
tabix $OUTPUT_VCF.gz

#clean intermediate files
#rm -r $VCF_FOLDER