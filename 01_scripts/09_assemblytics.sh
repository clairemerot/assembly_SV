#!/bin/bash
#SBATCH -J "09_assemblytics"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=1-00:00
#SBATCH --mem=1G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
module load bcftools/1.13

#after running Assemblytics on the web version for each chromosome of the ref (we plit the delta), I decompressed the folder and kept only the bed files
#note: there is a way to do assemblytics with command lines, try for next time 

#all of what is below  is really fast and can be run frontally in <5min

VCF_FOLDER=09_assemblytics/TMP_VCF
mkdir $VCF_FOLDER


###1-gather the bed by chromosome into a single bed
head -n 1 09_assemblytics/bed_by_chr/Whitefish_chr01.Assemblytics_structural_variants.bed > "$VCF_FOLDER"/header.bed
#initialise the bedfile
cp "$VCF_FOLDER"/header.bed "$VCF_FOLDER"/all_chr.bed

ls -1 09_assemblytics/bed_by_chr/*.bed | while read i 
do
echo "adding $i"

tail -n +2 "$i" > "$i".temp
cat 09_assemblytics/all_chr.bed "$i".temp >> "$VCF_FOLDER"/all_chr.bed
done

###2-reformat and extract info
#because we have run assemblytics by chromosome Id are not unique.
# we use a R script to modify the id 
#at the same time we extract the coordinates of the sequence to improve the vcf with actual sequences 

Rscript 01_scripts/Rscripts/extract_SV_info_assemblytics.r "$VCF_FOLDER"/all_chr.bed
##we  use a python script coded by Eric Normandeau to extract the sequence in the fasta
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/camer78/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/camer78/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/camer78/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/camer78/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
#<<< conda initialize <<<
python 01_scripts/utility/fasta_extract_scaffold_regions_sv_claire.py 02_genomes/normal_chrs.fasta  "$VCF_FOLDER"/all_chr.bed_ref_info.tsv "$VCF_FOLDER"/SV_ref_info.fasta
python 01_scripts/utility/fasta_extract_scaffold_regions_sv_claire.py 02_genomes/dwarf.fasta  "$VCF_FOLDER"/all_chr.bed_qry_info.tsv "$VCF_FOLDER"/SV_qry_info.fasta



###3- convert bedfile into vcf
SURVIVOR convertAssemblytics "$VCF_FOLDER"/all_chr.bed_formatted.bed 50 "$VCF_FOLDER"/SVassemblytics.vcf
(grep ^"#" "$VCF_FOLDER"/SVassemblytics.vcf; grep -v ^"#" "$VCF_FOLDER"/SVassemblytics.vcf | sort -k1,1 -k2,2n) > "$VCF_FOLDER"/SVassemblytics.sorted.vcf
rm "$VCF_FOLDER"/SVassemblytics.vcf
#Bed file from Assemblytics
#Min size to keep
#Output vcf file

###4-edit the vcf to add sequences
#we use R to  replace ALT and REF in the vcf
# I wanted to do it with bcftools annotate but it never worked for ALT/REF
#duplication are not so clearly duplication - we chose to label them as INS
grep ^"#" "$VCF_FOLDER"/SVassemblytics.sorted.vcf > "$VCF_FOLDER"/vcf_header
grep -v ^"#" "$VCF_FOLDER"/SVassemblytics.sorted.vcf > "$VCF_FOLDER"/vcf_content
sed -e 's/SVTYPE=DUP/SVTYPE=INS/g' "$VCF_FOLDER"/vcf_content > "$VCF_FOLDER"/vcf_content_edited

Rscript 01_scripts/Rscripts/join_SV_sequences_assemblytics.r "$VCF_FOLDER"/SV_ref_info.fasta "$VCF_FOLDER"/SV_qry_info.fasta "$VCF_FOLDER"/vcf_content_edited "$VCF_FOLDER"/vcf_content_with_seq
cat "$VCF_FOLDER"/vcf_header "$VCF_FOLDER"/vcf_content_with_seq > "$VCF_FOLDER"/SVassemblytics.withseq.vcf
bgzip "$VCF_FOLDER"/SVassemblytics.withseq.vcf
tabix "$VCF_FOLDER"/SVassemblytics.withseq.vcf.gz

#bcftools query -f '%CHROM %POS %ID %INFO/SVTYPE %INFO/SVLEN %ALT %REF\n' "$VCF_FOLDER"/SVassemblytics.withseq.vcf.gz > "$VCF_FOLDER"/SVassemblytics.withseq_info.txt
#head $VCF_FOLDER/raw_info.txt

