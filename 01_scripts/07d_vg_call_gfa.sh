#!/bin/bash
#SBATCH -J "vg_autoindex_gfa"
#SBATCH -o vg_autoindex_%j.out
#SBATCH -c 8 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=7-00:00
#SBATCH --mem=200G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Loading the htslib module for compressing and indexing
#module load htslib/1.10.2
NB_CPU=8


ref=02_genomes/normal_chrs.fasta
graph=07_minigraph/normal_chrsVSdwarf_giraffe.xg
snarls=07_minigraph/normal_chrsVSdwarf_giraffe.pb

#get snarls
vg snarls -t $NB_CPU $graph > $snarls

#we will map just one sample to be able to call the vcf
id=ID4
ind=$(grep $id ../genotyping_SV/02_info/fq.list)

echo $id
echo $ind

# Creating variables for the location of the input files
fq1=../wgs_sample_preparation_ALL/05_trimmed/${ind}_1.trimmed.fastq.gz
fq2=../wgs_sample_preparation_ALL/05_trimmed/${ind}_2.trimmed.fastq.gz


# Mapping the paired reads
#this requires the graph.xg and the graph.gcsa (-d is their BASENAME)
# -t threads -f for the fastq

vg giraffe -t $NB_CPU \
-H 07_minigraph/normal_chrsVSdwarf_giraffe.giraffe.gbwt \
-g 07_minigraph/normal_chrsVSdwarf_giraffe.gg \
-m 07_minigraph/normal_chrsVSdwarf_giraffe.min \
-d 07_minigraph/normal_chrsVSdwarf_giraffe.dist \
-f $fq1 -f $fq2 -N $id > 07_minigraph/"$id"_paired.gam
#pack
vg pack -t $NB_CPU -Q 5 -x $graph -g 07_minigraph/"$id"_paired.gam -o 07_minigraph/"$id".pack

#call
vg call -t $NB_CPU -a -k 07_minigraph/"$id".pack -r $snarls -f $ref $graph > 07_minigraph/normal_chrsVSdwarf.vcf

