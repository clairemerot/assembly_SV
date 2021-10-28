#!/bin/bash
#SBATCH -J "SVMUlastz"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=21-00:00
#SBATCH --mem=100G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR


#align one genome to the other
NB_CPU=1

INPUT_FOLDER=02_genomes
OUTPUT_FOLDER=10_svmu

REFERENCE=normal_chrs.fasta
QUERY=dwarf.fasta

echo "aligning" "${QUERY%.fasta}" "to" "${REFERENCE%.fasta}"
#using the parameters from Weissensteiner et al
# -l is the minimum length of the single match (default =20)
#-c is the minimum length of a cluster of match (defaut 65)
lastz "$INPUT_FOLDER"/"$REFERENCE"[multiple] "$INPUT_FOLDER"/"$QUERY"[multiple] --chain \
   --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > "$OUTPUT_FOLDER"/sam_lastz.txt