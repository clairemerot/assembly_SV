#!/bin/bash
#SBATCH -J "07_bubble"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=1-00:00
#SBATCH --mem=10G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR



#align one genome to the other
#variables 
NB_CPU=1
INPUT_FOLDER=02_genomes
OUTPUT_FOLDER=07_minigraph

TARGET=normal_chrs.fasta
QUERY=dwarf.fasta


echo "output bubble as SV" "${TARGET%.fasta}" "to" "${QUERY%.fasta}"

gfatools bubble "$OUTPUT_FOLDER"/"${TARGET%.fasta}"VS"${QUERY%.fasta}".gfa > "$OUTPUT_FOLDER"/"${TARGET%.fasta}"VS"${QUERY%.fasta}".SV.bed