#!/bin/bash
#SBATCH -J "07_minigraph"
#SBATCH -o log_%j
#SBATCH -c 8 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=7-00:00
#SBATCH --mem=100G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR



#align one genome to the other
#variables 
NB_CPU=8
INPUT_FOLDER=02_genomes
OUTPUT_FOLDER=07_minigraph

TARGET=normal_chrs.fasta
QUERY=dwarf.fasta


echo "generate a graph with" "${TARGET%.fasta}" "to" "${QUERY%.fasta}"
minigraph -xggs -t8 "$INPUT_FOLDER"/"$TARGET" "$INPUT_FOLDER"/"$QUERY" > "$OUTPUT_FOLDER"/"${TARGET%.fasta}"VS"${QUERY%.fasta}".gfa


