#!/bin/bash
#SBATCH -J "SVMU"
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
#NB_CPU=10

INPUT_FOLDER=02_genomes
OUTPUT_FOLDER=10_svmu

REFERENCE=normal_chrs.fasta
QUERY=dwarf.fasta

echo "running svmu on" "${QUERY%.fasta}" "to" "${REFERENCE%.fasta}"

#the lastz file is empty
svmu "$OUTPUT_FOLDER"/sam2ref.delta "$INPUT_FOLDER"/"$REFERENCE" "$INPUT_FOLDER"/"$QUERY" l "$OUTPUT_FOLDER"/sam_lastz.txt "$OUTPUT_FOLDER"/sam2ref_svmu 
