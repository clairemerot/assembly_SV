#!/bin/bash
#SBATCH -J "08_unimap"
#SBATCH -o log_%j
#SBATCH -c 8 
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=1-00:00
#SBATCH --mem=140G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR


#align one genome to the other
#variables 
NB_CPU=8
INPUT_FOLDER=02_genomes
OUTPUT_FOLDER=08_unimap

TARGET=normal_chrs.fasta
QUERY=dwarf.fasta


echo "aligning" "${TARGET%.fasta}" "to" "${QUERY%.fasta}"
unimap -b24 -d "$INPUT_FOLDER"/"$TARGET".umi "$INPUT_FOLDER"/"$TARGET"
unimap -cxasm5 -t$NB_CPU "$INPUT_FOLDER"/"$TARGET".umi "$INPUT_FOLDER"/"$QUERY" > "$OUTPUT_FOLDER"/"${TARGET%.fasta}"VS"${QUERY%.fasta}".asm5.paf



