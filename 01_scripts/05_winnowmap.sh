#!/bin/bash
#SBATCH -J "winnomap"
#SBATCH -o log_%j
#SBATCH -c 1 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=7-00:00
#SBATCH --mem=100G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR


#align one genome to the other
#variables 
NB_CPU=1
INPUT_FOLDER=02_genomes
OUTPUT_FOLDER=05_winnowmap

TARGET=normal_chrs.fasta
QUERY=dwarf.fasta

KMER_SIZE=15


echo "aligning" "${TARGET%.fasta}" "to" "${QUERY%.fasta}"

meryl count k=$KMER_SIZE output "$OUTPUT_FOLDER"/merylDB "$INPUT_FOLDER"/"$TARGET"
meryl print greater-than distinct=0.9998 05_winnowmap/merylDB > "$OUTPUT_FOLDER"/repetitive_k"$KMER_SIZE".txt

winnowmap -k $KMER_SIZE -W "$OUTPUT_FOLDER"/repetitive_k"$KMER_SIZE".txt -ax asm5 "$INPUT_FOLDER"/"$TARGET" "$INPUT_FOLDER"/"$QUERY" > "$OUTPUT_FOLDER"/"${TARGET%.fasta}"VS"${QUERY%.fasta}".paf

#use asm5 for very similar genoem <5% divergence, and asm20 for very divergent (20%)

