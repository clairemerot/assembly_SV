#!/bin/bash
#SBATCH -J "nucmerAMER"
#SBATCH -o log_%j
#SBATCH -c 10 
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=21-00:00
#SBATCH --mem=200G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR


#align one genome to the other
NB_CPU=10

INPUT_FOLDER=../genome_masked
OUTPUT_FOLDER=04_nucmer_symap

REFERENCE=normal.masked.chrsonly.fasta
QUERY=normal.masked.chrsonly.fasta


echo "aligning" "${QUERY%.fasta}" "to" "${REFERENCE%.fasta}"
#using the parameters from Weissensteiner et al
# -l is the minimum length of the single match (default =20)
#-c is the minimum length of a cluster of match (defaut 65)
#nucmer --maxmatch -l 100 -c 500 -t $NB_CPU \
#   "$INPUT_FOLDER"/"$REFERENCE" "$INPUT_FOLDER"/"$QUERY" \
#   -p "$OUTPUT_FOLDER"/"${QUERY%.fasta}"vs"${REFERENCE%.fasta}"
#
#show-coords -dlTH "$OUTPUT_FOLDER"/"${QUERY%.fasta}"vs"${REFERENCE%.fasta}".delta > "$OUTPUT_FOLDER"/"${QUERY%.fasta}"vs"${REFERENCE%.fasta}".mum

nucmer -t $NB_CPU \
   "$INPUT_FOLDER"/"$REFERENCE" "$INPUT_FOLDER"/"$QUERY" \
   -p "$OUTPUT_FOLDER"/"${QUERY%.fasta}"vs"${REFERENCE%.fasta}".all

show-coords -dlTH "$OUTPUT_FOLDER"/"${QUERY%.fasta}"vs"${REFERENCE%.fasta}".all.delta > "$OUTPUT_FOLDER"/"${QUERY%.fasta}"vs"${REFERENCE%.fasta}".all.mum


#Self synteny:
#Use nucmer. For the same chromosome against itself, use the parameter "--maxmatch".
#The output of nucmer is input to "show-coords -dlTH".
#Result files:
#The result files must have suffix ".mum".
#Put the ".mum" files in the directory data/seq_results/<proj1-to-proj2>/align.
#In the <proj1-to-proj2>/align directory, execute:
#touch all.done
#This creates a file, which indicates to SyMAP not to do any alignment, but to process the files in the directory ending with ".mum".
#When you run "Selected Pair", SyMAP will recognize the files and use them to build the synteny blocks.