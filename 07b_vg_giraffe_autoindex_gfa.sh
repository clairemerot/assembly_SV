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


gfa_graph=07_minigraph/normal_chrsVSdwarf.gfa
ref=02_genomes/normal_chrs.fasta

mkdir 07_minigraph/TMPDIR
#use the autoindex workflow to build the graph into vg format
vg autoindex --workflow giraffe \
--prefix 07_minigraph/normal_chrsVSdwarf_giraffe \
--tmp-dir 07_minigraph/TMPDIR \
--target-mem 190G --threads $NB_CPU \
-R XG \
--ref-fasta $ref \
--gfa $gfa_graph


