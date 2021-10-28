#!/bin/bash
#SBATCH -J "06_paftools"
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
INPUT_PAF=03_minimap/normal_chrsVSdwarf.sorted.paf
INPUT_PAF=03_minimap/normal_chrsVSdwarf.asm10.sorted.paf
INPUT_PAF=08_unimap/normal_chrsVSdwarf.asm5.paf
REF=02_genomes/normal_chrs.fasta

paftools=/home/camer78/Softwares/minimap2/misc/paftools.js

k8 $paftools call $INPUT_PAF -f $REF -L 50 -l 10 >  06_paftools/SV_paftools
k8 $paftools stat $INPUT_PAF


k8 $paftools sam2paf 13_svimAsm/normal_chrsVSdwarf.sam > 06_paftools/normal_chrsVSdwarf.paf


k8 $paftools call 06_paftools/normal_chrsVSdwarf.paf -f $REF -L 50 -l 10 > 06_paftools/normal_chrsVSdwarf.sv