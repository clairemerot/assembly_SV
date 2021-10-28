#!/bin/bash
#SBATCH -J "svim-asm"
#SBATCH -o log_%j
#SBATCH -c 5 
#SBATCH -p medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=7-00:00
#SBATCH --mem=100G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR


#align one genome to the other
#variables 
NB_CPU=5
INPUT_FOLDER=02_genomes

REF=02_genomes/normal_chrs.fasta
ASSEMBLY=02_genomes/dwarf.fasta

module load samtools


echo "aligning $REF and $ASSEMBLY"
minimap2 -a -x asm5 --cs -r2k -t $NB_CPU $REF $ASSEMBLY > 13_svimAsm/normal_chrsVSdwarf.sam
echo "sort and index"
samtools sort -m4G -@4 -o 13_svimAsm/normal_chrsVSdwarf.sorted.bam 13_svimAsm/normal_chrsVSdwarf.sam
samtools index 13_svimAsm/normal_chrsVSdwarf.sorted.bam


#load the program. It requires python 3
__conda_setup="$('/home/camer78/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/camer78/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/camer78/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/camer78/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate svimasm_env


echo "run svim-asm"
svim-asm haploid --min_sv_size 50 --max_sv_size 200000 --tandem_duplications_as_insertions --interspersed_duplications_as_insertions 13_svimAsm 13_svimAsm/normal_chrsVSdwarf.sorted.bam $REF 