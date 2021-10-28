#!/bin/bash
#SBATCH -J "ragtag"
#SBATCH -o log_%j
#SBATCH -c 10 
#SBATCH -p large
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=21-00:00
#SBATCH --mem=100G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
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
#<<< conda initialize <<<


#for Syri to work we need to rename teh chromosomes to be exactly the same
#variables
INPUT_GENOME=11_syri/ragtag.scaffold.fasta
OUTPUT_GENOME=11_syri/dwarf_ragtag_chr.fasta
MIN_SIZE_LG=10000000 #below this threshold, contigs are removed
CORR_LG=11_syri/ragtag.scaffold.rename #tab-delim file with two columns (current name, wanted name)

#rename Lg above a given size with a file of correpsondence
#uses Eric python scripts
# <program> input_fasta correspondence min_length output_fasta
python 01_scripts/utility/rename_scaffolds.py $INPUT_GENOME $CORR_LG $MIN_SIZE_LG $OUTPUT_GENOME
