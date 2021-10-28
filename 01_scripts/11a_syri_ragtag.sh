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


query=dwarf
ref=normal_chrs

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

 module load minimap2/2.17



#make a pseudo chromosome dwarf based on normal order

ragtag.py scaffold 02_genomes/normal_chrs.fasta 02_genomes/dwarf.fasta -t 10 -o 11_syri
