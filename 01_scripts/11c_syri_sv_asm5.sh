#!/bin/bash
#SBATCH -J "syri"
#SBATCH -o log_%j
#SBATCH -c 8 
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=claire.merot@gmail.com
#SBATCH --time=1-00:00
#SBATCH --mem=50G

# Important: Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

module load samtools

#align one genome to the other with minimap2
#variables 
NB_CPU=8
REF=02_genomes/normal_chrs.fasta
QUERY=11_syri/dwarf_ragtag_chr.fasta


#minimap2 -ax asm5 -t $NB_CPU --eqx $REF $QUERY > 11_syri/asm5/normalVSdwarf_pseudochrs_asm5.sam
#samtools view -b 11_syri/normalVSdwarf_pseudochrs.sam > 11_syri/normalVSdwarf_pseudochrs.bam


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
PATH_syri=/home/camer78/Softwares/syri/syri/bin/

#run zyri
python3 "$PATH_syri"/syri -c 11_syri/asm5/normalVSdwarf_pseudochrs_asm5.sam -r $REF -q $QUERY -k -F S --dir 11_syri/asm5 --prefix normalVSdwarf_pseudochrs_asm5_essai --nc $NB_CPU

#plot
python $PATH_syri/plotsr 11_syri/asm5/normalVSdwarf_pseudochrs_asm5_essaisyri.out 02_genomes/normal_chrs.fasta 11_syri/dwarf_ragtag_chr.fasta -o pdf
