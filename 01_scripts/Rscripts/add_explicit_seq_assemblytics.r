argv <- commandArgs(T)
INPUT<- argv[1]
OUTPUT<- argv[2]
GENOME<- argv[3]

source("01_scripts/Rscripts/fix_sniffles_assemblytics.R")

fix_sniffles(input_vcf=INPUT, output_vcf=OUTPUT, refgenome = GENOME)  
