#!/bin/bash
#SBATCH --job-name=isoquant
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=16
#SBATCH --mem=800G
#SBATCH --time=150:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

# Activate environment
source /gscratch/stergachislab/asedeno/tools/miniconda3/bin/activate
conda activate isoquant

isoquant.py --reference /gscratch/stergachislab/assemblies/simple-names/hg38.fa \
  --genedb /mmfs1/gscratch/stergachislab/asedeno/data/reference_files/gencode.v47.annotation.gtf \
  --bam /gscratch/stergachislab/asedeno/data/LRAA/all.sorted.bam \
  --data_type pacbio --threads 16 --high_memory --complete_genedb --output isoquant_out2 --sqanti_output



