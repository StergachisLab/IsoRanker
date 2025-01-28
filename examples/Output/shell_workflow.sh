#!/bin/bash
# Hank Cheng
# 2025/01

# Usage bash shell_workflow.sh > shell_workflow.output.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate

conda activate isoranker_testing

pip install --force-reinstall git+https://github.com/yhhc2/IsoRanker.git

run_analysis \
  --read_stat_path /mmfs1/gscratch/stergachislab/yhhc/projects/IsoRanker/examples/Input/read_stats.txt \
  --sample_info_path /mmfs1/gscratch/stergachislab/yhhc/projects/IsoRanker/examples/Input/Sample_info.csv \
  --classification_path /mmfs1/gscratch/stergachislab/yhhc/projects/IsoRanker/examples/Input/filtered_classification.txt \
  --genemap_path /mmfs1/gscratch/stergachislab/yhhc/projects/IsoRanker_testing/genemap2.txt \
  --output_dir ./Output