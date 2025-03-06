#!/bin/bash
# Hank Cheng
# 2025/01

# Usage bash shell_workflow.sh > shell_workflow.output.txt 2>&1

source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate

conda activate testing_isoranker

pip install git+https://github.com/yhhc2/IsoRanker.git

isoranker_run_analysis \
  --read_stat_path ./IsoRanker/examples/Input/read_stats.txt \
  --sample_info_path ./IsoRanker/examples/Input/Sample_info.csv \
  --classification_path ./IsoRanker/examples/Input/filtered_classification.txt \
  --genemap_path ./IsoRanker/examples/Input/genemap2_placeholder.txt \
  --output_dir ./Output
