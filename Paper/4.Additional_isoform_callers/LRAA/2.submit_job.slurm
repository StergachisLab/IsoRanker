#!/bin/bash
#SBATCH --job-name=lraa_test
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=20
#SBATCH --mem=700G
#SBATCH --time=240:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

# Activate environment
source /gscratch/stergachislab/asedeno/tools/miniconda3/bin/activate
conda activate lraa_env
export MINIWDL__RUNTIME_DEFAULTS__CONTAINER_ENGINE=singularity
export PYTHONIOENCODING=latin1

# Set folder
INPUT_DIR="$PWD/lraa_inputs"
INPUT_FILES=($INPUT_DIR/*.json)

# If submitting manually: detect number of files
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Please submit this script with --array=0-$((${#INPUT_FILES[@]}-1))"
    exit 1
fi

# Pick the correct JSON
JSON_FILE="${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}"

echo "Running LRAA for input JSON: $JSON_FILE"

miniwdl run /gscratch/stergachislab/asedeno/tools/LongReadAlignmentAssembler/WDL/LRAA.wdl \
  --input "$JSON_FILE" \
  --verbose
