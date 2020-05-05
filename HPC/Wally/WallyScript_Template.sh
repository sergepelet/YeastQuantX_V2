#!/bin/bash

#SBATCH --account=spelet1_single_cells
#SBATCH --array=1-NUM_FILES
#SBATCH --job-name=JOB_ID
#SBATCH --partition wally

#SBATCH --chdir OUT_DIR/
#SBATCH --output=%x_%A-%4a.out
#SBATCH --error=%x_%A-%4a.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=60
#SBATCH --export=NONE

# Add Matlab as Module
module add Development/Languages/Matlab_Compiler_Runtime/v96

# use sed to extract the content of the file list one by one based on Task_ID
VARLIST=VAR_FILE_PATH
VAR=$(sed -n ${SLURM_ARRAY_TASK_ID}p $VARLIST)
echo "Running analysis on" $VAR

# Perform analysis
/scratch/wally/FAC/FBM/DMF/spelet1/single_cells/D2c/YQ_VX/src/YQ_VX $VAR