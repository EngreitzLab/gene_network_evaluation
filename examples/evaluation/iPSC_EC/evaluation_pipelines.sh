#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --job-name=evaluation_pipelines
#SBATCH --output=logs/evaluation_pipelines_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-8%8

#####
# USAGE:
# sbatch evaluation_pipelines.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID"

# Configuring env
source activate /cellar/users/aklie/opt/miniconda3/envs/test_celloracle

# Input configs should be one per line with tab indentation
configs=(
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF_250/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF_80/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF_30/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF_100/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF_300/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF_60/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF_200/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/factor_analysis_100/evaluation_pipeline.yml
)
config=${configs[$SLURM_ARRAY_TASK_ID-1]}

# Run the command
echo "Running command: python /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/evaluation_pipeline.py --config $config"
python /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/evaluation_pipeline.py --config $config

# Date
date
