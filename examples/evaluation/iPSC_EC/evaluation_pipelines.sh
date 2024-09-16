#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --job-name=evaluation_pipelines
#SBATCH --output=logs/evaluation_pipelines_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-26

#####
# USAGE:
# sbatch evaluation_pipelines.sh <config>
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID"

# Configuring env
source activate /cellar/users/aklie/opt/miniconda3/envs/test_celloracle

# Input configs should be one per line with tab indentation
configs=(
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF/cNMF_200/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF/cNMF_80/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF/cNMF_30/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF/cNMF_300/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF/cNMF_60/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF/cNMF_100/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/cNMF/cNMF_250/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_13/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_40/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_17/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_50/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_12/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_16/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_45/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_18/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_35/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_10/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_14/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_25/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_20/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_11/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_15/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_30/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_19/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/Topyfic/Topyfic_5/evaluation_pipeline.yml
    /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/iPSC_EC/factor_analysis/factor_analysis_100/evaluation_pipeline.yml
)
config=${configs[$SLURM_ARRAY_TASK_ID-1]}

# Run the command
echo "Running command: python /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/evaluation_pipeline.py --config $config"
python /cellar/users/aklie/opt/gene_program_evaluation/examples/evaluation/evaluation_pipeline.py --config $config

# Date
date
