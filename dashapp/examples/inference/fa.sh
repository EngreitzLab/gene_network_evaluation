#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/bin/slurm_logs/%x.%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01-00:00:00
#SBATCH --array=1-10%10

#####
# USAGE:
# sbatch --job-name=factor_analysis fa.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Params
n_components=(5 10 15 20 25 30 35 40 45 50)
n_components=${n_components[$SLURM_ARRAY_TASK_ID-1]}
path_script=/cellar/users/aklie/opt/gene_program_evaluation/src/inference/program_models/factor_analysis/factor_analysis.py
path_data=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/analysis/condition/factor_analysis/fa.h5mu
path_config=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/bin/data_analysis/factor_analysis/configs/K${n_components}.gin
path_out=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/analysis/condition/factor_analysis/fa_K${n_components}.h5mu

# Echo inputs and parameters
echo -e "n_components: $n_components"
echo -e "path_script: $path_script"
echo -e "path_data: $path_data"
echo -e "path_config: $path_config"
echo -e "path_out: $path_out\n"

# Script
cmd="python $path_script \
$path_data \
--config_path $path_config \
--prog_key factor_analysis_K${n_components} \
--data_key rna \
--layer log1p_norm \
--output \
--path_out $path_out"
echo "Running command: $cmd"
eval $cmd

date
echo "Job completed."
