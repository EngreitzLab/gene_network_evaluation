#!/bin/bash

#####
# SLURM USAGE:
# sbatch \
# --partition=carter-compute 
# --job-name=factor_analysis \
# --output=%x.%A.out \
# --cpus-per-task=1 \
# --mem=16G \
# --time=01-00:00:00 \
# run_factor_analysis.sh \
# data.h5mu
# config.gin \
# rna \
# log1p_norm \
# path/to/output
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
path_script=path_script="../../../../src/inference/factor_analysis/factor_analysis.py"

# Inputs
path_data=$1
path_config=$2
data_key=$3
layer=$4
path_out=$5
n_components=$(grep -oP 'n_components = \K[0-9]+' $path_config)

# Echo inputs and parameters
echo -e "path_data: $path_data"
echo -e "path_config: $path_config"
echo -e "n_components: $n_components"
echo -e "data_key: $data_key"
echo -e "layer: $layer"
echo -e "path_out: $path_out\n"

# Make output directory
mkdir -p $(dirname $path_out)

# Script
cmd="python $path_script \
$path_data \
--config_path $path_config \
--prog_key factor_analysis_K${n_components} \
--data_key $data_key \
--layer $layer \
--output \
--path_out $path_out"
echo $cmd
eval $cmd

# Date
date
