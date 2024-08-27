#!/bin/bash

# Params
n_components=12
path_script=/cellar/users/aklie/opt/gene_program_evaluation/src/inference/program_models/factor_analysis/factor_analysis.py
path_data=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/analysis/condition/factor_analysis/fa.h5mu
path_config=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/bin/data_analysis/factor_analysis/configs/K${n_components}.gin
path_out=/cellar/users/aklie/data/datasets/igvf_sc-islet_10X-Multiome/analysis/condition/factor_analysis/fa_K${n_components}.h5mu

# Script
cmd="python $path_script \
$path_data \
--config_path $path_config \
--prog_key factor_analysis_K${n_components} \
--data_key rna \
--layer log1p_norm \
--output \
--path_out $path_out"
echo $cmd
eval $cmd
