import os

# Working directory from config
workdir: config['workdir']

# One rule to ring them all
rule all:
    input: 'resources/hocomoco_meme.meme'

# Load data in .h5mu format
# https://mudata.readthedocs.io/en/latest/
rule load_data:
    output: 'evaluation_mdata.h5mu'
    log: 'logs/0_load_data.log'
    script: 'scripts/0_load_data.py'

# Run technical evaluations
rule run_technical_evaluations:
    input: 'evaluation_mdata.h5mu'
    output: touch('touch/1_run_technical_evaluations.done')
    log: 'logs/1_run_technical_evaluations.log'
    script: 'scripts/1_run_technical_evaluations.py'

# Run gene-set enrichment

# Run motif enrichment
include: 'rules/3_run_motif_enrichment.smk'




    

