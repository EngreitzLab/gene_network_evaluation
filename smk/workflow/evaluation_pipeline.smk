import os

# Working directory from config
workdir: config['workdir']

# One rule to ring them all
rule all:
    input: 
        'touch/1_run_technical_evaluations.done',
        'touch/2_run_geneset_enrichment.done',
        'touch/3_1_run_enhancer_motif_enrichment.done',
        'touch/3_2_run_promoter_motif_enrichment.done',
        'resources/OpenTargets_L2G_Filtered.csv.gz'
        
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
rule run_geneset_enrichment:
    input: 'evaluation_mdata.h5mu'
    output: touch('touch/2_run_geneset_enrichment.done')
    log: 'logs/2_run_geneset_enrichment.log'
    script: 'scripts/2_run_geneset_enrichment.py'

# Run motif enrichment
include: 'rules/3_run_motif_enrichment.smk'

# Run open target analysis
include: 'rules/4_run_query_open_targets_analysis.smk'




    

