#specify that snakemake should run all steps necessary to
#produce the filtered OpenTargets query and run enrichment on cell programs
rule all:
  input:
      'resources/OpenTargets_L2G_Cell_Program_Enrichment_Results.csv'

#run full query on OpenTargets
rule run_opentargets_query:
    input:
        config['credentials']
    output: 'resources/OpenTargets_L2G_noQC.csv.gz'
    log: 'logs/4_0_run_full_query_on_OpenTargets.log'
    params:
        min_assoc_loci=10,
        min_n_cases=1000,
        min_l2g_score=0.2,
        study_ids_to_keep=None
    script:
        "../scripts/4_0_query_open_targets_genetics.py"

#filter the OpenTargets query to high-quality GWAS
rule run_filter_opentargets_query:
    input:
        input_file='resources/OpenTargets_L2G_noQC.csv.gz'
    output:
        output_file='resources/OpenTargets_L2G_Filtered.csv.gz'
    log: 'logs/4_1_run_filter_OpenTargets_Query.log'
    params:
        min_l2g_score=0.5,
        remove_mhc_region=True
    script:
        "../scripts/4_1_filter_open_targets_to_high_quality_gwas.py"
        
#using the filtered OpenTargets GWAS file, run GSEA enrichment for each program x GWAS
rule run_gwas_enrichment:
    input:
        gwas_data='resources/OpenTargets_L2G_Filtered.csv.gz'
    output:
        output_file='resources/OpenTargets_L2G_Cell_Program_Enrichment_Results.csv'
    log: 'logs/4_2_run_gwas_enrichment.log'
    params:
        mdata=config['input_loc'],
        prog_key=config['prog_key'],
        data_key=config['data_key'],
        n_jobs=config['n_jobs']
    script:
        "../scripts/4_2_run_gwas_enrichment.py"
