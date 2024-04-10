#filter the OpenTargets query to high-quality GWAS
rule run_filter_opentargets_query:
    input:
        input_file='../../resources/OpenTargets_L2G_noQC.csv.gz'
    output: '../../resources/OpenTargets_L2G_Filtered.csv.gz'
    log: '../../logs/4_1_run_filter_OpenTargets_Query.log'
    params:
        min_l2g_score=0.5,
        remove_mhc_region=True
    script:
        "../scripts/4_1_filter_open_targets_to_high_quality_gwas.py"