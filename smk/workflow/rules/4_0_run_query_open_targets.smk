#run full query on OpenTargets
rule run_opentargets_query:
    input:
        credentials_path='/home/robertg1/.ssh/test-bigquery-ot-956f8a01208f.json'
    output: '../../resources/OpenTargets_L2G_noQC.csv.gz'
    log: '../../logs/4_0_run_full_query_on_OpenTargets.log'
    params:
        min_assoc_loci=10,
        min_n_cases=1000,
        min_l2g_score=0.2,
        study_ids_to_keep=None
    script:
        "../scripts/4_0_query_open_targets_genetics.py"
