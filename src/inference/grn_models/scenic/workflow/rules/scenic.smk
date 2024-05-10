import os

# wildcards from config
outdir: config['outdir']
currdir = os.getcwd()

rule all:
    input:
        expand("{}/run{{run}}_pyscenic_output.loom".format(config['outdir']), run=range(1, config["num_runs"] + 1)),
        "{}/scenic.h5mu".format(config['outdir'])

rule pre:
    input:
        data=config['input_loc']
    output:
        path_out="{outdir}/rna.loom"
    params:
        layer=config['layer']
    log:
        "{outdir}/logs/pre.log"
    benchmark:
        "{outdir}/benchmarks/pre.benchmark"
    shell:
        """
        echo "Preprocessing data"
        python workflow/scripts/pre.py \
        -d {input.data} \
        -l {params.layer} \
        -o {output.path_out} > {log}
        """

rule grn:
    input:
        path_loom="{outdir}/rna.loom",
        path_tf_list=os.path.join(currdir, "resources", "tf_lists", os.path.basename(config['tf_list']))
    output:
        path_adj="{outdir}/run{run}_adj.tsv"
    params:
        method=config.get("inference_method", "grnboost2"),
        run="{run}"
    threads: config['n_jobs']
    log: 
        "{outdir}/logs/run{run}_grn.log"
    benchmark:
        "{outdir}/benchmarks/run{run}_grn.benchmark"
    shell:
        """
        echo "Running GRN inference for run {params.run}"
        pyscenic grn \
        {input.path_loom} \
        {input.path_tf_list} \
        -o {output.path_adj} \
        -m {params.method} \
        --num_workers {threads} > {log}
        """

rule prune:
    input:
        path_adj="{outdir}/run{run}_adj.tsv",
        path_ranking_db=os.path.join(currdir, "resources", "rankings_db", os.path.basename(config['rankings_db'])),
        path_motif_annotations=os.path.join(currdir, "resources", "motif_annotations", os.path.basename(config['motif_annotations'])),
        path_loom="{outdir}/rna.loom"
    output:
        path_reg="{outdir}/run{run}_reg.csv"
    params:
        run="{run}"
    threads: config['n_jobs']
    log:
        "{outdir}/logs/run{run}_prune.log"
    benchmark:
        "{outdir}/benchmarks/run{run}_prune.benchmark"
    shell:
        """
        echo "Running enrichment and pruning for run {params.run}"
        pyscenic ctx \
        {input.path_adj} \
        {input.path_ranking_db} \
        --annotations_fname {input.path_motif_annotations} \
        --expression_mtx_fname {input.path_loom} \
        --output {output.path_reg} \
        --mask_dropouts \
        --num_workers {threads} > {log}
        """

rule post:
    input:
        path_data=config['input_loc'],
        path_csvs=config['outdir'],
        path_loom="{}/rna.loom".format(config['outdir']),
    output:
        "{outdir}/scenic.h5mu"
    log:
        "{outdir}/logs/post.log"
    benchmark:
        "{outdir}/benchmarks/post.benchmark"
    shell:
        """
        echo "Postprocessing to grn.csv"
        python workflow/scripts/post.py \
        -i {input.path_data} \
        -l {input.path_loom} \
        -c {input.path_csvs} \
        -o {output} > {log}
        """

# Ensure rule execution order
ruleorder: pre > grn > prune > post
