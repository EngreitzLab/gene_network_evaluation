import os

# wildcards from config
outdir: config['outdir']
currdir = os.getcwd()

rule all:
    input:
        expand("{}/run{{run}}_pyscenic_output.loom".format(config['outdir']), run=range(1, config["num_runs"] + 1)),
        "{}/grn.csv".format(config['outdir'])

rule setup_outdir:
    output:
        directory(config['outdir'])
    shell:
        """
        echo "Setting up output directory"
        mkdir -p {config['outdir']}
        """

rule pre:
    input:
        data=config['input_loc']
    output:
        path_out="{outdir}/rna.loom"
    params:
        layer=lambda w: config['layer']
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
        path_tf_list=lambda w: os.path.join(currdir, "resources", "tf_lists", os.path.basename(config['tf_list'])),
    output:
        path_adj="{outdir}/run{run}_adj.tsv"
    params:
        run = "{run}",
        method = config.get("inference_method", "grnboost2")
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
        path_ranking_db=lambda w: os.path.join(currdir, "resources", "rankings", os.path.basename(config['rankings_db'])),
        path_motif_annotations=lambda w: os.path.join(currdir, "resources", "annotations", os.path.basename(config['motif_annotations'])),
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

rule aucell:
    input:
        path_loom="{outdir}/rna.loom",
        path_reg="{outdir}/run{run}_reg.csv"
    output:
        path_out="{outdir}/run{run}_pyscenic_output.loom"
    params:
        run="{run}"
    threads: 
        config['n_jobs']
    log:
        "{outdir}/logs/run{run}_aucell.log"
    benchmark:
        "{outdir}/benchmarks/run{run}_aucell.benchmark"
    shell:
        """
        echo "Calculating AUCell for run {params.run}"
        pyscenic aucell \
        {input.path_loom} \
        {input.path_reg} \
        --output {output.path_out} \
        --num_workers {threads} > {log}
        """

rule post:
    input:
        expand("{}/run{{run}}_reg.csv".format(config['outdir']), run=range(1, config["num_runs"] + 1)),
        path_loom="{}/rna.loom".format(config['outdir'])
    output:
        "{}/grn.csv".format(config['outdir'])
    params:
        script="/cellar/users/aklie/projects/igvf/topic_grn_links/opt/scenic_pipeline/scripts/scenic_output_adapter.py"
    shell:
        "python {params.script} --scenic_out_dir {input.outdir} --loom_file {input.path_loom}"

# Ensure rule execution order
ruleorder: setup_outdir > pre > grn > prune > aucell > post