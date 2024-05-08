import os

# Snakemake configuration for the pipeline
configfile: "config.yaml"

# Create the output directory
rule all:
    input:
        expand("{outdir}/run{run}_pyscenic_output.loom", run=range(1, config["num_runs"] + 1)),
        "{outdir}/grnboost2_net.tsv",
        "{outdir}/scenic_network.tsv"

rule setup_outdir:
    output:
        directory("{outdir}")
    shell:
        "mkdir -p {output}"

rule grn_inference:
    input:
        loom_in = config["loom_in"],
        tf_list = config["tf_list"],
        outdir = "{outdir}"
    output:
        adj = "{outdir}/run{run}_adj.tsv"
    params:
        run = "{run}",
        method = config.get("method", "grnboost2")
    threads: 4  # Adjust according to available resources
    shell:
        """
        echo "Running GRN inference for run {params.run}"
        pyscenic grn {input.loom_in} {input.tf_list} -o {output.adj} -m {params.method} --num_workers {threads}
        """

rule enrichment_pruning:
    input:
        adj = "{outdir}/run{run}_adj.tsv",
        ranking = config["ranking"],
        annotation = config["annotation"],
        loom_in = config["loom_in"]
    output:
        reg = "{outdir}/run{run}_reg.csv"
    params:
        run = "{run}"
    threads: 4
    shell:
        """
        echo "Running enrichment and pruning for run {params.run}"
        pyscenic ctx {input.adj} {input.ranking} --annotations_fname {input.annotation} --expression_mtx_fname {input.loom_in} --output {output.reg} --mask_dropouts --num_workers {threads}
        """

rule aucell_calculation:
    input:
        loom_in = config["loom_in"],
        reg = "{outdir}/run{run}_reg.csv"
    output:
        loom_out = "{outdir}/run{run}_pyscenic_output.loom"
    params:
        run = "{run}"
    threads: 4
    shell:
        """
        echo "Calculating AUCell for run {params.run}"
        pyscenic aucell {input.loom_in} {input.reg} --output {output.loom_out} --num_workers {threads}
        """

rule consolidate_grnboost2:
    input:
        expand("{outdir}/run{{run}}_adj.tsv", run=range(1, config["num_runs"] + 1)),
        loom_in = config["loom_in"]
    output:
        "{outdir}/grnboost2_net.tsv"
    params:
        script = "/cellar/users/aklie/projects/igvf/topic_grn_links/opt/grnboost2_pipeline/scripts/grnboost2_output_adapter.py"
    shell:
        "python {params.script} --grnboost2_out_dir {input.outdir} --out_file {output} --loom_file {input.loom_in}"

rule consolidate_scenic:
    input:
        expand("{outdir}/run{{run}}_reg.csv", run=range(1, config["num_runs"] + 1)),
        loom_in = config["loom_in"]
    output:
        "{outdir}/scenic_network.tsv"
    params:
        script = "/cellar/users/aklie/projects/igvf/topic_grn_links/opt/scenic_pipeline/scripts/scenic_output_adapter.py"
    shell:
        "python {params.script} --scenic_out_dir {input.outdir} --loom_file {input.loom_in}"
