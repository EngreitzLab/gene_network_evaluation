import os

# Output directory from config
outdir: config['outdir']
currdir = os.getcwd()

rule peak_corr:
    input:
        data=config['input_loc']
    singularity:
        "envs/celloracle.sif"
    output:
        path_all_peaks="{outdir}/all_peaks.csv",
        path_connections="{outdir}/cicero_connections.csv"
    params:
        organism=lambda w: config['organism']
    log:
        "{outdir}/logs/peak_corr.log"
    benchmark:
        "{outdir}/benchmarks/peak_corr.txt"
    shell:
        """
        Rscript workflow/scripts/peak_corr.R {input.data} {params.organism} {output.path_all_peaks} {output.path_connections}
        """

rule tss_annotation:
    input:
        all_peaks="{outdir}/all_peaks.csv",
        connections="{outdir}/cicero_connections.csv"
    singularity:
        "envs/celloracle.sif"
    output:
        "{outdir}/processed_peak_file.csv"
    params:
        organism=lambda w: config['organism'],
        thr_coaccess=lambda w: config['thr_coaccess']
    log:
        "{outdir}/logs/tss_annotation.log"
    benchmark:
        "{outdir}/benchmarks/tss_annotation.txt"
    shell:
         "python workflow/scripts/tss_annotation.py -a {input.all_peaks} -c {input.connections} -o {params.organism} -t {params.thr_coaccess} -p {output}"

rule tf_motif_scan:
    input:
        "{outdir}/processed_peak_file.csv"
    singularity:
        "envs/celloracle.sif"
    output:
        "{outdir}/motifs.celloracle.tfinfo"
    resources:
        mem_mb=32000
    params:
        organism=lambda w: config['organism'],
        fpr=lambda w: config['fpr']
    log:
        "{outdir}/logs/tf_motif_scan.log"
    benchmark:
        "{outdir}/benchmarks/tf_motif_scan.txt"
    shell:
        "python workflow/scripts/tf_motif_scan.py -p {input} -o {params.organism} -f {params.fpr} -t {output}"

rule build_base_grn:
    input:
        "{outdir}/motifs.celloracle.tfinfo"
    singularity:
        "envs/celloracle.sif"
    params:
        thr_motif_score=lambda w: config['thr_motif_score']
    log:
        "{outdir}/logs/build_base_grn.log"
    benchmark:
        "{outdir}/benchmarks/build_base_grn.txt"
    output:
        "{outdir}/base_GRN_dataframe.csv"
    shell:
        "python workflow/scripts/build_base_grn.py -i {input} -t {params.thr_motif_score} -g {output}"

rule build_grn:
    input:
        mdata=config['input_loc'],
        base_grn="{outdir}/base_GRN_dataframe.csv"
    singularity:
        "envs/celloracle.sif"
    log:
        "{outdir}/logs/build_grn.log"
    benchmark:
        "{outdir}/benchmarks/build_grn.txt"
    output:
        "{outdir}/grn.celloracle.links"
    shell:
        "python workflow/scripts/build_grn.py -m {input.mdata} -b {input.base_grn} -l {output}"

rule filter_grn:
    input:
        path_input=config['input_loc'],
        path_peaks="{outdir}/all_peaks.csv",
        path_cicero="{outdir}/cicero_connections.csv",
        path_motifs="{outdir}/motifs.celloracle.tfinfo",
        path_links="{outdir}/grn.celloracle.links",
        path_basegrn="{outdir}/base_GRN_dataframe.csv"
    singularity:
        "envs/celloracle.sif"
    params:
        organism=lambda w: config['organism'],
        thr_coaccess=lambda w: config['thr_coaccess'],
        thr_edge_pval=lambda w: config['thr_edge_pval'],
        thr_top_edges=lambda w: config['thr_top_edges']
    log:
        "{outdir}/logs/filter_grn.log"
    benchmark:
        "{outdir}/benchmarks/filter_grn.txt"
    output:
        path_r2g="{outdir}/r2g.csv",
        path_TF2r="{outdir}/TF2r.csv",
        path_grn="{outdir}/grn.csv",
        path_tri="{outdir}/tri.csv",
        path_mdata="{outdir}/celloracle.h5mu"
    shell:
        "python workflow/scripts/filter_grn.py \
        -i {input.path_input}  \
        -o {params.organism} \
        -k {input.path_peaks} \
        -c {input.path_cicero} \
        -q {params.thr_coaccess} \
        -j {output.path_r2g} \
        -f {input.path_motifs} \
        -e {output.path_TF2r} \
        -l {input.path_links} \
        -p {params.thr_edge_pval} \
        -t {params.thr_top_edges} \
        -g {output.path_grn} \
        -b {input.path_basegrn} \
        -r {output.path_tri} \
        -m {output.path_mdata}"
