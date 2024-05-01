import os

# organism from config
organism = config['organism']
if organism == "human":
    genome = "hg38"
elif organism == "mouse":
    genome = "mm10"
else:
    raise ValueError("Unknown organism: {organism}")
print(f"Organism: {organism}, Genome: {genome}")

# Output directory from config
outdir: config['outdir']
currdir = os.getcwd()

rule r2g_R:
    input:
        data=config['input_loc']
    singularity:
        "envs/celloracle.sif"
    output:
        path_all_peaks="{outdir}/all_peaks.csv",
        path_connections="{outdir}/cicero_connections.csv",
        path_cicero_output="{outdir}/cicero_output.rds"
    params:
        organism=lambda w: config['organism'],
        binarize=lambda w: config['binarize'],
        dim_reduction_key=lambda w: config['dim_reduction_key'],
        k=lambda w: config['k'],
        window=lambda w: config['window'],
        seed=lambda w: config['random_state']
    log:
        "{outdir}/logs/r2g.log"
    benchmark:
        "{outdir}/benchmarks/r2g.txt"
    shell:
        """
        Rscript workflow/scripts/r2g.R \
        {input.data} \
        {params.organism} \
        {params.binarize} \
        {params.dim_reduction_key} \
        {params.k} \
        {params.window} \
        {output.path_all_peaks} \
        {output.path_connections} \
        {output.path_cicero_output} \
        {params.seed}
        """

rule r2g_py:
    input:
        data=config['input_loc'],
        all_peaks="{outdir}/all_peaks.csv",
        connections="{outdir}/cicero_connections.csv"
    singularity:
        "envs/celloracle.sif"
    output:
        "{outdir}/r2g.csv"
    params:
        organism=lambda w: config['organism'],
        coaccess_thr=lambda w: config['coaccess_thr']
    log:
        "{outdir}/logs/r2g_py.log"
    benchmark:
        "{outdir}/benchmarks/r2g_py.txt"
    shell:
         "python workflow/scripts/p2g.py \
         -d {input.data} \
         -a {input.all_peaks} \
         -c {input.connections} \
         -o {params.organism} \
         -t {params.coaccess_thr} \
         -p {output}"

rule tf2r:
    input:
        data=config['input_loc'],
        r2g="{outdir}/r2g.csv",
        gdir="resources/genomes/{organism}"
    singularity:
        "envs/celloracle.sif"
    output:
        "{outdir}/tf2r.csv"
    resources:
        mem_mb=32000
    params:
        organism=lambda w: config['organism'],
        fpr=lambda w: config['fpr'],
        blen=lambda w: config['blen'],
        tfb_thr=lambda w: config['tfb_thr']
    log:
        "{outdir}/logs/tf2r.log"
    benchmark:
        "{outdir}/benchmarks/tf2r.txt"
    shell:
        "python workflow/scripts/tf2r.py \
        -d {input.data} \
        -r {input.r2g} \
        -o {params.organism} \
        -f {params.fpr} \
        -b {params.blen} \
        -t {params.tfb_thr} \
        -t {output}"

rule grn:
    input:
        data=config['input_loc'],
        path_r2g="{outdir}/r2g.csv",
        path_tf2r="{outdir}/tf2r.csv"
    singularity:
        "envs/celloracle.sif"
    params:
        cluster_key=lambda w: config['cluster_key'],
        layer=lambda w: config['layer'],
        alpha=lambda w: config['alpha'],
        bagging_number=lambda w: config['bagging_number'],
    log:
        "{outdir}/logs/grn.log"
    benchmark:
        "{outdir}/benchmarks/grn.txt"
    output:
        "{outdir}/grn.csv"
    shell:
        "python workflow/scripts/grn.py \
        -d {input.data} \
        -r {input.path_r2g} \
        -t {input.path_tf2r} \
        -c {params.cluster_key} \
        -l {params.layer} \
        -a {params.alpha} \
        -b {params.bagging_number} \
        -o {output}"

rule post:
    input:
        data=config['input_loc'],
        path_r2g="{outdir}/r2g.csv",
        path_tf2r="{outdir}/tf2r.csv",
        path_grn="{outdir}/grn.csv"
    singularity:
        "envs/celloracle.sif"
    log:
        "{outdir}/logs/post.log"
    benchmark:
        "{outdir}/benchmarks/post.txt"
    output:
        path_mdata="{outdir}/celloracle.h5mu"
    shell:
        "python workflow/scripts/post.py \
        -d {input.data} \
        -r {input.path_r2g} \
        -t {input.path_tf2r} \
        -g {input.path_grn} \
        -o {output.path_mdata}"
