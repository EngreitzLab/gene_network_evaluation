import os

# Output directory from config
outdir: config['outdir']
currdir = os.getcwd()

rule run_cisTopic:
    input:
        data=config['input_loc']
    singularity:
        "envs/scenic+.sif"
    output:
        path_cistopic_obj="{outdir}/cistopic_obj.pkl",
        path_models="{outdir}/models.pkl"
    params:
        organism=lambda w: config['organism'],
        n_topics=lambda w: config['n_topics'],
        n_iter=lambda w: config['n_iter'],
        alpha=lambda w: config['alpha'],
        alpha_by_topic=lambda w: config['alpha_by_topic'],
        eta=lambda w: config['eta'],
        eta_by_topic=lambda w: config['eta_by_topic'],
        ncpu=lambda w: config['n_jobs'],
        temp_dir=lambda w: config['scratchdir'],
        random_state=lambda w: config['random_state']
    log:
        "{outdir}/logs/cistopic.log"
    shell:
        """
        python workflow/scripts/run_cisTopic.py -i {input.data} -o {params.organism} -c {output.path_cistopic_obj} -n '{params.n_topics}' -t {params.n_iter} -a {params.alpha} -x {params.alpha_by_topic} -e {params.eta} -y {params.eta_by_topic} -m {output.path_models} -j {params.ncpu} -d {params.temp_dir} -r {params.random_state}
        """

rule candidate_enhancers:
    input:
        path_cistopic_obj="{outdir}/cistopic_obj.pkl"
    singularity:
        "envs/scenic+.sif"
    output:
        path_otsu_bin_topics="{outdir}/otsu_bin_topics.pkl",
        path_top_bin_topics="{outdir}/top_bin_topics.pkl",
        path_markers="{outdir}/markers.pkl"
    params:
        organism=lambda w: config['organism'],
        ntop=lambda w: config['ntop'],
        scale_factor_impute=lambda w: config['scale_factor_impute'],
        scale_factor_normalize=lambda w: config['scale_factor_normalize']
    log:
        "{outdir}/logs/candidate_enhancers.log"
    shell:
        """
        python workflow/scripts/candidate_enhancers.py -c {input.path_cistopic_obj} -o {output.path_otsu_bin_topics} -n {params.ntop} -t {output.path_top_bin_topics} -s {params.scale_factor_impute} -r {params.scale_factor_normalize} -m {output.path_markers}
        """

rule run_cisTarget:
    input:
        path_otsu_bin_topics="{outdir}/otsu_bin_topics.pkl",
        path_top_bin_topics="{outdir}/top_bin_topics.pkl",
        path_markers="{outdir}/markers.pkl"
    singularity:
        "envs/scenic+.sif"
    output:
        path_motif_enrichment="{outdir}/menr.pkl"
    params:
        organism=lambda w: config['organism'],
        tf_annotations=lambda w: os.path.join(currdir, "resources", "annotations", os.path.basename(config['tf_annotations'])),
        rankings=lambda w: os.path.join(currdir, "resources", "rankings", os.path.basename(config['rankings_db'])),
        scores=lambda w: os.path.join(currdir, "resources", "scores", os.path.basename(config['scores_db'])),
        run_without_promoters=lambda w: config['run_without_promoters'],
        ncpus=lambda w: config['n_jobs'],
        temp_dir=lambda w: config['scratchdir']
    log:
        "{outdir}/logs/cistarget.log"
    shell:
        """
        python workflow/scripts/run_cisTarget.py -o {params.organism} -p {input.path_otsu_bin_topics} -t {input.path_top_bin_topics} -m {input.path_markers} -d {output.path_motif_enrichment} -a {params.tf_annotations} -r {params.rankings} -s {params.scores} -w {params.run_without_promoters} -j {params.ncpus} -x {params.temp_dir}
        """

rule run_scenicplus:
    input:
        path_input=config['input_loc'],
        path_cistopic_obj="{outdir}/cistopic_obj.pkl",
        path_motif_enrichment="{outdir}/menr.pkl"
    singularity:
        "envs/scenic+.sif"
    output:
        path_scenic_obj="{outdir}/scplus_obj.pkl"
    params:
        organism=lambda w: config['organism'],
        ncpus=lambda w: config['n_jobs'],
        temp_dir=lambda w: config['scratchdir'],
        biomart_host=lambda w: config['biomart_host'],
        min_distance_upstream=lambda w: config['min_distance_upstream'],
        max_distance_upstream=lambda w: config['max_distance_upstream'],
        min_distance_downstream=lambda w: config['min_distance_downstream'],
        max_distance_downstream=lambda w: config['max_distance_downstream']
    log:
        "{outdir}/logs/scenicplus.log"
    shell:
        """
        python workflow/scripts/run_scenic+.py -i {input.path_input} -c {input.path_cistopic_obj} -m {input.path_motif_enrichment} -o {params.organism} -s {output.path_scenic_obj}  -j {params.ncpus} -d {params.temp_dir} -b {params.biomart_host} -u {params.min_distance_upstream} -v {params.max_distance_downstream} -w {params.min_distance_upstream} -z {params.max_distance_downstream}
        """

rule filter_links:
    input:
        path_input=config['input_loc'],
        path_scenic_obj="{outdir}/scplus_obj.pkl"
    singularity:
        "envs/scenic+.sif"
    output:
        path_cistromes="{outdir}/cistromes.pkl",
        path_grn="{outdir}/grn.csv",
        path_r2g="{outdir}/r2g.csv",
        path_tri="{outdir}/tri.csv",
        path_mdata="{outdir}/scenic+.h5mu"
    log:
        "{outdir}/logs/filter_links.log"
    shell:
        """
        python workflow/scripts/filter_links.py -i {input.path_input} -s {input.path_scenic_obj} -g {output.path_grn} -r {output.path_r2g} -t {output.path_tri} -m {output.path_mdata} -c {output.path_cistromes}
        """
