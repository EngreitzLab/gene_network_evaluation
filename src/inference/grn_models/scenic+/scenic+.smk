rule run_cisTopic:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu"
    singularity:
        "envs/scenic+.sif"
    benchmark:
        "benchmarks/scenic+/{dataset}.{trajectory}.run_cisTopic.txt"
    output:
        path_cistopic_obj="resources/{dataset}/{trajectory}/scenic+/cistopic_obj.pkl",
        path_models="resources/{dataset}/{trajectory}/scenic+/models.pkl"
    params:
        organism=lambda w: config[w.dataset]['organism'],
        n_topics=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['n_topics'],
        n_iter=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['n_iter'],
        alpha=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['alpha'],
        alpha_by_topic=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['alpha_by_topic'],
        eta=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['eta'],
        eta_by_topic=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['eta_by_topic']
    shell:
        """
        python scripts/scenic+/run_cisTopic.py -i {input.data} -o {params.organism} -c {output.path_cistopic_obj} -n '{params.n_topics}' -t {params.n_iter} -a {params.alpha} -x {params.alpha_by_topic} -e {params.eta} -y {params.eta_by_topic} -m {output.path_models}
        """

rule candidate_enhancers:
    input:
        path_cistopic_obj="resources/{dataset}/{trajectory}/scenic+/cistopic_obj.pkl"
    singularity:
        "envs/scenic+.sif"
    benchmark:
        "benchmarks/scenic+/{dataset}.{trajectory}.candidate_enhancers.txt"
    output:
        path_otsu_bin_topics="resources/{dataset}/{trajectory}/scenic+/otsu_bin_topics.pkl",
        path_top_bin_topics="resources/{dataset}/{trajectory}/scenic+/top_bin_topics.pkl",
        path_markers="resources/{dataset}/{trajectory}/scenic+/markers.pkl"
    params:
        organism=lambda w: config[w.dataset]['organism'],
        ntop=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['ntop'],
        scale_factor_impute=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['scale_factor_impute'],
        scale_factor_normalize=lambda w: config[w.dataset]['trajectories'][w.trajectory]['scenic+']['scale_factor_normalize']
    shell:
        """
        python scripts/scenic+/candidate_enhancers.py -c {input.path_cistopic_obj} -o {output.path_otsu_bin_topics} -n {params.ntop} -t {output.path_top_bin_topics} -s {params.scale_factor_impute} -r {params.scale_factor_normalize} -m {output.path_markers}
        """

rule run_cisTarget:
    input:
        path_otsu_bin_topics="resources/{dataset}/{trajectory}/scenic+/otsu_bin_topics.pkl",
        path_top_bin_topics="resources/{dataset}/{trajectory}/scenic+/top_bin_topics.pkl",
        path_markers="resources/{dataset}/{trajectory}/scenic+/markers.pkl"
    singularity:
        "envs/scenic+.sif"
    benchmark:
        "benchmarks/scenic+/{dataset}.{trajectory}.run_cisTarget.txt"
    output:
        path_motif_enrichment="resources/{dataset}/{trajectory}/scenic+/menr.pkl"
    params:
        organism=lambda w: config[w.dataset]['organism'],
    shell:
        """
        python scripts/scenic+/run_cisTarget.py -o {params.organism} -p {input.path_otsu_bin_topics} -t {input.path_top_bin_topics} -m {input.path_markers} -d {output.path_motif_enrichment}
        """

rule run_scenicplus:
    input:
        path_input="resources/{dataset}/{trajectory}/mdata.h5mu",
        path_cistopic_obj="resources/{dataset}/{trajectory}/scenic+/cistopic_obj.pkl",
        path_motif_enrichment="resources/{dataset}/{trajectory}/scenic+/menr.pkl"
    singularity:
        "envs/scenic+.sif"
    benchmark:
        "benchmarks/scenic+/{dataset}.{trajectory}.run_scenic+.txt"
    output:
        path_scenic_obj="resources/{dataset}/{trajectory}/scenic+/scplus_obj.pkl",
        path_grn="resources/{dataset}/{trajectory}/scenic+/grn.csv",
        path_r2g="resources/{dataset}/{trajectory}/scenic+/r2g.csv",
        path_tri="resources/{dataset}/{trajectory}/scenic+/tri.csv"
    params:
        organism=lambda w: config[w.dataset]['organism'],
    shell:
        """
        python scripts/scenic+/run_scenic+.py -i {input.path_input} -c {input.path_cistopic_obj} -m {input.path_motif_enrichment} -o {params.organism} -s {output.path_scenic_obj} -g {output.path_grn} -r {output.path_r2g} -t {output.path_tri}
        """