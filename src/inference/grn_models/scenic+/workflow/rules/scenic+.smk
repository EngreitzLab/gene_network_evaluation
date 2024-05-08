import os

# Output directory from,
outdir: config['outdir']
currdir = os.getcwd()
organism = config['organism']
if organism == "human":
    species = "homo_sapiens"
elif organism == "mouse":
    species = "mus_musculus"

rule pre:
    input:
        data=config['input_loc']
    singularity:
        "envs/scenic+.sif"
    output:
        path_rna="{outdir}/rna.h5ad",
        path_out="{outdir}/pre.h5mu"
    params:
        cluster_key=lambda w: config['cluster_key'],
        organism=lambda w: config['organism'],
        n_topics=lambda w: config['n_topics'],
        n_iter=lambda w: config['n_iter'],
        alpha=lambda w: config['alpha'],
        alpha_by_topic=lambda w: config['alpha_by_topic'],
        eta=lambda w: config['eta'],
        eta_by_topic=lambda w: config['eta_by_topic'],
        ntop=lambda w: config['ntop'],
        scale_factor_impute=lambda w: config['scale_factor_impute'],
        scale_factor_normalize=lambda w: config['scale_factor_normalize'],
        n_cpu=lambda w: config['n_jobs'],
        temp_dir=lambda w: config['scratchdir'],
        random_state=lambda w: config['random_state']
    log:
        "{outdir}/logs/pre.log"
    benchmark:
        "{outdir}/benchmarks/pre.benchmark"
    shell:
        """
        # Run pre for ATAC
        python workflow/scripts/pre_atac.py \
        -i {input.data} \
        -p {outdir} \
        -c {params.cluster_key} \
        -o {params.organism} \
        -n '{params.n_topics}' \
        -t {params.n_iter} \
        -a {params.alpha} \
        -x {params.alpha_by_topic} \
        -e {params.eta} \
        -y {params.eta_by_topic} \
        -nt {params.ntop} \
        -si {params.scale_factor_impute} \
        -sn {params.scale_factor_normalize} \
        -j {params.n_cpu} \
        -d {params.temp_dir} \
        -r {params.random_state} >> {log}

        # Run pre for RNA
        python workflow/scripts/pre_rna.py \
        -i {input.data} \
        -c {params.cluster_key} \
        -o {output.path_rna} >> {log}

        # Create the pre mudata
        scenicplus prepare_data prepare_GEX_ACC \
        --cisTopic_obj_fname {outdir}/cistopic_obj.pkl \
        --GEX_anndata_fname {output.path_rna} \
        --out_file {output.path_out} >> {log}
        """

rule motif_enrichment_cistarget:
    input:
        path_region_sets="{outdir}/region_sets",
        path_ctx_db=lambda w: os.path.join(currdir, "resources", "rankings", os.path.basename(config['rankings_db'])),
        path_motif_annotations=lambda w: os.path.join(currdir, "resources", "annotations", os.path.basename(config['tf_annotations'])),
    output:
        path_out="{outdir}/ctx_results.hdf5",
    params:
        species=lambda w: species,
        fraction_overlap_w_ctx_database=lambda w: config['fraction_overlap_w_ctx_database'],
        auc_threshold=lambda w: config['auc_threshold'],
        nes_threshold=lambda w: config['nes_threshold'],
        rank_threshold=lambda w: config['rank_threshold'],
        annotation_version=lambda w: config['annotation_version'],
        motif_similarity_fdr=lambda w: config['motif_similarity_fdr'],
        orthologous_identity_threshold=lambda w: config['orthologous_identity_threshold'],
        annotations_to_use=lambda w: config['annotations_to_use'],
        temp_dir=lambda w: config['scratchdir'],
    threads: config['n_jobs']
    log:
        "{outdir}/logs/motif_enrichment_cistarget.log"
    benchmark:
        "{outdir}/benchmarks/motif_enrichment_cistarget.benchmark"
    shell:
        """
        scenicplus grn_inference motif_enrichment_cistarget \
        --region_set_folder {input.path_region_sets} \
        --cistarget_db_fname {input.path_ctx_db} \
        --output_fname_cistarget_result {output.path_out} \
        --temp_dir {params.temp_dir} \
        --species {params.species} \
        --fr_overlap_w_ctx_db {params.fraction_overlap_w_ctx_database} \
        --auc_threshold {params.auc_threshold} \
        --nes_threshold {params.nes_threshold} \
        --rank_threshold {params.rank_threshold} \
        --path_to_motif_annotations {input.path_motif_annotations} \
        --annotation_version {params.annotation_version} \
        --motif_similarity_fdr {params.motif_similarity_fdr} \
        --orthologous_identity_threshold {params.orthologous_identity_threshold} \
        --annotations_to_use {params.annotations}
        --n_cpu {params.n_cpu} >> {log}
        """

rule motif_enrichment_dem
    input:
        path_region_sets="{outdir}/region_sets",
        path_dem_db=lambda w: os.path.join(currdir, "resources", "scores", os.path.basename(config['scores_db'])),
        path_genome_annotation=lambda w: os.path.join(currdir, "resources", "genomes", f"{organism}.tsv"),
        path_motif_annotations=lambda w: os.path.join(currdir, "resources", "annotations", os.path.basename(config['tf_annotations'])),
    output:
        path_out="{outdir}/dem_results.hdf5",
    params:
        species=lambda w: species,
        fraction_overlap_w_dem_database=lambda w: config['fraction_overlap_w_dem_database'],
        dem_max_bg_regions=lambda w: config['dem_max_bg_regions'],
        dem_balance_number_of_promoters: lambda w: config['dem_balance_number_of_promoters'],
        dem_promoter_space: lambda w: config['dem_promoter_space'],
        dem_adj_pval_thr: lambda w: config['dem_adj_pval_thr'],
        dem_log2fc_thr: lambda w: config['dem_log2fc_thr'],
        dem_mean_fg_thr: lambda w: config['dem_mean_fg_thr'],
        dem_motif_hit_thr: lambda w: config['dem_motif_hit_thr'],
        direct_annotation: lambda w: config['direct_annotation'],
        extended_annotation: lambda w: config['extended_annotation'],
        annotation_version: lambda w: config['annotation_version'],
        motif_similarity_fdr: lambda w: config['motif_similarity_fdr'],
        orthologous_identity_threshold: lambda w: config['orthologous_identity_threshold'],
        annotations_to_use: lambda w: config['annotations_to_use'],
        temp_dir=lambda w: config['scratchdir'],
        seed=lambda w: config['random_state']
    threads: config['n_jobs']
    log:
        "{outdir}/logs/motif_enrichment_dem.log"
    benchmark:
        "{outdir}/benchmarks/motif_enrichment_dem.benchmark"
    shell:
        """
        scenicplus grn_inference motif_enrichment_dem \
        --region_set_folder {input.path_region_sets} \
        --dem_db_fname {input.path_dem_db} \
        --output_fname_dem_result {output.path_out} \
        --temp_dir {params.temp_dir} \
        --species {params.species} \
        --fraction_overlap_w_dem_database {params.fraction_overlap_w_dem_database} \
        --max_bg_regions {params.max_bg_regions} \
        --genome_annotation {input.path_genome_annotation} \
        --balance_number_of_promoters \
        --promoter_space {params.dem_promoter_space} \
        --adjpval_thr {params.dem_adj_pval_thr} \
        --log2fc_thr {params.dem_log2fc_thr} \
        --mean_fg_thr {params.dem_mean_fg_thr} \
        --motif_hit_thr {params.dem_motif_hit_thr} \
        --path_to_motif_annotations {input.path_motif_annotations} \
        --annotation_version {params.annotation_version} \
        --motif_similarity_fdr {params.motif_similarity_fdr} \
        --orthologous_identity_threshold {params.orthologous_identity_threshold} \
        --annotations_to_use {params.annotations_to_use} \
        --seed {params.seed} \
        --n_cpu {threads} >> {log}
        """
    
rule prepare_menr:
    input:
        path_dem_result="{outdir}/dem_results.hdf5",
        path_ctx_result="{outdir}/ctx_results.hdf5",
        path_mdata="{outdir}/pre.h5mu",
    output:
        path_tf_names="{outdir}/tf_names.txt",
        path_cistromes_direct="{outdir}/cistromes_direct.h5ad",
        path_cistromes_extended="{outdir}/cistromes_extended.h5ad",
    params:
        direct_annotation=lambda w: config['direct_annotation'],
        extended_annotation=lambda w: config['extended_annotation'],
    log:
        "{outdir}/logs/prepare_menr.log"
    benchmark:
        "{outdir}/benchmarks/prepare_menr.benchmark"
    shell:
        """
        scenicplus prepare_data prepare_menr \
        --paths_to_motif_enrichment_results {input.path_dem_result} {input.path_ctx_result} \
        --multiome_mudata_fname {input.path_mdata} \
        --out_file_tf_names {output.tf_names} \
        --out_file_direct_annotation {output.cistromes_direct} \
        --out_file_extended_annotation {output.cistromes_extended} \
        --direct_annotation {params.direct_annotation} \
        --extended_annotation {params.extended_annotation} >> {log}
        """

rule get_search_space:
    input:
        path_pre="{outdir}/pre.h5mu",
        genome_annotation=lambda w: os.path.join(currdir, "resources", "genomes", f"{organism}.tsv"),
        chromsizes=lambda w: os.path.join(currdir, "resources", "chromsizes", f"{organism}.tsv"),
    output:
        search_space="{outdir}/search_space.tsv"
    params:
        upstream=lambda w: config["search_space_upstream"],
        downstream=lambda w: config["search_space_downstream"],
        extend_tss=lambda w: config["search_space_extend_tss"]
    log:
        "{outdir}/logs/get_search_space.log"
    benchmark:
        "{outdir}/benchmarks/get_search_space.benchmark"
    shell:
        """
        scenicplus prepare_data search_spance \
        --multiome_mudata_fname {input.path_pre} \
        --gene_annotation_fname {input.genome_annotation} \
        --chromsizes_fname {input.chromsizes} \
        --out_fname {output.search_space} \
        --upstream {params.upstream} \
        --downstream {params.downstream} \
        --extend_tss {params.extend_tss} >> {log}
        """

rule tf2g:
   input:
        path_pre="{outdir}/pre.h5mu",
        path_tf_names="{outdir}/tf_names.txt",
    output:
        path_tf2g="{outdir}/tf2g.tsv"
    params:
        temp_dir=lambda w: config["scratchdir"],
        method=lambda w: config["tf_to_gene_importance_method"],
        seed=lambda w: config["random_state"]
    threads: config["n_cpu"]
    shell:
        """
        scenicplus grn_inference TF_to_gene \
        --multiome_mudata_fname {input.path_pre} \
        --tf_names {input.path_tf_names} \
        --temp_dir {params.temp_dir} \
        --out_tf_to_gene_adjacencies {output.path_tf2g} \
        --method {params.method} \
        --n_cpu {threads} \
        --seed {params.seed}
        """ 

rule r2g:
    input:
        path_pre="{outdir}/pre.h5mu",
        search_space="{outdir}/search_space.tsv",
    output:
        path_r2g="{outdir}/r2g.tsv"
    params:
        temp_dir=lambda w: config["scratchdir"],
        method_importance=lambda w: config["region_to_gene_importance_method"],
        method_correlation=lambda w: config["region_to_gene_correlation_method"],
    threads: config["n_cpu"]
    shell:
        """
        scenicplus grn_inference region_to_gene \
        --multiome_mudata_fname {input.path_pre} \
        --search_space_fname {input.search_space} \
        --temp_dir {params.temp_dir} \
        --out_region_to_gene_adjacencies {output.path_r2g} \
        --importance_scoring_method {params.method_importance} \
        --correlation_scoring_method {params.method_correlation} \
        --n_cpu {threads}
        """

rule eGRN_direct:
    input:
        path_tf2g="{outdir}/tf2g.tsv",
        path_r2g="{outdir}/r2g.tsv",
        path_cistromes_direct="{outdir}/cistromes_direct.h5ad",
        path_ctx_db=lambda w: os.path.join(currdir, "resources", "rankings", os.path.basename(config['rankings_db'])),
    output:
        path_eRegulons_direct="{outdir}/eRegulon_direct.tsv"
    params:
        temp_dir=lambda w: config["temp_dir"],
        order_regions_to_genes_by=lambda w: config["order_regions_to_genes_by"],
        order_TFs_to_genes_by=lambda w: config["order_TFs_to_genes_by"],
        gsea_n_perm=lambda w: config["gsea_n_perm"],
        quantiles=lambda w: config["quantile_thresholds_region_to_gene"],
        top_n_regionTogenes_per_gene=lambda w: config["top_n_regionTogenes_per_gene"],
        top_n_regionTogenes_per_region=lambda w: config["top_n_regionTogenes_per_region"],
        min_regions_per_gene=lambda w: config["min_regions_per_gene"],
        rho_threshold=lambda w: config["rho_threshold"],
        min_target_genes=lambda w: config["min_target_genes"]
    threads: config["n_cpu"]
    shell:
        """
        scenicplus grn_inference eGRN \
        --TF_to_gene_adj_fname {input.path_tf2g} \
        --region_to_gene_adj_fname {input.path_r2g} \
        --cistromes_fname {input.path_cistromes_direct} \
        --ranking_db_fname {input.path_ctx_db} \
        --eRegulon_out_fname {output.path_eRegulons_direct} \
        --temp_dir {params.temp_dir} \
        --order_regions_to_genes_by {params.order_regions_to_genes_by} \
        --order_TFs_to_genes_by {params.order_TFs_to_genes_by} \
        --gsea_n_perm {params.gsea_n_perm} \
        --quantiles {params.quantiles} \
        --top_n_regionTogenes_per_gene {params.top_n_regionTogenes_per_gene} \
        --top_n_regionTogenes_per_region {params.top_n_regionTogenes_per_region} \
        --min_regions_per_gene {params.min_regions_per_gene} \
        --rho_threshold {params.rho_threshold} \
        --min_target_genes {params.min_target_genes} \
        --n_cpu {threads}
        """

rule eGRN_extended:
    input:
        path_tf2g="{outdir}/tf2g.tsv",
        path_r2g="{outdir}/r2g.tsv",
        path_cistromes_extended="{outdir}/cistromes_extended.h5ad",
        path_dem_db=lambda w: os.path.join(currdir, "resources", "scores", os.path.basename(config['scores_db'])),
    output:
        path_eRegulons_extended="{outdir}/eRegulon_extended.tsv"
    params:
        temp_dir=lambda w: config["temp_dir"],
        order_regions_to_genes_by=lambda w: config["order_regions_to_genes_by"],
        order_TFs_to_genes_by=lambda w: config["order_TFs_to_genes_by"],
        gsea_n_perm=lambda w: config["gsea_n_perm"],
        quantiles=lambda w: config["quantile_thresholds_region_to_gene"],
        top_n_regionTogenes_per_gene=lambda w: config["top_n_regionTogenes_per_gene"],
        top_n_regionTogenes_per_region=lambda w: config["top_n_regionTogenes_per_region"],
        min_regions_per_gene=lambda w: config["min_regions_per_gene"],
        rho_threshold=lambda w: config["rho_threshold"],
        min_target_genes=lambda w: config["min_target_genes"]
    threads: config["n_cpu"]
    shell:
        """
        scenicplus grn_inference eGRN \
        --TF_to_gene_adj_fname {input.path_tf2g} \
        --region_to_gene_adj_fname {input.path_r2g} \
        --cistromes_fname {input.path_cistromes_extended} \
        --ranking_db_fname {input.path_dem_db} \
        --eRegulon_out_fname {output.path_eRegulons_extended} \
        --temp_dir {params.temp_dir} \
        --order_regions_to_genes_by {params.order_regions_to_genes_by} \
        --order_TFs_to_genes_by {params.order_TFs_to_genes_by} \
        --gsea_n_perm {params.gsea_n_perm} \
        --quantiles {params.quantiles} \
        --top_n_regionTogenes_per_gene {params.top_n_regionTogenes_per_gene} \
        --top_n_regionTogenes_per_region {params.top_n_regionTogenes_per_region} \
        --min_regions_per_gene {params.min_regions_per_gene} \
        --rho_threshold {params.rho_threshold}

rule post