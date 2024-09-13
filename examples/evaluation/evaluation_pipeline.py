#!/usr/bin/env python
import os
import sys
import yaml
import logging
import argparse
import mudata
import pandas as pd

# Import evaluation functions
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from src.evaluation import (
    compute_categorical_association,
    compute_geneset_enrichment,
    compute_trait_enrichment,
    compute_perturbation_association,
    compute_explained_variance_ratio,
    compute_motif_enrichment
)
from src.evaluation.enrichment_trait import process_enrichment_data


def main(config_path):

    # Load the configuration file
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    # Load the input/output configuration
    io_config = config['io']

    # Load mdata
    path_mdata = io_config['path_mdata']
    mdata = mudata.read(path_mdata)

    # Set up output directory
    path_out = io_config['path_out']
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    # Set up logging
    log_path = os.path.join(path_out, 'evaluation_pipeline.log')
    if os.path.exists(log_path):
        os.remove(log_path)
    logging.basicConfig(filename=log_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Update cNMF key and save
    old_prog_key = io_config['prog_key']
    prog_key = os.path.basename(path_out)
    mdata.mod[prog_key] = mdata.mod.pop(old_prog_key)
    mdata.update()
    mdata.write(os.path.join(path_out, f'{prog_key}.h5mu'))

    # Get data key
    data_key = io_config['data_key']

    # Run categorical association
    if 'categorical_association' in config:
        logging.info("Running categorical association")
        categorical_association_config = config['categorical_association']
        categorical_keys = categorical_association_config['categorical_keys']
        categorical_association_config.pop('inplace')
        for key in categorical_keys:
            results_df, posthoc_df = compute_categorical_association(
                mdata, 
                prog_key=prog_key,
                categorical_key=key,
                inplace=False,
                **categorical_association_config,
            )
            results_df.to_csv(os.path.join(path_out, f'{prog_key}_{key}_association_results.txt'), sep='\t', index=False) 
            posthoc_df.to_csv(os.path.join(path_out, f'{prog_key}_{key}_association_posthoc.txt'), sep='\t', index=False)
    else:
        logging.info("Skipping categorical association, configuration not found")

    # Run perturbation association
    if 'perturbation_association' in config:
        logging.info("Running perturbation association")
        perturbation_association_config = config['perturbation_association']
        perturbation_association_config.pop('inplace')
        target_type = "gene" if perturbation_association_config["collapse_targets"] else "guide"
        if perturbation_association_config.get("groupby_key"):
            groupby_key = perturbation_association_config.pop("groupby_key")
            for group in mdata[data_key].obs[groupby_key].unique():
                mdata_ = mdata[mdata[data_key].obs[groupby_key] == group]
                test_stats_df = compute_perturbation_association(
                    mdata_, 
                    prog_key=prog_key,
                    inplace=False,
                    **perturbation_association_config,
                )
                test_stats_df.to_csv(os.path.join(path_out, f'{prog_key}_{target_type}_{groupby_key}_{group}_perturbation_association.txt'), sep='\t', index=False)
        else:
            perturbation_association_config.pop("groupby_key")
            perturbation_association_df = compute_perturbation_association(
                mdata, 
                prog_key=prog_key,
                inplace=False,
                **perturbation_association_config,
            )
            perturbation_association_df.to_csv(os.path.join(path_out, f'{prog_key}_{target_type}_perturbation_association.txt'), sep='\t', index=False)
    else:
        logging.info("Skipping perturbation association, configuration not found")

    # Run gene set enrichment analysis
    if 'gene_set_enrichment' in config:
        logging.info("Running gene set enrichment analysis")
        gene_set_enrichment_config = config['gene_set_enrichment']
        gene_set_enrichment_config.pop('inplace')
        libraries = gene_set_enrichment_config.pop('libraries')
        for library in libraries:
            res = compute_geneset_enrichment(
                mdata, 
                prog_key=prog_key,
                data_key=data_key,
                library=library,
                inplace=False,
                **gene_set_enrichment_config,
            )
            if gene_set_enrichment_config["method"] == "fisher":
                res = res.rename(columns={"Term": "term", "P-value": "pval", "Adjusted P-value": "adj_pval", "Odds Ratio": "enrichment", "Genes": "genes"})
            elif gene_set_enrichment_config["method"] == "gsea":
                res = res.rename(columns={"Term": "term", "NOM p-val": "pval", "FDR q-val": "adj_pval", "NES": "enrichment", "Lead_genes": "genes"})
    
            # Save results
            res.to_csv(os.path.join(path_out, f'{prog_key}_{library}_{gene_set_enrichment_config["method"]}_geneset_enrichment.txt'), sep='\t', index=False)
    else:
        logging.info("Skipping gene set enrichment analysis, configuration not found")

    # Run trait enrichment analysis
    if 'trait_enrichment' in config:
        logging.info("Running trait enrichment analysis")
        trait_enrichment_config = config['trait_enrichment']
        trait_enrichment_config.pop('inplace')
        res_trait = compute_trait_enrichment(
            mdata, 
            prog_key=prog_key,
            data_key=data_key,
            gwas_data=trait_enrichment_config['gwas_data'],
            prog_nam=trait_enrichment_config['prog_nam'],
            library=trait_enrichment_config['library'],
            n_jobs=trait_enrichment_config['n_jobs'],
            inplace=False,
            key_column=trait_enrichment_config['key_column'],
            gene_column=trait_enrichment_config['gene_column'],
            method=trait_enrichment_config['method'],
            loading_rank_thresh=trait_enrichment_config['loading_rank_thresh'],
        )
        if trait_enrichment_config["method"] == "fisher":
            res_trait = res_trait.rename(columns={"Term": "term", "P-value": "pval", "Adjusted P-value": "adj_pval", "Odds Ratio": "effect_size", "Genes": "genes"})
        elif trait_enrichment_config["method"] == "gsea":
            res_trait = res_trait.rename(columns={"Term": "term", "NOM p-val": "pval", "FDR q-val": "adj_pval", "NES": "effect_size", "Lead_genes": "genes"})
        data = process_enrichment_data(
            enrich_res=res_trait,
            metadata=trait_enrichment_config['metadata'],
            pval_col=trait_enrichment_config["pval_col"],
            enrich_geneset_id_col=trait_enrichment_config["enrich_geneset_id_col"],
            metadata_geneset_id_col=trait_enrichment_config["metadata_geneset_id_col"],
            color_category_col=trait_enrichment_config["color_category_col"],
            program_name_col=trait_enrichment_config["program_name_col"],
            annotation_cols=trait_enrichment_config["annotation_cols"],
        )
        data.to_csv(os.path.join(path_out, f"{prog_key}_{trait_enrichment_config['library']}_{trait_enrichment_config['method']}_trait_enrichment.txt"), sep='\t', index=False)
    else:
        logging.info("Skipping trait enrichment analysis, configuration not found")

    # Run motif enrichment analysis
    if 'motif_enrichment' in config:
        logging.info("Running motif enrichment analysis")
        motif_enrichment_config = config['motif_enrichment']
        motif_enrichment_config.pop('inplace')
        loci_files = motif_enrichment_config['loci_files']
        names = motif_enrichment_config['names']
        for loci_file, name in zip(loci_files, names):
            logging.info(f'Running motif enrichment analysis for {loci_file}')
            motif_match_df, motif_count_df, motif_enrichment_df = compute_motif_enrichment(
                mdata, 
                prog_key=prog_key,
                data_key=data_key,
                loci_file=loci_file,
                inplace=False,
                **motif_enrichment_config,
            )
            motif_match_df.to_csv(os.path.join(path_out, f'{prog_key}_enhancer_test_{motif_enrichment_config["correlation"]}_sample_{name}_motif_match.txt'), sep='\t', index=False)
            motif_count_df.to_csv(os.path.join(path_out, f'{prog_key}_enhancer_test_{motif_enrichment_config["correlation"]}_sample_{name}_motif_count.txt'), sep='\t', index=False)
            motif_enrichment_df.to_csv(os.path.join(path_out, f'{prog_key}_enhancer_test_{motif_enrichment_config["correlation"]}_sample_{name}_motif_enrichment.txt'), sep='\t', index=False)
    else:
        logging.info("Skipping motif enrichment analysis, configuration not found")

    # Run explained variance analysis
    if 'explained_variance' in config:
        logging.info("Running explained variance analysis")
        explained_variance_config = config['explained_variance']
        explained_variance_config.pop('inplace')
        explained_variance_ratio = compute_explained_variance_ratio(
            mdata, 
            prog_key=prog_key,
            data_key=data_key,
            inplace=False,
            **explained_variance_config,
        )
        explained_variance_ratio.index = mdata.mod[prog_key].var.index
        explained_variance_ratio.index.name = 'program_name'
        explained_variance_ratio.columns = ["variance_explained_ratio"]
        explained_variance_ratio.to_csv(os.path.join(path_out, f'{prog_key}_variance_explained_ratio.txt'), sep='\t', index=True)
    else:
        logging.info("Skipping explained variance analysis, configuration not found")

    # Save software versions
    import joblib
    import numpy as np
    import scipy
    import sklearn
    import statsmodels
    import scikit_posthocs as posthocs
    import gseapy
    import tangermeme

    versions = {
        "evaluation_pipeline_versions": {
            'gene_program_evaluation': '0.0.1',
            'numpy': np.__version__,
            'pandas': pd.__version__,
            'mudata': mudata.__version__,
            'scipy': scipy.__version__,
            'scikit-learn': sklearn.__version__,
            'scikit-posthocs': posthocs.__version__,
            'statsmodels': statsmodels.__version__,
            'gseapy': gseapy.__version__,  # gene set enrichment analysis
            'tangermeme': tangermeme.__version__,  # motif enrichment analysis
        }
    }

    with open(os.path.join(path_out, 'software_versions.yml'), 'w') as f:
        yaml.dump(versions, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run gene program evaluation pipeline.')
    parser.add_argument('--config', type=str, help='Path to the configuration file', required=True)
    args = parser.parse_args()

    main(args.config)
