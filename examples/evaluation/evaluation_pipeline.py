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


def setup_logging(path_out):
    log_path = os.path.join(path_out, 'evaluation_pipeline.log')
    if os.path.exists(log_path):
        os.remove(log_path)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[
        logging.FileHandler(log_path, mode='w'),
        logging.StreamHandler()
    ])


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
    setup_logging(path_out)
    
    # Update cNMF key and save
    old_prog_key = io_config['prog_key']
    prog_key = os.path.basename(path_out)
    mdata.mod[prog_key] = mdata.mod.pop(old_prog_key)
    mdata.update()
    mdata.write(os.path.join(path_out, 'eval.h5mu'))

    data_key = io_config['data_key']

    # Run categorical association
    if 'categorical_association' in config:
        logging.info("Running categorical association")
        categorical_association_config = config['categorical_association']
        results_df, posthoc_df = compute_categorical_association(
            mdata, 
            prog_key=prog_key,
            **categorical_association_config,
        )
        results_df.to_csv(os.path.join(path_out, 'categorical_association_results.txt'), sep='\t', index=False) 
        posthoc_df.to_csv(os.path.join(path_out, 'categorical_association_posthoc.txt'), sep='\t', index=False)
    else:
        logging.info("Skipping categorical association, configuration not found")

    # Run perturbation association
    if 'perturbation_association' in config:
        logging.info("Running perturbation association")
        perturbation_association_config = config['perturbation_association']
        if perturbation_association_config.get("groupby_key"):
            groupby_key = perturbation_association_config.pop("groupby_key")
            perturbation_association_df = pd.DataFrame()
            for group in mdata[data_key].obs[groupby_key].unique():
                mdata_ = mdata[mdata[data_key].obs[groupby_key] == group]
                test_stats_df = compute_perturbation_association(
                    mdata_, 
                    prog_key=prog_key,
                    **perturbation_association_config,
                )
                test_stats_df.to_csv(os.path.join(path_out, f'perturbation_association_results_{group}.txt'), sep='\t', index=False)
                perturbation_association_df = pd.concat([perturbation_association_df, test_stats_df])
                perturbation_association_df["group"] = group
            perturbation_association_df.to_csv(os.path.join(path_out, 'perturbation_association_results.txt'), sep='\t', index=False)
        else:
            perturbation_association_config.pop("groupby_key")
            perturbation_association_df = compute_perturbation_association(
                mdata, 
                prog_key=prog_key,
                **perturbation_association_config,
            )
            perturbation_association_df.to_csv(os.path.join(path_out, 'perturbation_association_results.txt'), sep='\t', index=False)
    else:
        logging.info("Skipping perturbation association, configuration not found")

    # Run gene set enrichment analysis
    if 'gene_set_enrichment' in config:
        logging.info("Running gene set enrichment analysis")
        gene_set_enrichment_config = config['gene_set_enrichment']
        libraries = gene_set_enrichment_config.pop('libraries')
        gene_set_enrichment_df = pd.DataFrame()
        for library in libraries:
            pre_res = compute_geneset_enrichment(
                mdata, 
                prog_key=prog_key,
                data_key=data_key,
                library=library,
                **gene_set_enrichment_config,
            )
            if gene_set_enrichment_config["method"] == "fisher":
                pre_res = pre_res.rename(columns={"Term": "term", "P-value": "pval", "Adjusted P-value": "adj_pval", "Odds Ratio": "effect_size", "Genes": "genes"})
            elif gene_set_enrichment_config["method"] == "gsea":
                pre_res = pre_res.rename(columns={"Term": "term", "NOM p-val": "pval", "FDR q-val": "adj_pval", "NES": "effect_size", "Lead_genes": "genes"})
            pre_res['library'] = library
            pre_res.to_csv(os.path.join(path_out, f'geneset_enrichment_{library}.txt'), sep='\t', index=False)
            gene_set_enrichment_df = pd.concat([gene_set_enrichment_df, pre_res])
        gene_set_enrichment_df.to_csv(os.path.join(path_out, 'geneset_enrichment.txt'), sep='\t', index=False)
    else:
        logging.info("Skipping gene set enrichment analysis, configuration not found")

    # Run trait enrichment analysis
    if 'trait_enrichment' in config:
        logging.info("Running trait enrichment analysis")
        trait_enrichment_config = config['trait_enrichment']
        pre_res_trait = compute_trait_enrichment(
            mdata, 
            prog_key=prog_key,
            data_key=data_key,
            **trait_enrichment_config,
        )
        pre_res_trait.to_csv(os.path.join(path_out, 'trait_enrichment.txt'), sep='\t', index=False)
        data = process_enrichment_data(
            enrich_res=pre_res_trait,
            **trait_enrichment_config,
        )
        data.to_csv(os.path.join(path_out, "trait_enrichment_processed.txt"), sep='\t', index=False)
    else:
        logging.info("Skipping trait enrichment analysis, configuration not found")

    # Run motif enrichment analysis
    if 'motif_enrichment' in config:
        logging.info("Running motif enrichment analysis")
        motif_enrichment_config = config['motif_enrichment']
        motif_match_df, motif_count_df, motif_enrichment_df = compute_motif_enrichment(
            mdata, 
            prog_key=prog_key,
            data_key=data_key,
            **motif_enrichment_config,
        )
        motif_match_df.to_csv(os.path.join(path_out, 'motif_enrichment_matches.txt'), sep='\t', index=False)
        motif_count_df.to_csv(os.path.join(path_out, 'motif_enrichment_counts.txt'), sep='\t', index=False)
        motif_enrichment_df.to_csv(os.path.join(path_out, 'motif_enrichment.txt'), sep='\t', index=False)
    else:
        logging.info("Skipping motif enrichment analysis, configuration not found")

    # Run explained variance analysis
    if 'explained_variance' in config:
        logging.info("Running explained variance analysis")
        explained_variance_config = config['explained_variance']
        explained_variance_ratio = compute_explained_variance_ratio(
            mdata, 
            prog_key=prog_key,
            data_key=data_key,
            **explained_variance_config,
        )
        explained_variance_ratio.to_csv(os.path.join(path_out, 'explained_variance_ratio.txt'), sep='\t')
    else:
        logging.info("Skipping explained variance analysis, configuration not found")

    # Save software versions
    import joblib
    import numpy as np
    import scipy
    import sklearn
    import scikit_posthocs as posthocs
    import gseapy
    import pymemesuite

    versions = {
        "evaluation_pipeline_versions": {
            'gene_program_evaluation': '0.0.1',
            'mudata': mudata.__version__,
            'joblib': joblib.__version__,
            'scipy': scipy.__version__,
            'numpy': np.__version__,
            'pandas': pd.__version__,
            'scikit-learn': sklearn.__version__,
            'scikit-posthocs': posthocs.__version__,
            'gseapy': gseapy.__version__,
            'pymemesuite': pymemesuite.__version__,
        }
    }

    with open(os.path.join(path_out, 'software_versions.yml'), 'w') as f:
        yaml.dump(versions, f)

    # Post-flight checks
    expected_files = [
        "trait_enrichment_processed.txt",
        "software_versions.yml",
        "trait_enrichment.txt",
        "perturbation_association_results.txt",
        "motif_enrichment_matches.txt",
        "motif_enrichment_counts.txt",
        "motif_enrichment.txt",
        "geneset_enrichment.txt",
        "explained_variance_ratio.txt",
        "categorical_association_results.txt",
        "categorical_association_posthoc.txt"
    ]
    for file in expected_files:
        if not os.path.exists(os.path.join(path_out, file)):
            logging.error(f"Missing file: {file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run gene program evaluation pipeline.')
    parser.add_argument('--config', type=str, help='Path to the configuration file', required=True)
    args = parser.parse_args()

    main(args.config)
