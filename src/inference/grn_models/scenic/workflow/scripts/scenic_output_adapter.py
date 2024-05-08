import os
import glob
import logging
import argparse
import time
import pandas as pd
from scipy.stats import ttest_1samp

logging.basicConfig(level=logging.INFO)


def calc_p_value(importances):
    _, p_value = ttest_1samp(importances, 0)
    return p_value


def regulon2sadj(
    regulons,
):
    net_lst = []
    for tf in regulons:
        tf_name = tf.name.split("(")[0]
        tf_targets = tf.gene2weight
        for target, weight in tf_targets.items():
            net_lst.append([tf_name, target, weight])
    net = pd.DataFrame(net_lst, columns=["TF", "target", "importance"])
    return net


def main(args):
    start_time = time.time()
    logging.info("Starting consolidation...")
    reg_csvs = sorted(glob.glob(os.path.join(args.scenic_out_dir, "*.csv")))

    logging.info("Reading regulons...")
    import pandas as pd
    from pyscenic.cli.utils import load_signatures
    all_edges = pd.DataFrame()
    for reg_csv in reg_csvs:
        regulons = load_signatures(reg_csv)
        adj_df = regulon2sadj(regulons)
        all_edges = pd.concat([all_edges, adj_df])
    all_edges.head()
    logging.info(f"Total edges: {len(all_edges)}")

    logging.info("Grouping by source and target and filtering edges...")
    grouped = all_edges.groupby(['TF', 'target'])
    filtered = grouped.filter(lambda x: len(x) > 1)
    logging.info(f"{len(all_edges) - len(filtered)} edges dropped")

    logging.info("Calculating mean importance for each edge...")
    mean_importance = filtered.groupby(['TF', 'target'])['importance'].mean()
    logging.info(f"Total unique edges: {len(mean_importance)}")

    logging.info("Calculating empirical p-value...")
    from tqdm import tqdm
    tqdm.pandas()
    import numpy as np
    TINY = np.finfo(np.float32).tiny
    p_values_series = filtered.groupby(['TF', 'target'])['importance'].progress_apply(calc_p_value)
    p_values = p_values_series.values + TINY

    logging.info("Transforming values...")
    neg_log_p = -np.log10(p_values)
    normalized_importance = (mean_importance - mean_importance.min()) / (mean_importance.max() - mean_importance.min())

    if args.loom_file is not None:
        logging.info("Adding correlation...")
        import scanpy as sc
        from pyscenic.utils import add_correlation
        adata = sc.read_loom(args.loom_file, sparse=True)
        filtered = add_correlation(filtered, adata.to_df())
        mean_corr = filtered.groupby(['TF', 'target'])['rho'].mean()

    consolidated = pd.DataFrame({
        'source': mean_importance.index.get_level_values('TF'),
        'target': mean_importance.index.get_level_values('target'),
        'weight_signed': np.nan,
        'weight_unsigned': mean_importance.values,
        'weight_minmax_normalized': normalized_importance.values,
        'p': p_values,
        '-logp': neg_log_p,
        'description': np.nan,
        'corr': mean_corr.values if args.loom_file is not None else np.nan
    }).reset_index(drop=True)

    # Remove self-loops
    logging.info("Removing self-loops...")
    consolidated = consolidated[consolidated["source"] != consolidated["target"]]
            
    # Save
    output_path = os.path.join(args.scenic_out_dir, args.out_file)
    consolidated.to_csv(output_path, sep="\t", index=False)
    
    end_time = time.time()
    logging.info(f"Process completed. Results saved at {output_path}. Duration: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Consolidate SCENIC regulons into a single network.")
    parser.add_argument('--scenic_out_dir', type=str, required=True, help='Path to the output directory.')
    parser.add_argument('--loom_file', type=str, default=None, help='Path to the loom file for correlation addition.')
    parser.add_argument('--out_file', type=str, default="net.tsv", help='Name of the output file.')
    
    args = parser.parse_args()
    main(args)
