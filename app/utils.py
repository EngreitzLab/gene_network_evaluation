from typing import List, Dict
import re
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_config(config_path):
    with open(config_path, 'r') as config_file:
        return yaml.safe_load(config_file)


def count(categorical_var, count_var, dataframe):
    counts_df = dataframe.value_counts([categorical_var, count_var])
    counts_df = counts_df.groupby(categorical_var).sum()
    counts_df = counts_df.sort_values(ascending=False)
    counts_df = pd.DataFrame(counts_df.reset_index().values,
                             columns=[categorical_var, count_var])
    return counts_df


def count_unique(categorical_var, count_var, dataframe, cummul=False, unique=False):
    counts_df = count(categorical_var, count_var, dataframe)
    new_df = []
    terms = []
    for prog in counts_df[categorical_var].unique():
        terms_ = dataframe.loc[dataframe[categorical_var] == prog, count_var].unique()
        unique_terms = [term for term in terms_ if term not in terms]
        terms.extend(unique_terms)
        new_df.append([prog, len(unique_terms)])
    new_df = pd.DataFrame(new_df, columns=[categorical_var, count_var])
    if cummul:
        new_df[count_var] = new_df[count_var].cumsum()
    if unique:
        return counts_df
    else:
        return new_df


def infer_dashboard_type(keys: List[str]) -> str:
    if len(keys) == 1:
        return "single_run"
    
    base_names = {}
    for key in keys:
        match = re.match(r'([a-zA-Z]+)_?(\d+)?', key)
        if match:
            base_name, num = match.groups()
            if base_name not in base_names:
                base_names[base_name] = []
            if num:
                base_names[base_name].append(int(num))
    
    if len(base_names) == 1:
        return "cross_k"
    elif all(len(nums) == 0 for nums in base_names.values()):
        return "cross_method"
    else:
        return "mixed"


def process_enrichment_data(
    enrich_res,
    metadata,
    pval_col="FDR q-val",
    enrich_geneset_id_col="Term",
    metadata_geneset_id_col="trait_efos",
    color_category_col="trait_category",
    program_name_col="program_name",
    annotation_cols=["trait_reported", "Lead_genes", "study_id", "pmid"]
):

    # Read in enrichment results
    if isinstance(enrich_res, str):
        enrich_df = pd.read_csv(enrich_res)
    elif isinstance(enrich_res, pd.DataFrame):
        enrich_df = enrich_res
    else:
        raise ValueError("enrich_res must be either a pandas DataFrame or a file path to a CSV file.")

    if isinstance(metadata, str):
        metadata_df = pd.read_csv(metadata, compression='gzip', low_memory=False)
    elif isinstance(metadata, pd.DataFrame):
        metadata_df = metadata
    else:
        raise ValueError("metadata must be either a pandas DataFrame or a file path to a CSV file.")

    # Join the enrichment results and the metadata
    enrich_ps = enrich_df.merge(metadata_df, left_on=enrich_geneset_id_col, right_on=metadata_geneset_id_col, how="left")

    # Only keep the relevant columns
    keep_cols = list([enrich_geneset_id_col, pval_col, metadata_geneset_id_col, color_category_col, program_name_col] + annotation_cols)
    enrich_ps = enrich_ps[keep_cols]

    # Sort by P-value
    enrich_ps = enrich_ps.drop_duplicates().sort_values(by=[color_category_col, pval_col])

    # If the input P-value == 0, then replace it with the lowest non-zero P-value in the dataframe
    min_value = enrich_ps.query(f"`{pval_col}` > 0")[pval_col].min()

    # Compute the -log(10) P-value and deal with edge-cases (e.g. P=0, P=1)
    enrich_ps.loc[enrich_ps[pval_col] == 0, pval_col] = min_value  # Replace P=0 with min non-0 p-value
    enrich_ps['-log10(p-value)'] = abs(-1 * np.log10(enrich_ps[pval_col]))

    enrich_ps.reset_index(drop=True, inplace=True)

    return enrich_ps


def filter_and_count(
    data, 
    categorical_var, 
    count_var,
    sig_var, 
    sig_threshold
):
    filtered_data = data[data[sig_var] < sig_threshold]
    count_df = count(categorical_var=categorical_var, count_var=count_var, dataframe=filtered_data)
    unique_data = filtered_data.sort_values(by=sig_var)
    unique_data = unique_data.drop_duplicates(subset=count_var)
    unique_df = count_unique(categorical_var=categorical_var, count_var=count_var, dataframe=unique_data)
    unique_df = unique_df.sort_values(count_var, ascending=False)
    return count_df, unique_df


# Generate a larger categorical colormap using matplotlib
def generate_large_colormap(num_colors):
    cmap = plt.get_cmap('tab20b', num_colors)
    colors = [f'rgba({int(r*255)}, {int(g*255)}, {int(b*255)}, 0.8)' for r, g, b, _ in cmap(np.linspace(0, 1, num_colors))]
    return colors


# Helper function to map categories to colors
def map_categories_to_colors(categories):
    unique_categories = sorted(categories.unique())
    colors = generate_large_colormap(len(unique_categories))
    color_map = {category: color for category, color in zip(unique_categories, colors)}
    return categories.map(color_map).tolist(), color_map