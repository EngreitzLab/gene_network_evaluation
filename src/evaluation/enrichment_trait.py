import os
import argparse

import mudata

import numpy as np
import pandas as pd
from .enrichment_geneset import compute_geneset_enrichment, create_geneset_dict

import logging
logging.basicConfig(level = logging.INFO)


def run_opentargets_query(credentials_path, output_file, min_assoc_loci=10, 
                          min_n_cases=1000, min_l2g_score=0.2, study_ids_to_keep=None):
    """Query Open Targets using BigQuery to obtain GWAS to Gene mappings across many traits"""

    from google.cloud import bigquery
    import pandas as pd
    from scipy.sparse import coo_matrix

    try:
        # Authenticate with BigQuery
        client = bigquery.Client.from_service_account_json(credentials_path)

        # Construct parameterized query
        query = '''
        WITH ranked_genes AS (
            SELECT locus2gene.study_id, 
                locus2gene.chrom, locus2gene.pos, locus2gene.ref, locus2gene.alt, 
                study_metadata.*, 
                genes.gene_name, locus2gene.y_proba_full_model,
                lead_variants.pval,
                ROW_NUMBER() OVER(PARTITION BY locus2gene.study_id, genes.gene_name ORDER BY lead_variants.pval) AS rn
            
            FROM `bigquery-public-data.open_targets_genetics.locus2gene` AS locus2gene
            
            -- Get GWAS metadata
            INNER JOIN `bigquery-public-data.open_targets_genetics.studies` AS study_metadata
            ON locus2gene.study_id = study_metadata.study_id
            
            -- Get HGNC IDs
            INNER JOIN `bigquery-public-data.open_targets_genetics.genes` AS genes
            ON locus2gene.gene_id = genes.gene_id
            
            -- Get lead variant P-values
            INNER JOIN `bigquery-public-data.open_targets_genetics.variant_disease` AS lead_variants
            ON locus2gene.pos = lead_variants.lead_pos
                AND locus2gene.chrom = lead_variants.lead_chrom
                AND locus2gene.study_id = lead_variants.study_id
            '''

        # Add filter conditions
        query += f'''
            WHERE
                -- Remove the "raw" Neale lab results -- I'm not sure what this is
                locus2gene.study_id NOT LIKE '%raw%'
                
                -- Filter to a l2g score threshold
                AND locus2gene.y_proba_full_model > {min_l2g_score}
                
                -- Filter to number of associated loci of at least {min_assoc_loci}
                AND study_metadata.num_assoc_loci >= {min_assoc_loci}
                
                -- Filter to n_cases of at least {min_n_cases}
            '''

        # Optionally filter by study_ids_to_keep if provided
        if study_ids_to_keep and isinstance(study_ids_to_keep, list) and any(isinstance(x, str) for x in study_ids_to_keep):
            query += f'''
                AND locus2gene.study_id IN UNNEST(@study_ids_to_keep)
                '''

        query += '''
        )
        SELECT * FROM ranked_genes WHERE rn = 1;
        '''

        # Set query parameters
        job_config = bigquery.QueryJobConfig()

        if study_ids_to_keep:
            job_config.query_parameters = [bigquery.ArrayQueryParameter("study_ids_to_keep", "STRING", study_ids_to_keep)]

        # Run the query
        query_job = client.query(query, job_config=job_config)

        # Convert the query results to a Pandas dataframe
        l2g = query_job.to_dataframe()
        l2g.sort_values(by=['study_id', 'pval'])

        # Save dataframe to output file
        l2g.to_csv(output_file, compression="gzip", index=False)

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


def process_json_format_l2g_columns(row, column_name):
    """Filter Open Targets query results to a set of high-quality GWAS and L2G scores
    
    Extracts a comma-separated list of IDs from a string representation of a list of dictionaries in a DataFrame row.

    Parameters:
        row (pandas.Series): A row of a pandas DataFrame.
        column_name (str): The name of the column containing the string representation.

    Returns:
        str: A comma-separated list of IDs extracted from the string representation.
             Returns None if the string representation is not properly formatted or if an error occurs during extraction.
    """
    try:
        elements_str = row[column_name]
        start_index = elements_str.find("[{")  # Find the starting index of the list
        end_index = elements_str.find("}]") + 2  # Find the ending index of the list and include the closing bracket
        
        if start_index != -1 and end_index != -1:  # Ensure both start and end indices are found
            elements_list_str = elements_str[start_index:end_index]
            
            # Remove newline characters and convert to a list of dictionaries
            elements_list = eval(elements_list_str.replace('\n', '').replace('array', 'list'))

            ids = [elem['element'] for elem in elements_list]
            ids.sort()
            return ', '.join(ids)
        else:
            return None  # Return None if start or end index not found
    except Exception as e:
        print(f"Error: {e}, Row: {row}")
        return None


def filter_open_targets_gwas_query(input_file, output_file, min_l2g_score=None, remove_mhc_region=True):
    l2g = pd.read_csv(input_file, compression='gzip', low_memory=False)

    # Convert the trait_efos into just a flat list of EFO IDs rather than json formatted
    l2g['trait_efos'] = l2g.apply(lambda row: process_json_format_l2g_columns(row, "trait_efos"), axis=1)

    # Remove GWAS with lots of EFO IDs
    filtered_l2g = l2g[l2g['trait_efos'].fillna('').apply(lambda x: x.count(',')) <= 2]

    # Remove double quotes present in some trait_reported rows, but not others
    filtered_l2g.loc[:, 'trait_reported'] = filtered_l2g['trait_reported'].str.replace('"', '')

    # Remove unusual traits
    filtered_l2g = filtered_l2g.query(
        "not trait_reported.str.contains(' or |conditional| and | x |pleiotropy|interaction|eg:', case=False)"
        "and not trait_reported.str.contains('EA', case=True)"
        "and trait_category != 'Uncategorised'"
        "and not trait_reported.str.contains('intake', case=False)"
        "and not trait_reported.str.contains(' ms]', case=False)"
        "and not trait_reported.str.contains('Protein quantitative trait', case=False)"
        "and not trait_reported.str.contains('adjusted for', case=False)"
        "and not trait_reported.str.contains('diet', case=False)"
    )

    # Retain the GWAS with the largest sample size by number of cases by EFO group
    filtered_l2g.loc[:, 'n_cases'] = filtered_l2g['n_cases'].fillna(filtered_l2g['n_initial'])
    filtered_l2g = (filtered_l2g.assign(rank=filtered_l2g.groupby(["trait_efos"])["n_cases"].rank(method="min", ascending=False))
                 .query("rank == 1")
                 .drop(columns=["rank"])
                 .reset_index(drop=True))

    # Rename the y_proba_full_model to L2G
    filtered_l2g.rename(columns={'y_proba_full_model': 'L2G'}, inplace=True)

    # Filter to positions outside the extended MHC region, which can have poor L2G mapping
    if remove_mhc_region:
        mhc_start = 28510120
        mhc_end = 33480577
        filtered_l2g = filtered_l2g[~((filtered_l2g['chrom'] == '6') & 
                                              (filtered_l2g['pos'] >= mhc_start) & 
                                              (filtered_l2g['pos'] <= mhc_end))]

    # Filter rows based on min L2G score
    if min_l2g_score is not None:
        filtered_l2g = filtered_l2g[filtered_l2g['L2G'] >= min_l2g_score]

    # Filter down the columns to make this DataFrame a bit easier to work with
    desired_columns = ['trait_category', 'trait_efos', 'trait_reported',
                       'gene_name', 'L2G',
                       'chrom', 'pos', 'ref', 'alt', 'pval',
                       'source', 'study_id', 'pmid',
                       'pub_date', 'num_assoc_loci']
    filtered_l2g = filtered_l2g[desired_columns]

    # Sort the DataFrame
    filtered_l2g = filtered_l2g.sort_values(by=['trait_category', 'study_id', 'pval', 'L2G'],
                                            ascending=[True, True, True, False])

    # Save the filtered DataFrame to a CSV file
    filtered_l2g.to_csv(output_file, index=False)
    

def process_enrichment_data(
    enrich_res,
    metadata,
    pval_col="P-value",
    enrich_geneset_id_col="Term",
    metadata_geneset_id_col="trait_efos",
    color_category_col="trait_category",
    program_name_col="program_name",
    annotation_cols=["trait_reported", "Genes", "study_id", "pmid"]
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
    enrich_ps[f'-log10({pval_col})'] = abs(-1 * np.log10(enrich_ps[pval_col]))

    enrich_ps.reset_index(drop=True, inplace=True)

    return enrich_ps


def compute_trait_enrichment(
    mdata, 
    gwas_data, 
    prog_key='prog', 
    prog_name=None, 
    data_key='rna', 
    library='OT_GWAS', 
    n_jobs=1, 
    inplace=False, 
    key_column='trait_efos',
    gene_column='gene_name', 
    method='fisher', 
    min_size=0,
    max_size=2000,
    n_top=500, 
    low_cutoff=-np.inf, 
    **kwargs
):
    """Compute Trait enrichment using open-targets GWAS"""
    #read in gwas data
    if isinstance(gwas_data, str):
        df = pd.read_csv(gwas_data, compression='gzip', low_memory=False)
    elif isinstance(gwas_data, pd.DataFrame):
        df = gwas_data
    else:
        raise ValueError("gwas_data must be either a pandas DataFrame or a file path to a CSV file.")
    
    gmt = create_geneset_dict(df, key_column=key_column, gene_column=gene_column)

    return (compute_geneset_enrichment(
        mdata=mdata, 
        prog_key=prog_key, 
        data_key=data_key, 
        prog_name=prog_name, 
        method=method, 
        organism="human",
        library=library, 
        database=None, 
        user_geneset=gmt,
        min_size=min_size,
        max_size=max_size,
        low_cutoff=low_cutoff,
        n_top=n_top,
        n_jobs=n_jobs, 
        inplace=inplace, 
        **kwargs
    ))


if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('GWAS_data_path')    
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-dk', '--data_key', default='rna', type=str)
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    compute_trait_enrichment(mdata, args.GWAS_data_path, 
                             prog_key=args.prog_key, 
                             data_key=args.data_key, 
                             n_jobs=args.n_jobs, 
                             inplace=args.output)
    