from google.cloud import bigquery
import pandas as pd
from scipy.sparse import coo_matrix

### Query Open Targets using BigQuery to obtain GWAS to Gene mappings across many traits

def run_opentargets_query(credentials_path, output_file, min_assoc_loci=10, min_n_cases=1000, min_l2g_score=0.2, study_ids_to_keep=None):
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



### Filter Open Targets query results to a set of high-quality GWAS and L2G scores

def process_json_format_l2g_columns(row, column_name):
    """
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
    