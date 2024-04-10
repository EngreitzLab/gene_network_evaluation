import pandas as pd

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
    l2g = pd.read_csv(input_file, compression='gzip')
    
    # Convert the trait_efos into just a flat list of EFO IDs rather than json formatted
    l2g['trait_efos'] = l2g.apply(lambda row: process_json_format_l2g_columns(row, "trait_efos"), axis=1)

    # Remove GWAS with lots of EFO IDs
    filtered_l2g = l2g[l2g['trait_efos'].fillna('').apply(lambda x: x.count(',')) <= 2]

    # Remove double quotes present in some trait_reported rows, but not others
    filtered_l2g.loc[:, 'trait_reported'] = filtered_l2g['trait_reported'].str.replace('"', '')

    # Remove unusual traits
    filtered_l2g = filtered_l2g.query(
        "not trait_reported.str.contains(' or |conditional| and | x |pleiotropy', case=False) "
        "and not trait_reported.str.contains('EA', case=True)"
        "and trait_category != 'Uncategorised'"
        "and trait_category != 'phenotype'"
    )

    # Retain the GWAS with the largest sample size by number of cases by EFO group
    filtered_l2g['n_cases'] = filtered_l2g['n_cases'].fillna(filtered_l2g['n_initial'])
    filtered_l2g = (filtered_l2g.assign(rank=filtered_l2g.groupby(["trait_efos"])["n_cases"].rank(method="min", ascending=False))
                 .query("rank == 1")
                 .drop(columns=["rank"])
                 .reset_index(drop=True))

    # Rename the y_proba_full_model to L2G
    filtered_l2g.rename(columns={'y_proba_full_model': 'L2G'}, inplace=True)

    # Filter positions outside the extended MHC region
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

    
#assuming snakemake  
with open(snakemake.log[0], 'a') as f:
    f.write('\nstderr:\n')
    sys.stderr = sys.stdout = f

    if __name__=='__main__':
        filter_open_targets_gwas_query(
            input_file=snakemake.input[0],
            output_file=snakemake.output[0],
            min_l2g_score=snakemake.params[0],
            remove_mhc_region=snakemake.params[1]
        )
