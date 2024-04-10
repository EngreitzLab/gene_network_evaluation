from google.cloud import bigquery
import pandas as pd

#setup basic logging
import sys
import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)

def run_opentargets_query(credentials_path, output_file, min_assoc_loci, min_n_cases, min_l2g_score, study_ids_to_keep=None):
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
    
with open(snakemake.log[0], 'a') as f:
    f.write('\nstderr:\n')
    sys.stderr = sys.stdout = f

    if __name__=='__main__':
        run_opentargets_query(credentials_path = snakemake.input[0],
                      output_file = snakemake.output[0],
                      min_assoc_loci = snakemake.params[0],
                      min_n_cases = snakemake.params[1],
                      min_l2g_score = snakemake.params[2],
                      study_ids_to_keep = snakemake.params[3])
