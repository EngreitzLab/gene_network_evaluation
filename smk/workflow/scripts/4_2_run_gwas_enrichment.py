import sys
sys.path.insert(1, snakemake.config['repodir'])
from src.evaluation import compute_trait_enrichment

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)
    
#assuming snakemake  
with open(snakemake.log[0], 'a') as f:
    f.write('\nstderr:\n')
    sys.stderr = sys.stdout = f

    if __name__=='__main__':
        gwas_df = compute_trait_enrichment(
                        gwas_data=snakemake.input[0],
                        mdata=snakemake.params[0],
                        prog_key=snakemake.params[1],
                        data_key=snakemake.params[2],
                        prog_nam=None,
                        library='OT_GWAS',
                        n_jobs=snakemake.params[3],
                        inplace=False,
                        key_column='trait_efos',
                        gene_column='gene_name',
                        method='fisher'
        )
        # Create output folder
        output_dir = os.path.join(config['workdir'], 
                                  'evaluations/4_trait_enrichment/')
        os.makedirs(output_dir, exists_ok=True)

        gwas_df.to_csv(os.path.join(output_dir, 'trait_enrichment.txt'), sep='\t')


