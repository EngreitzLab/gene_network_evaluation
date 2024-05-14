import sys
sys.path.insert(1, snakemake.config['repodir'])
from src.evaluation import compute_geneset_enrichment

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)

import mudata as mu

def run_geneset_enrichment(config, input_):

    # Read mudata
    mdata = mu.read(input_)

    # Create output folder
    output_dir = os.path.join(config['workdir'], 
                              'evaluations/2_geneset_enrichment/')
    os.makedirs(output_dir, exists_ok=True)

    # Compute geneset enrichment
    if config['library'] is not None:
        gsea_df = compute_geneset_enrichment(mdata, 
                                   prog_key=config['prog_key'], 
                                   data_key=config['rna_key'], 
                                   organism=config['organism'], 
                                   library=config['library'], 
                                   database=config['database'], 
                                   n_jobs=-1, inplace=False)
    else:
        logging.info('No geneset library was provided.')
    
    gsea_df.to_csv(os.path.join(output_dir_, 
                                'geneset_enrichment_{}.txt'.format(config['library'])), 
                   sep='\t')

# Execution (assumes Snakemake)
with open(snakemake.log[0], 'a') as f:
    f.write('\nstderr:\n')
    sys.stderr = sys.stdout = f

    if __name__=='__main__':
        run_geneset_enrichment(snakemake.config, 
                               snakemake.input[0])

