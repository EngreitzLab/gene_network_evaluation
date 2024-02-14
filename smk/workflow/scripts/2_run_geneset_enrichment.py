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

    # Compute geneset enrichment
    if config['library'] is not None:
        compute_geneset_enrichment(mdata, 
                                   prog_key=config['prog_key'], 
                                   data_key=config['rna_key'], 
                                   organism=config['organism'], 
                                   library=config['library'], 
                                   database=config['database'], 
                                   n_jobs=-1, inplace=True)
    else:
        logging.info('No geneset library was provided.')

# Execution (assumes Snakemake)
with open(snakemake.log[0], 'a') as f:
    f.write('\nstderr:\n')
    sys.stderr = sys.stdout = f

    if __name__=='__main__':
        run_geneset_enrichment(snakemake.config, 
                               snakemake.input[0])

