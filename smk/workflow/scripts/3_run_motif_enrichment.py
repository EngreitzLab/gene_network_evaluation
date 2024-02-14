import sys
sys.path.insert(1, snakemake.config['repodir'])
from src.evaluation import compute_motif_enrichment

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)

import mudata as mu

def run_motif_enrichment(config, input_, coord_file_loc, 
                         seq_class, sig, output_loc):

    # Read mudata
    mdata = mu.read(input_)

    # Compute motif enrichment
    if config['{}_coordinates'.format(seq_class)] is not None:
        compute_motif_enrichment(mdata, 
                                motif_file=coord_file_loc, 
                                seq_file=config['genome_fasta'], 
                                coords_file=config['{}_coordinates'.format(seq_class)], 
                                prog_key=config['prog_key'],
                                data_key=config['data_key'],
                                output_loc=output_loc,
                                sig=sig, 
                                num_genes=None, 
                                n_jobs=-1, 
                                inplace=True)
    else:
        logging.info('No coordinate file was provided.')

# Execution (assumes Snakemake)
with open(snakemake.log[0], 'a') as f:
    f.write('\nstderr:\n')
    sys.stderr = sys.stdout = f

    if __name__=='__main__':
        run_motif_enrichment(snakemake.config, 
                             snakemake.input[0],
                             snakemake.input[1],
                             snakemake.params[0],
                             float(snakemake.params[1]),
                             snakemake.params[2])

