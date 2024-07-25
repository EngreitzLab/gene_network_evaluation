import sys
sys.path.insert(1, snakemake.config['repodir'])
from src.evaluation import compute_motif_enrichment

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)

import mudata as mu

def run_motif_enrichment(config, input_, 
                         coord_file_loc, 
                         seq_class, sig):

    # Read mudata
    mdata = mu.read(input_)

    # Create output folder
    output_dir = os.path.join(config['workdir'], 
                              'evaluations/3_motif_enrichment/')
    os.makedirs(output_dir, exists_ok=True)

    # Compute motif enrichment
    if config['{}_coordinates'.format(seq_class)] is not None:
        motif_match_df, motif_count_df, motif_enrichment_df = compute_motif_enrichment(mdata, 
                                motif_file=coord_file_loc, 
                                seq_file=config['genome_fasta'], 
                                coords_file=config['{}_coordinates'.format(seq_class)], 
                                prog_key=config['prog_key'],
                                data_key=config['data_key'],
                                output_loc='evaluations/3_motif_enrichment/{}/'.format(seq_class),
                                sig=sig, 
                                num_genes=None, 
                                n_jobs=-1, 
                                inplace=False)

        motif_match_df.to_csv(os.path.join(output_dir, 'motif_match_{}.txt'.format(seq_class)),
                            sep='\t')
        motif_count_df.to_csv(os.path.join(output_dir, 'motif_count_{}.txt'.format(seq_class)),
                            sep='\t')
        motif_enrichment_df.to_csv(os.path.join(output_dir, 'motif_enrichment_{}.txt'.format(seq_class)),
                                sep='\t')
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
                             float(snakemake.params[1]))

