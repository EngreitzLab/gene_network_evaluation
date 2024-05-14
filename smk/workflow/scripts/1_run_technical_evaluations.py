import sys
sys.path.insert(1, snakemake.config['repodir'])
from src.evaluation import compute_explained_variance_ratio, compute_categorical_association

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)

import mudata as mu

def run_technical_evaluations(config, input_):

    # Read mudata
    mdata = mu.read(input_)

    # Create output folder
    output_dir = os.path.join(config['workdir'], 
                              'evaluations/1_technical_evaluations/')
    os.makedirs(output_dir, exists_ok=True)

    # Compute explained variance ratio
    # explained_variance_ratios = \
    # compute_explained_variance_ratio(mdata, n_jobs=1, 
    #                                  prog_key=config['prog_key'], 
    #                                  data_key=config['data_key'], 
    #                                  inplace=False)
    # explained_variance_ratios.to_csv(os.path.join(output_dir, 
    #                                  'explained_variance_ratio_{}_{}.txt'.format(config['prog_key'],
    #                                                                              config['data_key'])),
    #                                  sep='\t')

    # Compute batch association
    for key in config['categorical_keys']:
        results_df, posthoc_df = compute_categorical_association(mdata, 
                                        prog_key=config['prog_key'],                                
                                        categorical_key=key, 
                                        pseudobulk_key=config['pseudobulk_key'],
                                        n_jobs=mdata[config['prog_key']].shape[1],
                                        inplace=False)
        results_df.to_csv(os.path.join(output_dir, 'association_{}.txt'.format(key)), sep='\t')
        posthoc_df.to_csv((os.path.join(output_dir, 'association_posthoc_{}.txt'.format(key))), sep='\t')
    

# Execution (assumes Snakemake)
with open(snakemake.log[0], 'a') as f:
    f.write('\nstderr:\n')
    sys.stderr = sys.stdout = f

    if __name__=='__main__':
        run_technical_evaluations(snakemake.config, 
                                  snakemake.input[0])

