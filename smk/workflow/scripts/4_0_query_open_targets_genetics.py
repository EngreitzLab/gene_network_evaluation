import sys
sys.path.insert(1, snakemake.config['repodir'])
from src.evaluation import run_opentargets_query

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)
    
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
