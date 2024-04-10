import sys
sys.path.insert(1, snakemake.config['repodir'])
from src.evaluation import create_opentargets_gwas_benchmarking_data

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)
    
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
