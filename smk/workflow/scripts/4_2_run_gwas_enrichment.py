import sys
sys.path.insert(1, snakemake.config['repodir'])
from src.evaluation import compute_geneset_enrichment_ot_gwas

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)
    
#assuming snakemake  
with open(snakemake.log[0], 'a') as f:
    f.write('\nstderr:\n')
    sys.stderr = sys.stdout = f

    if __name__=='__main__':
        compute_geneset_enrichment_ot_gwas(
            gwas_data=snakemake.input[0],
            mdata=snakemake.params[0],
            prog_key=snakemake.params[1],
            data_key=snakemake.params[2],
            output_file=snakemake.output[0],
            prog_nam=None,
            library='OT_GWAS',
            n_jobs=snakemake.params[3],
            inplace=False,
            key_column='trait_efos',
            gene_column='gene_name',
            method='fisher'
        )

