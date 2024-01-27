import os
import gin
import argparse

import mudata
import anndata
from tqdm.auto import tqdm

import logging
logging.basicConfig(level = logging.INFO)

# https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.FactorAnalysis.html
@gin.configurable
def run_factor_analysis_(exp_data, n_components=10, random_state=0, **kwargs):
  
    from scipy import sparse
    from sklearn.decomposition import FactorAnalysis

    # Init FA class
    logging.info('Normalising and log transforming count data is reccomended before running FA.')
    fa = FactorAnalysis(n_components=n_components, random_state=random_state, **kwargs)

    if sparse.issparse(exp_data):
        exp_data = exp_data.toarray()
    
    # Compute decomposition
    decomp = fa.fit_transform(exp_data)

    return fa, decomp

def run_factor_analysis(mdata, prog_key='factor_analysis', 
                        data_key='rna', layer='X', config_path=None, 
                        inplace=True):
    
    # Load method specific parameters
    try: gin.parse_config_file(config_path)
    except: raise ValueError('gin config file could not be found')

    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata

    if layer=='X':
        exp_data = mdata[data_key].X
    else:
        exp_data = mdata[data_key].layers[layer]

    # Compute cNMF and create prog anndata
    fa, decomp = run_factor_analysis_(exp_data)

    # Create new anndata object
    mdata.mod[prog_key]  = anndata.AnnData(X=decomp, obs=mdata[data_key].obs)
    mdata[prog_key].varm['loadings'] = fa.components_
    mdata[prog_key].uns['loadings_genes'] = mdata[data_key].var_names.tolist()

    if not inplace: return mdata[prog_key]

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-pk', '--prog_key', default='factor_analysis', typ=str) 
    parser.add_argument('-dk', '--data_key', default='rna', typ=str) 
    parser.add_argument('--layer', default='X', type=str)
    parser.add_argument('--config_path', default='./factor_analysis_config.gin', type=str)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    import mudata
    mdata = mudata.read(args.mudataObj_path)

    run_factor_analysis(mdata, layer=args.layer, prog_key=args.prog_key, 
                        data_key=args.data_key, inplace=args.output, 
                        config_path=args.config_path)