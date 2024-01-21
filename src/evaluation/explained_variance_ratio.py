import os
import argparse

import numpy as np
from scipy import sparse
from sklearn.metrics import explained_variance_score

from joblib import Parallel, delayed
from tqdm.auto import tqdm

def _compute_explained_variance_ratio(mdata, prog_nam=None, prog_key=None, 
                                      data_key=None, layer='X', **kwargs):

    # FIXME: This takes a lot of memory (200G single core on TeloHAEC)
    if layer=='X':
        recons = mdata[prog_key][:,prog_nam].X.dot(\
                sparse.csr_matrix(mdata[prog_key][:,prog_nam].varm['loadings']))
    else:
        recons = mdata[prog_key][:,prog_nam].layers[layer].dot(\
                sparse.csr_matrix(mdata[prog_key][:,prog_nam].varm['loadings']))        

    mdata[prog_key].var.loc[prog_nam, 'explained_variance_ratio_{}'.format(layer)] = \
                  explained_variance_score(mdata[data_key].X.toarray(), recons.toarray(),
                                           **kwargs)

# For explained variance vs r2_score see
# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html
def compute_explained_variance_ratio(mdata, n_jobs=1, prog_key='prog', data_key='rna', 
                                     layer='X', inplace=True, **kwargs):
    
    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata

    # Check if same number of vars are present
    try: assert mdata[prog_key].varm['loadings'].shape[1]==mdata[data_key].var.shape[0]
    except: raise ValueError('Different number of features present in data and program loadings')

    if layer=='X':
        if not sparse.issparse(mdata[data_key].X):
            mdata[data_key].X = sparse.csr_matrix(mdata[data_key].X)
    else:
        if not sparse.issparse(mdata[data_key].layers[layer]):
            mdata[data_key].layers[layer] = sparse.csr_matrix(mdata[data_key].layers[layer])
    if not sparse.issparse(mdata[prog_key].X):
        mdata[prog_key].X = sparse.csr_matrix(mdata[prog_key].X)
  
    mdata[prog_key].var['explained_variance_ratio_{}'.format(layer)] = None

    # Run in parallel (max n_jobs=num_progs)
    Parallel(n_jobs=n_jobs, backend='threading')(delayed(_compute_explained_variance_ratio)(mdata, 
                                                       prog_nam=prog_nam, 
                                                       prog_key=prog_key,
                                                       data_key=data_key,
                                                       layer=layer
                                                       ) \
                                                       for prog_nam in tqdm(mdata[prog_key].var_names,
                                                       desc='Computing explained variance', unit='programs'))

    if not inplace: return mdata[prog_key].var.loc[:, 'explained_variance_ratio_{}'.format(layer)]
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-dk', '--data_key', default='rna', type=str) 
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    import mudata
    mdata = mudata.read(args.mudataObj_path)

    compute_explained_variance_ratio(mdata, n_jobs=args.n_jobs, 
                                     prog_key=args.prog_key, data_key=args.data_key, 
                                     inplace=args.output)