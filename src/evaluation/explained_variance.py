import os
import argparse

import mudata

import numpy as np
from scipy import sparse
from sklearn.metrics import explained_variance_score, r2_score

from joblib import Parallel, delayed
from tqdm.auto import tqdm

import logging
logging.basicConfig(level = logging.INFO)

def _compute_explained_variance_ratio(mdata, prog_key=None, data_key=None, 
                                      prog_nam=None,layer='X', **kwargs):

    # FIXME: This takes a lot of memory (200G single core on TeloHAEC)
    # https://lightning.ai/docs/torchmetrics/stable/regression/explained_variance.html
    if layer=='X':
        recons = mdata[prog_key][:,prog_nam].X.dot(\
                sparse.csr_matrix(mdata[prog_key][:,prog_nam].varm['loadings']))
    else:
        recons = mdata[prog_key][:,prog_nam].layers[layer].dot(\
                sparse.csr_matrix(mdata[prog_key][:,prog_nam].varm['loadings']))        

    mdata[prog_key].var.loc[prog_nam, 'explained_variance_ratio_{}'.format(layer)] = \
    r2_score(mdata[data_key].X.toarray(), recons.toarray(), **kwargs)

# For explained variance vs r2_score see
# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html
def compute_explained_variance_ratio(mdata, prog_key='prog', data_key='rna', 
                                     layer='X', n_jobs=1, inplace=True, **kwargs):

    """
    Computes the proportion of variance in the data explained by each program.

    ARGS
        mdata : MuData
            mudata object containing anndata of program scores and cell-level metadata.
        prog_key: 
            index for the anndata object (mdata[prog_key]) in the mudata object.
        data_key: str
            index of the genomic data anndata object (mdata[data_key]) in the mudata object.
        layer: str (default: X)
            anndata layer (mdata[data_key].layers[layer]) where the data is stored.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy
       
    RETURNS 
        if not inplace:
            mdata[prog_key].var['explained_variance_ratio_{}'.format(layer)]
            
    """
    # Read in mudata if it is provided as a path
    frompath=False
    if isinstance(mdata, str):
        if os.path.exists(mdata):
            mdata = mudata.read(mdata)
            if inplace:
                logging.warning('Changed to inplace=False since path was provided')
                inplace=False
            frompath=True
        else: raise ValueError('Incorrect mudata specification.')
    
    if not inplace and not frompath:
        mdata = mudata.MuData({prog_key: mdata[prog_key].copy(),
                               data_key: mdata[data_key].copy()})

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
    Parallel(n_jobs=n_jobs, 
             backend='threading')(delayed(_compute_explained_variance_ratio)(mdata, 
                                                                             prog_key=prog_key,
                                                                             data_key=data_key,
                                                                             prog_nam=prog_nam, 
                                                                             layer=layer
                                                       ) \
                                                       for prog_nam in tqdm(mdata[prog_key].var_names,
                                                       desc='Computing explained variance', unit='programs'))

    if not inplace: return (mdata[prog_key].var.loc[:, ['explained_variance_ratio_{}'.format(layer)]])
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('--layer', default='X', type=str)
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-dk', '--data_key', default='rna', type=str) 
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    compute_explained_variance_ratio(mdata, prog_key=args.prog_key, data_key=args.data_key, 
                                     layer=args.layer, n_jobs=args.n_jobs, inplace=args.output)