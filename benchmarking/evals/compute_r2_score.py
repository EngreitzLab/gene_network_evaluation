import os
import argparse

import numpy as np
from scipy import sparse
from sklear.metrics import r2_score

from joblib import Parallel, delayed
from tqdm.auto import tqdm

# R2_score in sklearn accounts for non-zero mean residuals
# This is expected for component wise reconstruction
def _compute_r2_score(mudata, prog_nam=None, 
                      prog_key=None data_key=None,
                      **kwargs):

    recons = np.dot(mudata[prog_key][prog_nam].X, 
                    mudata[prog_key][prog_nam].varm['loadings'])

    mudata[prog_key].var.loc[prog_nam, 'r2_score'] = r2_score(mudata[data_key].X, 
                                                              recons, kwargs)

# For explained variance vs r2_score see
# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html
def compute_r2_score(mudata, prog_key='prog', data_key='rna', 
	                 inplace=True, **kwargs):
    
    #TODO: Don't copy entire mudata only relevant Dataframe
    mudata = mudata.copy() if not inplace else mudata

    if sparse.issparse(mudata[data_key].X):
        mudata[data_key].X = mudata[data_key].X.toarray()
  
    mudata[prog_key].var['r2_score'] = None

    # Run in parallel (max n_jobs=num_progs)
    Parallel(n_jobs=n_jobs)(delayed(_compute_r2_score)(mudata, 
                                                       prog_nam=prog_nam, 
                                                       prog_key=prog_key,
                                                       data_key=data_key,
                                                       **kwargs
                                                       ) \
                                                       for prog_nam in tqdm(mudata[prog_key].var_names,
                                                       desc='Computing r2', unit='programs'))

    if not inplace: return mudata[prog_key].var.loc[:, 'r2_score']
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

        parser.add_argument('mudataObj')
        parser.add_argument('-n', '--n_jobs', default=1, typ=int)
        parser.add_argument('-pk', '--prog_key', default='prog', typ=str) 
        parser.add_argument('-dk', '--data_key', default='rna', typ=str) 
        parser.add_argument('--output', action='store_false') 

        args = parser.parse_args()
    
        compute_r2_score(args.mudataObj, prog_key=args.prog_key, 
                         data_key=args.data_key, inplace=args.output)