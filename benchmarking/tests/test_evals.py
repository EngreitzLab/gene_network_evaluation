from mudata import MuData
from anndata import AnnData

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from ..evals import *

import logging
logging.basicConfig(level = logging.INFO)
from tqdm.auto import tqdm

# Generate multi-sample test dataset
def create_test_data(num_batches=None, num_obs=None, 
                     num_vars=None, means=None, cov=None):

    # Params for data gen
    if num_batches is None:
        num_batches=5
        logging.info('Setting num_batches to {}'.format(num_batches))
    if num_obs is None:
        num_obs=100
        logging.info('Setting num_obs to {}'.format(num_obs))
    if num_vars is None:
        num_vars=15
        logging.info('Setting num_vars to {}'.format(num_vars))

    # Generate data
    data = []
    obs = []   
    for batch in tqdm(range(num_batches), desc='Generating batch', unit='batches'):

        means = np.random.uniform(-5, 5, num_vars)
    
        cov_ = np.random.uniform(0.05, 0.5, 
                                 int((num_vars**2-num_vars)/2))
        cov = np.zeros((num_vars, num_vars))

        cov[np.tril_indices(num_vars, -1)] = cov_
        cov[np.triu_indices(num_vars, 1)] = cov_
        cov[np.diag_indices(num_vars)] = 1

        data.append(np.random.multivariate_normal(means, cov, num_obs))
        obs.append(pd.DataFrame([batch]*num_obs, 
                           index=np.arange(num_obs), 
                           columns=['batch']))

    data = np.concatenate(data)
    obs = pd.concat(obs).reset_index(drop=True)
    obs['batch'] = obs['batch'].astype(str)
    exp_data = AnnData(X=data, obs=obs)

    pca = PCA(n_components=num_vars)
    pca_data = pca.fit_transform(data)
    prog_data = AnnData(X=pca_data, obs=obs)
    prog_data.varm['loadings'] = pca.components_

    mdata = MuData({'rna': exp_data, 'prog': prog_data})

    return mdata
        
def test_explained_variance_ratio(mdata):
    assert round(mdata['prog'].var['explained_variance_ratio'].sum(),3)==1

# def test_batch_association(mdata):
#     assert

if __name__=='__main__':

    mdata = create_test_data()

    test_explained_variance_ratio(mdata)
    # test_batch_association(mdata)