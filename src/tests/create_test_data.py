from mudata import MuData
from anndata import AnnData

import numpy as np
import pandas as pd
from sklearn.decomposition import FactorAnalysis

import logging
logging.basicConfig(level = logging.INFO)
from tqdm.auto import tqdm

from pathlib import Path

def load_motif_names():
    path = (Path(__file__).parent / './test_data/motif_list.txt').resolve()
    with open(path, 'r') as fil:
        motif_names = fil.readlines()
        motif_names = [nam.strip() for nam in motif_names]
    return motif_names

def load_gene_names():
    path = (Path(__file__).parent / './test_data/gene_list.txt').resolve()
    with open(path, 'r') as fil:
        gene_names = fil.readlines()
        gene_names = [nam.strip() for nam in gene_names]
    return gene_names

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

        means = np.random.uniform(1, 10, num_vars)
    
        poisson_data_ = []
        for mean in means:
            poisson_data_.append(np.random.poisson(mean, size=num_obs))
        
        poisson_data = []
        for i in range(len(poisson_data_)):
            rand_idx = np.random.randint(10, size=3)
            poisson_data.append(np.sum((poisson_data_[rand_idx[0]],
                                        poisson_data_[rand_idx[1]],
                                        poisson_data_[rand_idx[2]]),0))
        data.append(np.array(poisson_data).T.astype(float))
        obs.append(pd.DataFrame([batch]*num_obs, 
                           index=np.arange(num_obs), 
                           columns=['batch']))

    data = np.concatenate(data)
    obs = pd.concat(obs).reset_index(drop=True)
    var = pd.DataFrame(index=load_gene_names()[:num_vars])
    obs['batch'] = obs['batch'].astype(str)
    exp_data = AnnData(X=data, obs=obs, var=var)

    fa = FactorAnalysis(n_components=num_vars)
    fa_data = fa.fit_transform(data)
    prog_data = AnnData(X=fa_data, obs=obs)
    prog_data.varm['loadings'] = fa.components_
    prog_data.uns['var_names'] = exp_data.var_names.tolist()

    mdata = MuData({'rna': exp_data, 'prog': prog_data})

    return mdata