from mudata import MuData
from anndata import AnnData

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from .. import *
from ..compute_motif_enrichment import read_motif_file, read_coords_file, read_sequence_file

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
    var = pd.DataFrame(index=load_gene_names()[:num_vars])
    obs['batch'] = obs['batch'].astype(str)
    exp_data = AnnData(X=data, obs=obs, var=var)

    pca = PCA(n_components=num_vars)
    pca_data = pca.fit_transform(data)
    prog_data = AnnData(X=pca_data, obs=obs)
    prog_data.varm['loadings'] = pca.components_

    mdata = MuData({'rna': exp_data, 'prog': prog_data})

    return mdata

# Functional tests        
def test_explained_variance_ratio(mdata):
    try: compute_explained_variance_ratio(mdata)
    except: raise RuntimeError('Explained variance')
    
    try: assert round(mdata['prog'].var['explained_variance_ratio'].sum(),3)==1
    except: raise AssertionError('Explained variance')

# def test_batch_association(mdata):
#    assert

def test_gene_enrichment(mdata):
    try: compute_geneset_enrichment(mdata)
    except: raise RuntimeError('Gene set enrichment')

def test_motif_enrichment(mdata, coords_file_loc, seq_file_loc):
    
    try: motif_file = read_motif_file((Path(__file__).parent /\
                                       '../../smk/resources/hocomoco_meme.meme').resolve())[:3]
    except: raise RuntimeError('Reading motif file')

    try: coords_file = read_coords_file(coords_file_loc)
    except: raise RuntimeError('Reading coordinate file')

    try: seq_file = read_sequence_file(seq_file_loc)
    except: raise RuntimeError('Reading sequence file')

    try: compute_motif_enrichment(mdata, motif_file, seq_file,
                                  coords_file, output_loc=None)
    except: raise RuntimeError('Motif enrichment')

if __name__=='__main__':

    mdata = create_test_data()

    test_explained_variance_ratio(mdata)
    # test_batch_association(mdata)
    test_gene_enrichment(mdata)
    # test_motif_enrichment(mdata, coords_file, seq_file_loc)