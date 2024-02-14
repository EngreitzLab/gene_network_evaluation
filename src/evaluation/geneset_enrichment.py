import os
import argparse

import mudata
import numpy as np
import pandas as pd
from anndata import AnnData

import gseapy as gp
from gseapy import Msigdb
from gseapy import Biomart

from tqdm.auto import tqdm

import logging
logging.basicConfig(level = logging.INFO)

# Convert everything to human format (upper-case)
# Maybe better to shift this to the data loader
def get_idconversion(var_names):

    bm = Biomart()
    # Check top 10
    for i in range(10):
        if var_names[i].lower().startswith('ens'):
            queries ={'ensembl_gene_id': list(var_names)}
            gene_names = \
            queries['external_gene_name'].apply(lambda x: x.upper()).values
        else:
            gene_names = [nam.upper() for nam in var_names]

    return gene_names

# Either MsigDB or Enrichr
# MsigDB https://www.gsea-msigdb.org/gsea/msigdb
# Enrichr https://maayanlab.cloud/Enrichr/#libraries
def get_geneset(organism='human', 
                library='h.all', 
                database='msigdb'):

    if database=='msigdb':
        msig = Msigdb()
        if organism=='human':
            gmt = msig.get_gmt(category=library, dbver='2023.2.Hs')
        elif organism=='mouse':
            gmt = msig.get_gmt(category=library, dbver='2023.1.Mm')

        if gmt is None:
            raise ValueError('Library does not exist')

    elif database=='enrichr':
            gmt = gp.get_library(name=library, 
                                 organism=organism.capitalize())

    return gmt

# TODO: ssGSEA using loading matrix -> get topic wise enrichment scores
def perform_ssGSEA():
    raise NotImplementedError()

# pre-ranked GSEA using loadings
def perform_prerank(mdata, prog_key='prog', data_key='rna',
                    prog_nam=None, geneset=None, library=None,
                    n_jobs=1, **kwargs):

    # Gene names in human (upper case) format
    if 'var_names' in mdata[prog_key].uns.keys():
        gene_names = get_idconversion(mdata[prog_key].uns['var_names'])
    else:
        try: assert mdata[prog_key].varm['loadings'].shape[1]==mdata[data_key].var.shape[0]
        except: raise ValueError('Different number of genes present in data and program loadings')
        gene_names = get_idconversion(mdata[data_key].var_names)

    # Gene scores for each program
    loadings = pd.DataFrame(data=mdata[prog_key][:, 
                                       prog_nam].varm['loadings'].flatten(),
                            index=gene_names)     

    # Defaults - pass kwargs if you want this done differently
    #threads=4, min_size=5, max_size=1000, permutation_num=1000, 
    #outdir=None, seed=0, verbose=True  

    # Compute enrichment using loadings
    # Parallelize over genesets
    pre_res = gp.prerank(rnk=loadings,
                         gene_sets=geneset,
                         threads=n_jobs,
                         **kwargs).res2d

    pre_res.set_index('Term', inplace=True)

    pre_res['Gene %'] = pre_res['Gene %'].apply(lambda x:x[:-1]).astype(float)

    # Store reports
    prog_idx = mdata[prog_key].var.index.get_loc(prog_nam)
    set_idxs = [mdata[prog_key].uns['genesets_{}'.format(library)].index(idx) \
                for idx in pre_res.index.values]
    for key, value in mdata[prog_key].uns['gsea_varmap_{}'.format(library)].items():
        if value == 'Tag %':
            if key=='tag_before_{}'.format(library):
                mdata[prog_key].varm[key][prog_idx, set_idxs] = \
                pre_res[value].apply(lambda x: x.split('/')[0])\
                .astype(int).values.reshape(1,-1)

            elif key=='tag_after_{}'.format(library):
                mdata[prog_key].varm[key][prog_idx, set_idxs] = \
                pre_res[value].apply(lambda x: x.split('/')[-1])\
                .astype(int).values.reshape(1,-1)

        else:
            mdata[prog_key].varm[key][prog_idx, set_idxs] = \
            pre_res[value].values.reshape(1,-1)

# Compute gene set enrichment analysis using loadings
def compute_geneset_enrichment(mdata, prog_key='prog', data_key='rna', 
                               organism='human', library='Reactome_2022', 
                               database='enrichr', n_jobs=1, inplace=True, 
                               **kwargs):

    """
    Perform GSEA using loadings are preranked gene lists.

    ARGS
        mdata : MuData
            mudata object containing anndata of program scores and cell-level metadata.
        prog_key: 
            index for the anndata object (mdata[prog_key]) in the mudata object.
        data_key: str
            index of the genomic data anndata object (mdata[data_key]) in the mudata object.
        organism: {'human', 'mouse'} (default: 'human')
            species to which the sequencing data was aligned to.
        library: str (default: Reactome_2022)
            gene-set library to use for computing enrichment.
            MsigDB https://www.gsea-msigdb.org/gsea/msigdb
            Enrichr https://maayanlab.cloud/Enrichr/#libraries
        database: {'msigdb', 'enrichr'} (default: 'enrichr')
            database of gene-set libraries to use.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy
       
    RETURNS 
        if not inplace:
            mdata[prog_key].uns['gsea_varmap_{}'.format(library)].keys(), 
            mdata[prog_key].varm,
            mdata[prog_key].uns['genesets_{}'.format(library)]        

    """
    
    if not inplace:
        mdata = mudata.MuData({prog_key: mdata[prog_key].copy(),
                               data_key: mdata[data_key].copy()})

    # Get geneset
    geneset = get_geneset(organism, library, database)
    mdata[prog_key].uns['genesets_{}'.format(library)] = list(geneset.keys())

    # Store results in a new anndata
    X_ = np.zeros((mdata[prog_key].shape[1], 
                  len(geneset.keys())))
    X_[:] = np.nan

    logging.info('If genes are tied in rankings there order is arbitrary. Check log for warnings.')                     
    mdata[prog_key].uns['gsea_varmap_{}'.format(library)] = {'ES_{}'.format(library):'ES',
                                                             'NES_{}'.format(library):'NES',
                                                             'p_values_{}'.format(library):'NOM p-val',
                                                             'FDR_{}'.format(library):'FDR q-val',
                                                             'FWER_{}'.format(library):'FWER p-val',
                                                             'tag_before_{}'.format(library):'Tag %',
                                                             'tag_after_{}'.format(library):'Tag %',
                                                             'percent_gene_{}'.format(library):'Gene %'}
    for key in mdata[prog_key].uns['gsea_varmap_{}'.format(library)].keys():
        mdata[prog_key].varm[key] = X_.copy()

    if 'prerank_kwargs' not in kwargs.keys():
        kwargs['prerank_kwargs'] = {}
    # Run GSEA for each program (prerank)
    for prog_nam in tqdm(mdata[prog_key].var_names, 
    desc='Computing GSEA per program', unit='programs'):
        perform_prerank(mdata, 
                        prog_nam=prog_nam, 
                        geneset=geneset,
                        library=library,
                        n_jobs=n_jobs,
                        prog_key=prog_key,
                        data_key=data_key,
                        **kwargs['prerank_kwargs']
                        ) 

    # TODO: Run ssGSEA using loadings
    # Conceptually programs represent reduced sample dimensionality
    if not inplace: return (mdata[prog_key].uns['gsea_varmap_{}'.format(library)].keys(), 
                            mdata[prog_key].varm,
                            mdata[prog_key].uns['genesets_{}'.format(library)])
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-og', '--organism', default='Human', choices={'human', 'mouse'}) 
    parser.add_argument('-gs', '--library', default='Reactome_2022', type=str) 
    parser.add_argument('-db', '--database', default='enrichr', choices={'msigdb', 'enrichr'}) 
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-dk', '--data_key', default='rna', type=str)
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    compute_geneset_enrichment(mdata, prog_key=args.prog_key, data_key=args.data_key, 
                               organism=args.organism, library=args.library, database=arg.database, 
                               n_jobs=args.n_jobs, inplace=args.output)

