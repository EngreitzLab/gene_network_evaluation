import os
import argparse

import numpy as np
import pandas as pd
from anndata import AnnData

import gseapy as gp
from gseapy import Msigdb
from gseapy import Biomart

from joblib import Parallel, delayed
from tqdm.auto import tqdm

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
            gmt = msig.get_gmt(category=library, dbver="2023.2.Hs")
        elif organism=='mouse':
            gmt = msig.get_gmt(category=library, dbver="2023.1.Mm")

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
def perform_prerank(mdata, prog_nam=None, geneset=None,
                    n_jobs=1, prog_key='prog', data_key='rna', 
                    **kwargs):

    # Gene names in human (upper case) format
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
    set_idxs = [mdata[prog_key].uns['genesets'].index(idx) \
                for idx in pre_res.index.values]
    for key, value in mdata[prog_key].uns['gsea_varmap'].items():
        if value == 'Tag %':
            if key=='tag_before':
                mdata[prog_key].varm[key][prog_idx, set_idxs] = \
                pre_res[value].apply(lambda x: x.split('/')[0])\
                .astype(int).values.reshape(1,-1)

            elif key=='tag_after':
                mdata[prog_key].varm[key][prog_idx, set_idxs] = \
                pre_res[value].apply(lambda x: x.split('/')[-1])\
                .astype(int).values.reshape(1,-1)

        else:
            mdata[prog_key].varm[key][prog_idx, set_idxs] = \
            pre_res[value].values.reshape(1,-1)

# Compute gene set enrichment analysis using loadings
def compute_geneset_enrichment(mdata, organism='human', library='h.all', database='msigdb',
                               n_jobs=1, prog_key='prog', data_key='rna', inplace=True, **kwargs):
    
    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata

    # Get geneset
    geneset = get_geneset(organism, library, database)
    mdata[prog_key].uns['genesets'] = list(geneset.keys())

    # Store results in a new anndata
    X_ = np.zeros((mdata[prog_key].shape[1], 
                  len(geneset.keys())))
    X_[:] = np.nan
                                  
    mdata[prog_key].uns['gsea_varmap'] = {'ES':'ES',
                                        'NES':'NES',
                                        'p_values':'NOM p-val',
                                        'FDR':'FDR q-val',
                                        'FWER':'FWER p-val',
                                        'tag_before':'Tag %',
                                        'tag_after':'Tag %',
                                        'percent_gene':'Gene %'}
    for key in mdata[prog_key].uns['gsea_varmap'].keys():
        mdata[prog_key].varm[key] = X_.copy()

    if 'prerank_kwargs' not in kwargs.keys():
        kwargs['prerank_kwargs'] = {}
    # Run GSEA for each program (prerank)
    for prog_nam in tqdm(mdata[prog_key].var_names, 
    desc='Computing GSEA per program', unit='programs'):
        perform_prerank(mdata, 
                        prog_nam=prog_nam, 
                        geneset=geneset,
                        n_jobs=n_jobs,
                        prog_key=prog_key,
                        data_key=data_key,
                        **kwargs['prerank_kwargs']
                        ) 

    # TODO: Run ssGSEA using loadings
    # Conceptually programs represent reduced sample dimensionality

    if not inplace: return mdata[prog_key].uns['gsea']
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj')
    parser.add_argument('-n', '--n_jobs', default=1, typ=int)
    parser.add_argument('-pk', '--prog_key', default='prog', typ=str) 
    parser.add_argument('-dk', '--data_key', default='rna', typ=str)
    parser.add_argument('-og', '--organism', default='Human', choices={'human', 'mouse'}) 
    parser.add_argument('-gs', '--library', default='Reactome_2022', typ=str) 
    parser.add_argument('-gs', '--database', default='enrichr', choices={'msigdb', 'enrichr'}) 
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    compute_geneset_enrichment(args.mudataObj, organism=args.organism, library=args.library,
                               database=arg.database, n_jobs=args.n_jobs, prog_key=args.prog_key, 
                               data_key=args.data_key, inplace=args.output)

