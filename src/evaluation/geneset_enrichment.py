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


def create_geneset_dict(df, key_column="trait_efos", gene_column="gene_name"):
    """
    Create a geneset dictionary from a pandas DataFrame.

    This function iterates over the rows of the DataFrame and constructs a geneset dictionary,
    where keys are unique values from the specified key column, and values are lists of gene names
    associated with each key.

    Parameters:
    - df : pandas DataFrame
        A pandas dataframe that you wish to convert to a dictionary for geneset enrichment.
        Useful for taking in OpenTargets GWAS L2G query results.
    - key_column : str, optional (default: "trait_efos")
        The name of the column in the DataFrame to use as the key for the geneset dictionary.
    - gene_column : str, optional (default: "gene_name")
        The name of the column in the DataFrame containing the gene names.

    Returns:
    - geneset_dict : dict
        A dictionary where keys are unique values from the key column (the name of each geneset),
        and values are lists of gene names.
    """
    geneset_dict = {}
    for index, row in df.iterrows():
        key = row[key_column]
        gene = row[gene_column]
        if key not in geneset_dict:
            geneset_dict[key] = []
        geneset_dict[key].append(gene)
    return geneset_dict


# Convert everything to human format (upper-case)
# Maybe better to shift this to the data loader
def get_idconversion(var_names):
    bm = Biomart()

    # Initialize gene_names list
    gene_names = []

    # Check top 10
    for i in range(min(10, len(var_names))):
        if var_names[i].lower().startswith('ens'):
            # Assuming you are querying Biomart here, but you didn't provide the code
            queries = {'ensembl_gene_id': list(var_names)}
            gene_names = queries['external_gene_name'].apply(lambda x: x.upper()).values
            break  # Exit loop since we found 'ens' prefix
        
        elif ':ens' in var_names[i].lower():
            # Modify gene names to keep only the part before ':' in cases where the format is
            # GENE:ENSG0000000000
            gene_names = [name.split(':')[0].upper() for name in var_names]
            break  # Exit loop since we found ':ens' substring
        
        else:
            # Convert all names to uppercase
            gene_names = [name.upper() for name in var_names]

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
                               user_geneset=None,
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
        user_geneset: dict (default: None)
            a user-defined geneset to be used instead of downloading from a library 
            where keys are the names of each geneset and values are capitalized, HGNC-style
            gene-names.
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
    if user_geneset is not None:
        geneset = user_geneset
    else:
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
    
### Add a special function to run gene set enrichment for OpenTargets GWAS data
def compute_geneset_enrichment_ot_gwas(mdata, gwas_data, prog_key='cNMF', data_key='rna', 
                                       organism='Human', library='GWAS', n_jobs=1, inplace=True,
                                       key_column='trait_efos', gene_column="gene_name", 
                                       **kwargs):
    """
    Perform gene set enrichment analysis using GWAS data.

    This function computes gene set enrichment using GWAS data. The GWAS data can be provided
    as a pandas DataFrame or a file path to a CSV file.

    Parameters:
    - mdata : MuData
        Mudata object containing anndata of program scores and cell-level metadata.
    - gwas_data : pandas DataFrame or str
        GWAS data provided as a pandas DataFrame or a file path to a CSV file.
    - prog_key : str, optional (default: 'cNMF')
        Cell program key name in the MuData object (e.g. cNMF)
    - data_key : str, optional (default: 'rna')
        Genomic data name in the MuData object (e.g. rna). This is where gene names are pulled from.
    - organism : str, optional (default: 'Human')
        Annotation of the species used for enrichment. Almost always Human for GWAS data.
    - library : str, optional (default: 'GWAS')
        How to name the user-defined gene set library
    - n_jobs : int, optional (default: 1)
        Number of threads to run processes on.
    - inplace : bool, optional (default: True)
        Update the mudata object inplace or return a copy.
    - key_column : str, optional (default: 'trait_efos')
        Name of the column in the GWAS data DataFrame to use as the key for the geneset dictionary.
    - gene_column : str, optional (default: 'gene_name')
        Name of the column in the GWAS data DataFrame containing the gene names.
    - **kwargs : dict
        Additional keyword arguments to pass to the underlying function.

    """
    # Read GWAS data from disc if provided as a path
    if isinstance(gwas_data, str):
        df = pd.read_csv(gwas_data)
    elif isinstance(gwas_data, pd.DataFrame):
        df = gwas_data
    else:
        raise ValueError("gwas_data must be either a pandas DataFrame or a file path to a CSV file.")
    
    # Create geneset dictionary from GWAS data
    gmt = create_geneset_dict(df, key_column=key_column, gene_column=gene_column)
    
    # Perform gene set enrichment using the created geneset dictionary
    compute_geneset_enrichment(mdata=mdata, prog_key=prog_key, data_key=data_key, 
                               organism=organism, library=library, 
                               database=None, n_jobs=n_jobs, inplace=inplace, user_geneset=gmt)
