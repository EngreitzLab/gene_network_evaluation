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
logging.basicConfig(level=logging.INFO)

# Create custom geneset
def create_geneset_dict(dataframe: pd.DataFrame, 
                        key_column='trait_efos', 
                        gene_column='gene_name'):

    geneset_dict = {}
    for _, row in dataframe.iterrows():
        key = row[key_column]
        gene = row[gene_column]
        geneset_dict.setdefault(key, []).append(gene)
    return geneset_dict

# Match gene IDs with database
def get_idconversion(var_names, organism='human'):

    bm = Biomart()
    gene_names = []
    for i in range(min(10, len(var_names))):
        if var_names[i].lower().startswith('ens'):
            queries = {'ensembl_gene_id': list(var_names)}

            if organism=='human':
                id_dataset='hsapiens_gene_ensembl'
            elif organism=='mouse':
                id_dataset='mmusculus_gene_ensembl'

            results = bm.query(dataset=id_dataset,
                               attributes=['ensembl_gene_id', 'external_gene_name'],
                               filters=queries)

            if type(results) is str:
                raise RuntimeError('Gene name query request did not suceed. Consider converting ENS IDs to gene names manually')
            gene_names = results['external_gene_name'].apply(lambda x: x.upper()).values
            break
        elif ':ens' in var_names[i].lower():
            gene_names = [name.split(':')[0].upper() for name in var_names]
            break
        else:
            gene_names = [name.upper() for name in var_names]
    return gene_names

# Extract program loadings
def get_program_gene_loadings(mdata, prog_key='prog', prog_nam=None, data_key='rna', organism='human'):

    if 'var_names' in mdata[prog_key].uns.keys():
        gene_names = get_idconversion(mdata[prog_key].uns['var_names'], organism=organism)
    else:
        assert mdata[prog_key].varm['loadings'].shape[1] == mdata[data_key].var.shape[0]
        gene_names = get_idconversion(mdata[data_key].var_names, organism=organism)

    if prog_nam:
        loadings = pd.DataFrame(data=mdata[prog_key][:, prog_nam].varm['loadings'].flatten(), index=gene_names)
        loadings.columns = [prog_nam]
    else:
        loadings = pd.DataFrame(data=mdata[prog_key].varm['loadings'], index=mdata[prog_key].var.index).T
        loadings["gene_names"] = gene_names
        loadings.set_index("gene_names", inplace=True)

    with open('var_names.txt', 'w') as fil:
        for nam in mdata[data_key].var_names:
            fil.write(nam+'\n')

    return loadings

# Download geneset
def get_geneset(organism='human', library='h.all', database='msigdb'):

    if database == 'msigdb':
        msig = Msigdb()
        dbver = '2023.2.Hs' if organism == 'human' else '2023.1.Mm'
        gmt = msig.get_gmt(category=library, dbver=dbver)
        if gmt is None:
            raise ValueError('Library does not exist')
    elif database == 'enrichr':
        gmt = gp.get_library(name=library, organism=organism.capitalize())
    return gmt

# Run GSEAS
def perform_prerank(loadings, geneset, n_jobs=1, **kwargs):

    # Defaults - pass kwargs if you want this done differently
    #threads=4, min_size=5, max_size=1000, permutation_num=1000, 
    #outdir=None, seed=0, verbose=True

    # Run GSEA prerank for each column of loadings (each cell program)
    pre_res = pd.DataFrame()
    for i in tqdm(loadings.columns, desc='Running GSEA', unit='programs'):
        temp_res = gp.prerank(rnk=loadings[i], gene_sets=geneset, threads=n_jobs, **kwargs).res2d

        temp_res['Gene %'] = temp_res['Gene %'].apply(lambda x: float(x[:-1]))
        temp_res['tag_before'] = temp_res['Tag %'].apply(lambda x: int(x.split('/')[0]))
        temp_res['tag_after'] = temp_res['Tag %'].apply(lambda x: int(x.split('/')[1]))
        temp_res.drop(columns=['Tag %'], inplace=True)
    
        if 'Name' in temp_res.columns and temp_res['Name'][0] == "prerank":
            temp_res['Name'] = i
    
        temp_res.rename(columns={'Name': 'program_name'}, inplace=True)
        temp_res = temp_res.sort_values(['program_name', 'FDR q-val'])
        pre_res = pd.concat([pre_res, temp_res], ignore_index=True)
    
    return pre_res

# TODO: ssGSEA using loading matrix -> get topic wise enrichment scores
def perform_ssGSEA():
    raise NotImplementedError()

def perform_fisher_enrich(loadings, geneset, loading_rank_thresh=500, **kwargs):
    
    # Find the intersection of genes present in the mudata object and in the library
    background_genes = set(value for sublist in geneset.values() for value in sublist)
    
    enr_res = pd.DataFrame()
    # TODO: Parallelize
    for i in tqdm(loadings.columns, desc='Running Fisher enrichment', unit='programs'):
        gene_list = list(loadings[i].sort_values(ascending=False).head(loading_rank_thresh).index)
        temp_res = gp.enrich(gene_list=gene_list,
                             gene_sets=geneset, background=background_genes).res2d
        temp_res["program_name"] = i
        enr_res = pd.concat([enr_res, temp_res], ignore_index=True)
    enr_res['overlap_numerator'] = enr_res['Overlap'].apply(lambda x: int(x.split('/')[0]))
    enr_res['overlap_denominator'] = enr_res['Overlap'].apply(lambda x: int(x.split('/')[1]))
    enr_res.drop(columns=['Overlap'], inplace=True)
    
    return enr_res

# Insert geneset enrichment into mudata
def insert_enrichment(mdata, df, library="GSEA", prog_key="prog",
                      geneset_index="Term", program_index="program_name",
                      varmap_name_prefix="gsea_varmap"):
    
    # Create a mudata key to column name mapping dictionary
    mudata_keys_dict = {}
    for col in df.columns:
        if col not in [geneset_index, program_index]:
            key = f"{col}_{library}"
            key = key.replace(' ', '_').replace('%', 'percent')
            mudata_keys_dict[key] = col

    # Insert the values from the dataframe into the array for each key
    for key, colname in mudata_keys_dict.items():
        # Create an empty dataframe with the right dimensions
        all_progs_df = pd.DataFrame(index=df[geneset_index].unique(), 
                                    columns=mdata[prog_key].var.index)
        
        # Pivot the dataframe for gene sets and programs
        pivot_df = df[[geneset_index, program_index, colname]].pivot(index=geneset_index, 
                                                                     columns=program_index, 
                                                                     values=colname)
        
        # Update the empty dataframe with new values
        all_progs_df[pivot_df.columns] = pivot_df
        
        # Convert dataframe to a numpy array
        all_progs_array = all_progs_df.T.to_numpy()
        
        # Add the array into the MuData object
        mdata[prog_key].varm[key] = all_progs_array
        
    # Add the varmap to the mudata object
    varmap_name = f"{varmap_name_prefix}_{library}"
    mdata[prog_key].uns[varmap_name] = all_progs_df.index

def compute_geneset_enrichment(mdata, prog_key='prog', data_key='rna', prog_nam=None,
                               organism='human', library='Reactome_2022', method="gsea",
                               database='enrichr', n_jobs=1, inplace=False, user_geneset=None, 
                               loading_rank_thresh=500, **kwargs):

    """
    Perform GSEA using loadings are preranked gene lists.

    ARGS
        mdata : MuData
            mudata object containing anndata of program scores and cell-level metadata.
        prog_key: 
            index for the anndata object (mdata[prog_key]) in the mudata object.
        data_key: str
            index of the genomic data anndata object (mdata[data_key]) in the mudata object.
        prog_nam: 
            Compute enrichment for a particular program.
        organism: {'human', 'mouse'} (default: 'human')
            species to which the sequencing data was aligned to.
        library: str (default: Reactome_2022)
            gene-set library to use for computing enrichment.
            MsigDB https://www.gsea-msigdb.org/gsea/msigdb
            Enrichr https://maayanlab.cloud/Enrichr/#libraries
        database: {'msigdb', 'enrichr'} (default: 'enrichr')
            database of gene-set libraries to use.
        method: {'gsea', 'fisher'}
            Run gseas or fisher gene set enrichment.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy
        user_geneset: gmt
            User specified gmt file for GSEA.
        loading_rank_thresh: int (default: 500)
            Threshold loadings to compute fishers exact test.
       
    RETURNS 
        if not inplace:
            return pre_res       

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
    
    #get the geneset
    if user_geneset is not None:
        geneset = user_geneset
    else:
        geneset = get_geneset(organism, library, database)
     
    #get the gene loadings for each program
    loadings = get_program_gene_loadings(mdata, prog_key=prog_key, prog_nam=prog_nam, 
                                         data_key=data_key, organism=organism)
    
    #run enrichment
    if method == "gsea":
        pre_res = perform_prerank(loadings=loadings, geneset=geneset, n_jobs=n_jobs)
    elif method == "fisher":
        pre_res = perform_fisher_enrich(loadings=loadings, geneset=geneset, loading_rank_thresh=loading_rank_thresh)
        
    #return the result depending on whether or not we want to do it inplace   
    if inplace:
        mdata=insert_enrichment(mdata, df=pre_res, library=library, prog_key=prog_key,
                                geneset_index="Term", program_index="program_name",
                                varmap_name_prefix="gsea_varmap")
    if not inplace:
        return(pre_res)

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

  