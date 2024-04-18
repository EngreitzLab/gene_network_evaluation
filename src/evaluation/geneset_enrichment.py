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

def create_geneset_dict(df, key_column='trait_efos', gene_column='gene_name'):
    geneset_dict = {}
    for _, row in df.iterrows():
        key = row[key_column]
        gene = row[gene_column]
        geneset_dict.setdefault(key, []).append(gene)
    return geneset_dict


def get_idconversion(var_names):
    bm = Biomart()
    gene_names = []
    for i in range(min(10, len(var_names))):
        if var_names[i].lower().startswith('ens'):
            queries = {'ensembl_gene_id': list(var_names)}
            gene_names = queries['external_gene_name'].apply(lambda x: x.upper()).values
            break
        elif ':ens' in var_names[i].lower():
            gene_names = [name.split(':')[0].upper() for name in var_names]
            break
        else:
            gene_names = [name.upper() for name in var_names]
    return gene_names


def get_program_gene_loadings(mdata, prog_key='cNMF', prog_nam=None, data_key='rna'):
    if 'var_names' in mdata[prog_key].uns.keys():
        gene_names = get_idconversion(mdata[prog_key].uns['var_names'])
    else:
        assert mdata[prog_key].varm['loadings'].shape[1] == mdata[data_key].var.shape[0]
        gene_names = get_idconversion(mdata[data_key].var_names)

    if prog_nam:
        loadings = pd.DataFrame(data=mdata[prog_key][:, prog_nam].varm['loadings'].flatten(), index=gene_names)
        loadings.columns = [prog_nam]
    else:
        loadings = pd.DataFrame(data=mdata[prog_key].varm['loadings'], index=mdata[prog_key].var.index).T
        loadings["gene_names"] = gene_names
        loadings.set_index("gene_names", inplace=True)
    
    return loadings

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


def perform_prerank(loadings, geneset, n_jobs=1, **kwargs):

    
    # Run GSEA prerank for each column of loadings (each cell program)
    pre_res = gp.prerank(rnk=loadings, gene_sets=geneset, threads=n_jobs, **kwargs).res2d
    pre_res['Gene %'] = pre_res['Gene %'].apply(lambda x: float(x[:-1]))
    pre_res['tag_before'] = pre_res['Tag %'].apply(lambda x: int(x.split('/')[0]))
    pre_res['tag_after'] = pre_res['Tag %'].apply(lambda x: int(x.split('/')[1]))
    pre_res.drop(columns=['Tag %'], inplace=True)
    
    if 'Name' in pre_res.columns and pre_res['Name'][0] == "prerank":
        pre_res['Name'] = loadings.columns[0]
    
    pre_res.rename(columns={'Name': 'program_name'}, inplace=True)
    pre_res = pre_res.sort_values(['program_name', 'FDR q-val'])
    
    return pre_res

def perform_fisher_enrich(loadings, geneset, loading_rank_thresh=500, **kwargs):
    
    #find the intersection of genes present in the mudata object and in the library
    background_genes = set(value for sublist in geneset.values() for value in sublist)
    
    enr_res = pd.DataFrame()
    for i in loadings.columns:
        gene_list = list(loadings[i].sort_values(ascending=False).head(loading_rank_thresh).index)
        temp_res = gp.enrich(gene_list=gene_list,
                             gene_sets=geneset, background=background_genes).res2d
        temp_res["program_name"] = i
        enr_res = pd.concat([enr_res, temp_res], ignore_index=True)
    enr_res['overlap_numerator'] = enr_res['Overlap'].apply(lambda x: int(x.split('/')[0]))
    enr_res['overlap_denominator'] = enr_res['Overlap'].apply(lambda x: int(x.split('/')[1]))
    enr_res.drop(columns=['Overlap'], inplace=True)
    
    return enr_res

def insert_df_into_mudata(mdata, df, library="GSEA", prog_key="cNMF",
                                               geneset_index="Term", program_index="program_name",
                                               varmap_name_prefix="gsea_varmap"):
    
    # Create a mudata key to column name mapping dictionary
    mudata_keys_dict = {}
    for col in df.columns:
        if col not in [geneset_index, program_index]:
            key = f"{col}_{library}"
            key = key.replace(' ', '_').replace('%', 'percent')
            mudata_keys_dict[key] = col

    #cp the mudata object
    mdata2=mdata

    # Insert the values from the dataframe into the array for each key
    for key, colname in mudata_keys_dict.items():
        # Create an empty dataframe with the right dimensions
        all_progs_df = pd.DataFrame(index=df[geneset_index].unique(), columns=mdata2[prog_key].var.index)
        
        # Pivot the dataframe for gene sets and programs
        pivot_df = df[[geneset_index, program_index, colname]].pivot(index=geneset_index, columns=program_index, values=colname)
        
        # Update the empty dataframe with new values
        all_progs_df[pivot_df.columns] = pivot_df
        
        # Convert dataframe to a numpy array
        all_progs_array = all_progs_df.T.to_numpy()
        
        # Add the array into the MuData object
        mdata2[prog_key].varm[key] = all_progs_array
        
    # Add the varmap to the mudata object
    varmap_name = f"{varmap_name_prefix}_{library}"
    mdata2[prog_key].uns[varmap_name] = all_progs_df.index
        
    return(mdata2)


def compute_geneset_enrichment(mdata, prog_key='cNMF', data_key='rna', organism='human', library='Reactome_2022',
                               method="gsea", prog_nam=None, database='enrichr', n_jobs=1, inplace=False,
                               user_geneset=None, loading_rank_thresh=500, output_file=None, **kwargs):

    #read in mudata if it is provided as a path
    if isinstance(mdata, str):
        mdata = mudata.read(mdata)
    
    #get the geneset
    if user_geneset is not None:
        geneset = user_geneset
    else:
        geneset = get_geneset(organism, library, database)
     
    #get the gene loadings for each program
    loadings = get_program_gene_loadings(mdata, prog_key=prog_key, prog_nam=prog_nam, data_key=data_key)
    
    #run enrichment
    if method == "gsea":
        pre_res = perform_prerank(loadings=loadings, geneset=geneset, n_jobs=n_jobs)
    elif method == "fisher":
        pre_res = perform_fisher_enrich(loadings=loadings, geneset=geneset, loading_rank_thresh=500)
        
    #return the result depending on whether or not we want to do it inplace
    if output_file:
        pre_res.to_csv(output_file, index=False)
    
    if inplace:
        mdata=insert_df_into_mudata(mdata, df=pre_res, library=library, prog_key=prog_key,
                                               geneset_index="Term", program_index="program_name",
                                               varmap_name_prefix="gsea_varmap")
    if not inplace and not output_file:
        return(pre_res)


def compute_geneset_enrichment_ot_gwas(mdata, gwas_data, prog_key='cNMF', prog_nam=None, data_key='rna',
                                       library='OT_GWAS', n_jobs=1, inplace=False, key_column='trait_efos',
                                       gene_column='gene_name', method='fisher', output_file=None, **kwargs):
    #read in gwas data
    if isinstance(gwas_data, str):
        df = pd.read_csv(gwas_data, compression='gzip', low_memory=False)
    elif isinstance(gwas_data, pd.DataFrame):
        df = gwas_data
    else:
        raise ValueError("gwas_data must be either a pandas DataFrame or a file path to a CSV file.")
    
    gmt = create_geneset_dict(df, key_column=key_column, gene_column=gene_column)

    return(compute_geneset_enrichment(mdata=mdata, prog_key=prog_key, data_key=data_key, library=library,
                               database=None, n_jobs=n_jobs, inplace=inplace, user_geneset=gmt, prog_nam=prog_nam, method=method,
                                     output_file=output_file))
