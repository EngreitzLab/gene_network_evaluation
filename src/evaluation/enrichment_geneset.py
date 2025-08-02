import argparse
import logging
import os
from typing import Dict, List, Optional, Tuple, Union, Literal

import gseapy as gp
import mudata
import numpy as np
import pandas as pd
from anndata import AnnData
from gseapy import Biomart, Msigdb
from tqdm.auto import tqdm

logging.basicConfig(level=logging.INFO)


def create_geneset_dict(
    dataframe: pd.DataFrame, 
    key_column='trait_efos', 
    gene_column='gene_name'
):
    """Create custom geneset"""

    geneset_dict = {}
    for _, row in dataframe.iterrows():
        key = row[key_column]
        gene = row[gene_column]
        geneset_dict.setdefault(key, []).append(gene)
    return geneset_dict

def get_idconversion(var_names, organism='human', chunk_size=100):
    """Match gene IDs with database, chunking large queries to avoid HTML errors."""
    # Ensure var_names is a simple list
    var_names = list(var_names)

    bm = Biomart()
    # Default to uppercase fallback
    gene_names = [v.upper() for v in var_names]

    # Detect Ensembl IDs by looking at the first element safely
    is_ensembl = (len(var_names) > 0 and var_names[0].lower().startswith('ens'))
    if is_ensembl:
        dataset = (
            'hsapiens_gene_ensembl'
            if organism=='human'
            else 'mmusculus_gene_ensembl'
        )

        lookup = {}
        for i in range(0, len(var_names), chunk_size):
            chunk = var_names[i:i+chunk_size]
            try:
                res = bm.query(
                    dataset=dataset,
                    attributes=['ensembl_gene_id', 'external_gene_name'],
                    filters={'ensembl_gene_id': chunk}
                )
                if not hasattr(res, 'columns'):
                    raise RuntimeError('BioMart returned HTML or invalid format')
                df = res[['ensembl_gene_id', 'external_gene_name']].astype(str)
                df['external_gene_name'] = df['external_gene_name'].str.upper()
                lookup.update(dict(zip(df['ensembl_gene_id'], df['external_gene_name'])))
            except Exception as e:
                logging.info(f"BioMart chunk lookup failed for IDs {i}-{i+len(chunk)} ({e}); skipping these.")

        gene_names = [lookup.get(v, v.upper()) for v in var_names]

    elif any(':ens' in v.lower() for v in var_names):
        gene_names = [v.split(':', 1)[0].upper() for v in var_names]

    return gene_names
        
def get_program_gene_loadings(mdata, prog_key='prog', prog_nam=None, data_key='rna', organism='human'):
    """Get gene loadings for each program in the mudata object."""

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

    # with open('var_names.txt', 'w') as fil:
    #     for nam in mdata[data_key].var_names:
    #         fil.write(nam+'\n')

    return loadings


def get_geneset(organism='human', library='h.all', database='msigdb', min_size: int = 0, max_size: int = 2000):
    """Download gene set from MsigDB or Enrichr."""
    if database == 'msigdb':
        msig = Msigdb()
        dbver = '2023.2.Hs' if organism == 'human' else '2023.1.Mm'
        gmt = msig.get_gmt(category=library, dbver=dbver)
        if gmt is None:
            raise ValueError('Library does not exist')
    elif database == 'enrichr':
        gmt = gp.get_library(name=library, organism=organism.capitalize(), min_size=min_size, max_size=max_size)
    return gmt


def perform_prerank(
    loadings: pd.DataFrame,
    geneset: Union[List[str], str, Dict[str, str]],
    n_jobs: int = 1,
    low_cutoff: float = -np.inf,
    n_top: int = None,
    **kwargs
) -> pd.DataFrame:
    """Run GSEA prerank on each gene program in the loadings matrix.
    
    Parameters
    ----------
    loadings : pd.DataFrame
        DataFrame of feature loadings for each program.
    geneset : str
        Name of the set to run GSEA on.
    n_jobs : int
        Number of parallel jobs to run.
    low_cutoff : float
        Remove features with loadings at or below this value.
    n_top : int
        Take the top n features with the highest loadings.
        
    Returns
    -------
    pd.DataFrame
        DataFrame of GSEA results sorted by program name and FDR q-value. Includes the following columns:
        - program_name: name of the program
        - Term: gene set name
        - ES: enrichment score
        - NES: normalized enrichment score
        - NOM p-val: nominal p-value (from the null distribution of the gene set)
        - FDR q-val: adjusted False Discovery Rate
        - FWER p-val: Family wise error rate p-values
        - Gene %: percent of gene set before running enrichment peak (ES)
        - Lead_genes: leading edge genes (gene hits before running enrichment peak)
        - tag_before: number of genes in gene set
        - tag_after: number of genes matched to the data
    """

    # Run GSEA prerank for each column of loadings (each cell program)
    pre_res = pd.DataFrame()
    for i in tqdm(loadings.columns, desc='Running GSEA', unit='programs'):

        # Filter out low loadings
        temp_loadings = loadings[i][(loadings[i] > low_cutoff)]

        # Take top n features if specified
        if n_top is not None:
            temp_loadings = temp_loadings.sort_values(ascending=False).head(n_top)
            if len(temp_loadings) < n_top:
                logging.warning(f"Program {i} has less than {n_top} features after filtering. Only {len(temp_loadings)} features will be used.")

        # Run GSEA prerank
        temp_res = gp.prerank(
            rnk=temp_loadings, 
            gene_sets=geneset, 
            threads=n_jobs, 
            **kwargs
        ).res2d

        # Post-process results
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


def perform_fisher_enrich(
    loadings, 
    geneset, 
    n_top=500,
    **kwargs
):
    """Run Fisher enrichment on each gene program in the loadings matrix.

    Parameters
    ----------
    loadings : pd.DataFrame
        DataFrame of feature loadings for each program.
    geneset : dict
        Dictionary of gene sets.
    n_top : int
        Number of top features to take.
    
    Returns
    -------
    ['Gene_set', 'Term', 'P-value', 'Adjusted P-value', 'Odds Ratio',
       'Combined Score', 'Genes', 'program_name', 'overlap_numerator',
       'overlap_denominator'],
    pd.DataFrame
        DataFrame of Fisher enrichment results sorted by program name and adjusted p-value. Includes the following columns:
        - program_name: name of the program
        - Term: gene set name
        - P-value: Fisher's exact test p-value
        - Adjusted P-value: adjusted p-value
        - Odds Ratio: odds ratio
        - Combined Score: combined score
        - Genes: genes in the gene set
        - overlap_numerator: number of overlapping genes
        - overlap_denominator: number of genes in the gene set
    TODO
    ----
    - Parallelize across programs
    """

    # Find the intersection of genes present in the mudata object and in the library
    background_genes = set(value for sublist in geneset.values() for value in sublist)
    
    enr_res = pd.DataFrame()
    for i in tqdm(loadings.columns, desc='Running Fisher enrichment', unit='programs'):
        gene_list = list(loadings[i].sort_values(ascending=False).head(n_top).index)
        temp_res = gp.enrich(
            gene_list=gene_list,
            gene_sets=geneset, 
            background=background_genes
        ).res2d
        temp_res["program_name"] = i
        enr_res = pd.concat([enr_res, temp_res], ignore_index=True)
    enr_res['overlap_numerator'] = enr_res['Overlap'].apply(lambda x: int(x.split('/')[0]))
    enr_res['overlap_denominator'] = enr_res['Overlap'].apply(lambda x: int(x.split('/')[1]))
    enr_res.drop(columns=['Overlap', 'Gene_set'], inplace=True)
    enr_res = enr_res[["program_name"] + [col for col in enr_res.columns if col != "program_name"]]
    enr_res = enr_res.sort_values(['program_name', 'Adjusted P-value'])
    
    return enr_res


def insert_enrichment(
    mdata: mudata.MuData,
    df: pd.DataFrame,
    library="GSEA", 
    prog_key="prog",
    geneset_index="Term", 
    program_index="program_name",
    varmap_name_prefix="gsea_varmap"
) -> None:
    """Insert geneset enrichment into mudata
    
    Parameters
    ----------
    mdata : mudata.MuData
        MuData object.
    df : pd.DataFrame
        DataFrame of geneset enrichment results.
    library : str
        name of the library used for enrichment.
    prog_key : str
        key for the program in the mdata object.
    geneset_index : str
        index for the gene set in the DataFrame.
    program_index : str
        index for the program in the DataFrame.
    varmap_name_prefix : str
        prefix for the varmap name.
    
    Returns
    -------
    None
    """
    
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


def compute_geneset_enrichment(
    mdata: Union[str, mudata.MuData],
    prog_key: str = 'prog',
    data_key: str = 'rna',
    prog_name: Optional[str] = None,
    method: Literal['gsea', 'fisher'] = 'gsea',
    organism: Literal['human', 'mouse'] = 'human',
    library: str = 'Reactome_2022', 
    database: Literal['msigdb', 'enrichr'] = 'enrichr',
    user_geneset: Optional[Dict[str, List[str]]] = None,
    min_size: int = 0,
    max_size: int = 2000,
    low_cutoff: float = -np.inf,
    n_top: int = 2000,
    n_jobs: int = 1,
    inplace: bool = True,
    **kwargs
) -> Optional[pd.DataFrame]:

    """
    Wrapper function to compute gene set enrichment for each program in the MuData object.

    Parameters
    ----------
    mdata : Union[str, mudata.MuData]
        Path to the MuData object or the MuData object itself.
    prog_key : str
        index for the anndata object (mdata[prog_key]) in the mudata object.
    data_key : str
        index of the genomic data anndata object (mdata[data_key]) in the mudata object.
    prog_name : str (default: None)
        Compute enrichment for a particular program.
    method : {'gsea', 'fisher'} (default: 'gsea')
        Run GSEA or Fisher exact test gene set enrichment.
    organism : {'human', 'mouse'} (default: 'human')
        species to which the sequencing data was aligned to.
    library : str (default: Reactome_2022)
        gene-set library to use for computing enrichment.
        MsigDB libraries: https://www.gsea-msigdb.org/gsea/msigdb
        Enrichr libraries: https://maayanlab.cloud/Enrichr/#libraries
    database : {'msigdb', 'enrichr'} (default: 'enrichr')
        database of gene-set libraries to use. Should match the library.
    user_geneset : dict
        user-defined gene set to use for enrichment. If provided, library and database are ignored.
    min_size: int (default: 0)
        minimum size of gene sets to consider.
    max_size: int (default: 2000)
        maximum size of gene sets to consider.
    low_cutoff : float (default: -np.inf)
        Remove features with loadings at or below this value.
    n_top : int
        Take the top n features with the highest loadings.
    n_jobs : int (default: 1)
        number of threads to run processes on.
    inplace : bool (default: True)
        whether to insert the results back into the mudata object.
    
    Returns
    -------
    if not inplace:
        return pre_res
    else:
        inserts the results back into the mudata object
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
    
    # get the geneset
    if user_geneset is not None:
        geneset = user_geneset
    else:
        geneset = get_geneset(
            organism=organism,
            library=library,
            database=database,
            max_size=max_size,
            min_size=min_size
        )
     
    # get the gene loadings for each program
    loadings = get_program_gene_loadings(
        mdata, 
        prog_key=prog_key, 
        prog_nam=prog_name,
        data_key=data_key, 
        organism=organism,
    )
    
    # run enrichment
    if method == "gsea":
        pre_res = perform_prerank(
            loadings=loadings, 
            geneset=geneset,
            n_jobs=n_jobs,
            low_cutoff=low_cutoff,
            n_top=n_top,
            **kwargs
        )
    elif method == "fisher":
        pre_res = perform_fisher_enrich(
            loadings=loadings,
            geneset=geneset,
            n_top=n_top,
            **kwargs
        )
        
    # insert results back into the mudata object if inplace
    if inplace:
        mdata=insert_enrichment(
            mdata, df=pre_res, 
            library=library, 
            prog_key=prog_key,
            geneset_index="Term", 
            program_index="program_name",
            varmap_name_prefix="gsea_varmap"
        )
    else:
        return(pre_res)
    