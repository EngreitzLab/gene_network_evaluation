import os
import gin
import argparse

import mudata
import anndata
from tqdm.auto import tqdm

import numpy as np

# consensus NMF described -> https://github.com/dylkot/cNMF
def init_cnmf_obj(output_dir=None, name=None):
    
    from cnmf import cNMF

    # Compute cNMF and create prog anndata
    cnmf_obj = cNMF(output_dir=output_dir, name=name)

    return cnmf_obj

def run_cnmf_factorization(output_dir=None, name=None, counts_fn=None,
                           components=[7,8,9,10], n_iter=10, seed=14,
                           total_workers=-1, num_highvar_genes=2000, 
                           solver='halsvar', beta_loss='frobenius', use_gpu=False, 
                           mode='original'):

    if mode=='torch':

        # Compute cNMF and create prog anndata
        cnmf_obj = init_cnmf_obj(output_dir=output_dir, name=name)
        cnmf_obj.prepare(counts_fn=counts_fn, components=components, n_iter=n_iter, 
                        seed=seed, num_highvar_genes=num_highvar_genes, 
                        solver=solver, beta_loss=beta_loss,
                        total_workers=total_workers, use_gpu=use_gpu)
        cnmf_obj.factorize()
    
    elif mode=='original':

        # Compute cNMF and create prog anndata
        cnmf_obj = init_cnmf_obj(output_dir=output_dir, name=name)
        cnmf_obj.prepare(counts_fn=counts_fn, components=components, n_iter=n_iter, 
                        seed=seed, num_highvar_genes=num_highvar_genes, beta_loss=beta_loss)
        # FIXME: Doesn't seem to work in multithreading
        cnmf_obj.factorize(total_workers=total_workers)

    return cnmf_obj

def run_cnmf_consensus(cnmf_obj=None, output_dir=None, name=None, 
                       components=[7,8,9,10], density_thresholds=[0.01, 0.05, 2.0]):
    
    #TODO: Usually this step has to be rerun with a manually chosen density threshold.

    if cnmf_obj is None:
        cnmf_obj = init_cnmf_obj(output_dir=output_dir, name=name)

    for k in tqdm(components, desc='Running cNMF'):
        for thresh in density_thresholds:
            cnmf_obj.consensus(k=k, density_threshold=thresh, show_clustering=True)

@gin.configurable
def run_consensus_NMF_(K=10, output_dir=None, name=None, counts_fn=None,
                       components=[7,8,9,10], n_iter=10, seed=14,
                       total_workers=-1, density_thresholds=[0.01, 2.0],
                       num_highvar_genes=2000, beta_loss='frobenius',
                       output_all_k=True, output_all_thresh=True,
                       use_gpu=False, mode='original'):

    cnmf_obj = run_cnmf_factorization(output_dir=output_dir, name=name, counts_fn=counts_fn,
                                      components=components, n_iter=n_iter, seed=seed,
                                      total_workers=total_workers, num_highvar_genes=num_highvar_genes, 
                                      beta_loss=beta_loss, use_gpu=use_gpu, mode=mode)
    cnmf_obj.combine()
    cnmf_obj.k_selection_plot()

    # Plot & store for many 
    run_cnmf_consensus(cnmf_obj=cnmf_obj, 
                       components=components, 
                       density_thresholds=density_thresholds)
    
    return cnmf_obj, K, components, density_thresholds, output_all_k, output_all_thresh

def run_consensus_NMF(mdata, work_dir='./', scratch_dir=None,  
                      prog_key='consensus_NMF', data_key='rna', 
                      layer='X', config_path=None, n_jobs=1, 
                      inplace=True, use_gpu=False, mode='original'):

    """
    Perform gene program inference using consensus NMF.

    https://github.com/dylkot/cNMF

    ARGS:
        mdata : MuData
            mudata object containing anndata of data and cell-level metadata.
        work_dir: str
            path to directory to store outputs.
        scratch_dir: str
            path to directory to store intermediate outputs.
        prog_key: 
            index for the anndata object (mdata[prog_key]) in the mudata object.
        data_key: str
            index of the genomic data anndata object (mdata[data_key]) in the mudata object.
        layer: str (default: X)
            anndata layer (mdata[data_key].layers[layer]) where the data is stored.
        config_path: str
            path to gin configurable config file containing method specific parameters.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        use_gpu: bool (default: False)
            use GPU if available. Only works with mode='torch'.
        mode: str (default: 'original')
            Use cNMF_torch ('torch') or original ('original').
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy

    RETURNS:
        if not inplace:
            mdata
    
    """

    # Load method specific parameters
    try: gin.parse_config_file(config_path)
    except: raise ValueError('gin config file could not be found')

    if not inplace:
        mdata = mudata.MuData({data_key: mdata[data_key].copy()})

    # Create output directory for cNMF results
    if work_dir is not None:
        try: os.makedirs(work_dir, exist_ok=True)
        except: raise ValueError('Work directory location does not exist.')

    # Store temporary anndata
    if scratch_dir==None:
        scratch_dir=work_dir

    # Create temp anndata 
    if layer=='X':
        counts_path = os.path.join(scratch_dir, '{}_temp.h5ad'.format(data_key))
        temp_data = mdata[data_key].copy()
        temp_data.var_names_make_unique()
        temp_data.write(counts_path)
    else:
        temp_data = anndata.AnnData(data=mdata[data_key].layers[layer], 
                                    obs=mdata[data_key].obs,
                                    var=mdata[data_key].var)
        counts_path = os.path.join(scratch_dir, '{}_{}_temp.h5ad'.format(data_key, layer))
        temp_data.var_names_make_unique()
        temp_data.write(counts_path)
  
    # Compute cNMF and create prog anndata
    cnmf_object, K, components, density_thresholds, \
    output_all_k, output_all_thresh = \
    run_consensus_NMF_(output_dir=work_dir, 
                       counts_fn=counts_path,
                       total_workers=n_jobs,
                       use_gpu=use_gpu,
                       mode=mode)

    # Create new anndata object
    usage, spectra_scores, spectra_tpm, top_genes = \
    cnmf_obj.load_results(K=K, density_threshold=min(density_thresholds))

    mdata.mod[prog_key] = anndata.AnnData(X=usage, obs=mdata[data_key].obs)
    mdata[prog_key].varm['loadings'] = spectra_tpm
    mdata[prog_key].varm['loadings_zscore'] = spectra_scores
    mdata[prog_key].uns['loadings_genes'] = top_genes

    # Store outputs in mdata
    if not output_all_k:
        components = [K]
    if not output_all_thresh:
        density_thresholds = [min(density_thresholds), 
                              max(density_thresholds)]
    adatas = {}
    for k in tqdm(components, desc='Storing output'):
        for thresh in density_thresholds:
            usage, spectra_scores, spectra_tpm, top_genes = \
            cnmf_obj.load_results(K=k, density_threshold=thresh)

            adata_ = anndata.AnnData(data=usage, obs=mdata[data_key].obs)
            adata_.varm['loadings'] = spectra_tpm
            adata_.varm['loadings_zscore'] = spectra_scores
            adata_.uns['loadings_genes'] = top_genes

            adatas[prog_key+'_{}_{}'.format(k, thresh)] = adata_
    
    adatas[data_key] = mdata[data_key]
    adatas[prog_key] = mdata[prog_key]

    mdata = mudata.MuData(adatas)

    if not inplace: return mdata

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('--work_dir', default='./', type=str)
    parser.add_argument('--scratch_dir', default='./', type=str)
    parser.add_argument('-pk', '--prog_key', default='consensus_NMF', typ=str) 
    parser.add_argument('-dk', '--data_key', default='rna', typ=str) 
    parser.add_argument('--layer', default='X', type=str)
    parser.add_argument('--config_path', default='./consensus_NMF_config.gin', type=str)
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--use_gpu', default=False, store_action=True)
    parser.add_argument('--mode', default='original', type=str, choices=['torch', 'original'])
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    run_consensus_NMF(mdata, work_dir=args.work_dir, scratch_dir=args.scratch_dir, 
                      prog_key=args.prog_key, data_key=args.data_key, layer=args.layer, 
                      config_path=args.config_path, n_jobs=args.n_jobs, inplace=args.output,
                      use_gpu=args.use_gpu, mode=args.mode)