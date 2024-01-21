import os
import gin
import argparse

import anndata
from tqdm.auto import tqdm

# consensus NMF described -> https://github.com/dylkot/cNMF
@gin.configurable
def run_consensus_NMF_(K=10, output_dir=None, name=None, counts_fn=None,
                       components=[7,8,9,10], n_iter=100, seed=14,
                       total_workers=-1, density_thresholds=[0.01, 2],
                       num_highvar_genes=2000, beta_loss='frobenius',
                       output_all_k=True, output_all_thresh=True):
  
    from cnmf import cNMF

    # Compute cNMF and create prog anndata
    cnmf_obj = cNMF(output_dir=output_dir, name=name)
    cnmf_obj.prepare(counts_fn=counts_fn, components=components, 
                     n_iter=n_iter, seed=seed, num_highvar_genes=num_highvar_genes)
    cnmf_obj.factorize(total_workers=total_workers)
    cnmf_obj.combine()
    cnmf_obj.k_selection_plot()
    
    # Plot & store for many 
    for k in tqdm(components, desc='Running cNMF'):
        for thresh in density_thresholds:
            cnmf_obj.consensus(k=k, density_threshold=thresh, show_clustering=True)
    
    return cnmf_obj, K, components, density_thresholds, output_all_k, output_all_thresh

def run_consensus_NMF(mdata, work_dir='./', scratch_dir=None, n_jobs=-1, 
                      prog_key='consensus_NMF', data_key='rna', layer='X', 
                      config_path=None, inplace=True):
    
    # Load method specific parameters
    try: gin.parse_config_file(config_path)
    except: raise ValueError('gin config file could not be found')

    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata

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
                       total_workers=n_jobs)

    # Create new anndata object
    usage, spectra_scores, spectra_tpm, top_genes = \
    cnmf_obj.load_results(K=K, density_threshold=min(density_thresholds))

    mdata[prog_key] = anndata.AnnData(data=usage, 
                                      obs=mdata[data_key].obs)
    mdata[prog_key].varm['loadings'] = spectra_tpm
    mdata[prog_key].varm['loadings_zscore'] = spectra_scores
    mdata[prog_key].uns['loadings_genes'] = top_genes

    # Store outputs in mdata
    if not output_all_k:
        components = [K]
    if not output_all_thresh:
        density_thresholds = [min(density_thresholds), 
                              max(density_thresholds)]
    for k in tqdm(components, desc='Storing output'):
        for thresh in density_thresholds:
            usage, spectra_scores, spectra_tpm, top_genes = \
            cnmf_obj.load_results(K=k, density_threshold=thresh)

            mdata[prog_key+'_{}_{}'.format(k, thresh)] = \
            anndata.AnnData(data=usage, obs=mdata[data_key].obs)

            mdata[prog_key+'_{}_{}'.format(k, thresh)].varm['loadings'] = \
            spectra_tpm
            mdata[prog_key+'_{}_{}'.format(k, thresh)].varm['loadings_zscore'] = \
            spectra_scores
            mdata[prog_key+'_{}_{}'.format(k, thresh)].uns['loadings_genes'] = \
            top_genes

    if not inplace: return mdata[prog_key]

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('-pk', '--prog_key', default='consensus_NMF', typ=str) 
    parser.add_argument('-dk', '--data_key', default='rna', typ=str) 
    parser.add_argument('--layer', default='X', type=str)
    parser.add_argument('--work_dir', default='./', type=str)
    parser.add_argument('--config_path', default='./consensus_NMF_config.gin', type=str)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    import mudata
    mdata = mudata.read(args.mudataObj_path)

    run_consensus_NMF(mdata, work_dir=args.work_dir, n_jobs=args.n_jobs, layer=args.layer,
                      prog_key=args.prog_key, data_key=args.data_key, inplace=args.output,
                      config_path=args.config_path)