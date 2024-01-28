import os
import gin
import argparse

import mudata
import anndata
from tqdm.auto import tqdm

import logging
logging.basicConfig(level = logging.INFO)

# https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html
@gin.configurable
def run_linear_SCVI_(adata, layer='X', batch_key=None, labels_key=None, 
                     devices=-1, n_latent=10, max_epochs=250):
  
    import scvi

    if batch_key is not None and labels_key is not None:
        scvi.model.LinearSCVI.setup_anndata(adata, layer=layer,
                                            batch_key=batch_key,
                                            labels_key=labels_key)
    elif batch_key is not None:
        scvi.model.LinearSCVI.setup_anndata(adata, layer=layer,
                                            batch_key=batch_key)
    elif labels_key is not None:
        scvi.model.LinearSCVI.setup_anndata(adata, layer=layer,
                                            labels_key=labels_key)
    else:
        scvi.model.LinearSCVI.setup_anndata(adata, layer=layer)
                        
    model = scvi.model.LinearSCVI(adata, n_latent=n_latent)

    model.train(max_epochs=max_epochs,
                devices=devices, 
                plan_kwargs={"lr": 5e-3}, 
                check_val_every_n_epoch=10)

    train_elbo = model.history["elbo_train"][1:]
    test_elbo = model.history["elbo_validation"]

    ax = train_elbo.plot()
    test_elbo.plot(ax=ax)

    Z_hat = model.get_latent_representation()
    loadings = model.get_loadings()

    return model, ax, Z_hat, loadings

def run_linear_SCVI(mdata, n_jobs=-1, prog_key='linear_SCVI', data_key='rna',  
                    batch_key=None, labels_key=None, layer='X', config_path=None, 
                    inplace=True):
    
    # Load method specific parameters
    try: gin.parse_config_file(config_path)
    except: raise ValueError('gin config file could not be found')

    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata

    if layer=='X':
        layer_=None
    else:
        layer_=layer

    # Compute cNMF and create prog anndata
    # FIXME: Multi-processing crashes on HPC
    if n_jobs==-1:
        n_jobs=os.cpu_count()
    model, ax, Z_hat, loadings = run_linear_SCVI_(mdata[data_key], layer=layer_,
                                                  batch_key=batch_key, labels_key=labels_key,
                                                  devices=n_jobs)

    # Create new anndata object
    mdata.mod[prog_key]  = anndata.AnnData(X=Z_hat, obs=mdata[data_key].obs)
    mdata[prog_key].varm['loadings'] = loadings.values.T
    mdata[prog_key].uns['loadings_genes'] = loadings.index.values

    for col in mdata[data_key].obs.columns:
        if '_scvi' in col:
            mdata[prog_key].obs[col] = mdata[data_key].obs[col]
    for key in mdata[data_key].uns.keys():
        if '_scvi' in key:
            mdata[prog_key].uns[key] = mdata[data_key].uns[key]

    if not inplace: return mdata[prog_key]

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('-pk', '--prog_key', default='linear_SCVI', typ=str) 
    parser.add_argument('-dk', '--data_key', default='rna', typ=str) 
    parser.add_argument('-bk', '--batch_key', default=None, typ=str) 
    parser.add_argument('-lk', '--labels_key', default=None, typ=str) 
    parser.add_argument('--layer', default='X', type=str)
    parser.add_argument('--config_path', default='./linear_SCVI_config.gin', type=str)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    import mudata
    mdata = mudata.read(args.mudataObj_path)

    run_linear_SCVI(mdata, n_jobs=args.n_jobs, layer=args.layer, prog_key=args.prog_key, 
                    data_key=args.data_key, batch_key=args.batch_key, labels_key=args.labels_key, 
                    inplace=args.output, config_path=args.config_path)