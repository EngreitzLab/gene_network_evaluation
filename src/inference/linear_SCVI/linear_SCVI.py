import os
import gin
import argparse

import mudata
import anndata
from tqdm.auto import tqdm

import logging
logging.basicConfig(level = logging.INFO)

def plot_training(model, ax):

    train_elbo = model.history["elbo_train"][1:]
    test_elbo = model.history["elbo_validation"]

    if ax is not None:
        train_elbo.plot(ax=ax)
    else:
        ax = train_elbo.plot()
    test_elbo.plot(ax=ax)

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

    Z_hat = model.get_latent_representation()
    loadings = model.get_loadings()

    return model, Z_hat, loadings

def run_linear_SCVI(mdata, batch_key=None, labels_key=None,
                    prog_key='linear_SCVI', data_key='rna',  
                    layer='X', config_path=None, n_jobs=1, 
                    inplace=True):

    """
    Perform gene program inference using linear SCVI.
    
    https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html

    ARGS:
        mdata : MuData
            mudata object containing anndata of data and cell-level metadata.
        batch_key: str
            cell metadata key (mdata[data_key].obs[batch_key]) under which batch labels are stored.
        labels_key: str
            cell metadata key (mdata[data_key].obs[labels_key]) under which celltype labels are stored.
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

    if layer=='X':
        layer_=None
    else:
        layer_=layer

    # Compute cNMF and create prog anndata
    # FIXME: Multi-processing crashes on HPC
    if n_jobs==-1:
        n_jobs=os.cpu_count()
    model, Z_hat, loadings = run_linear_SCVI_(mdata[data_key], layer=layer_,
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

    if not inplace: return mdata

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-bk', '--batch_key', default=None, typ=str) 
    parser.add_argument('-lk', '--labels_key', default=None, typ=str) 
    parser.add_argument('-pk', '--prog_key', default='linear_SCVI', typ=str) 
    parser.add_argument('-dk', '--data_key', default='rna', typ=str) 
    parser.add_argument('--layer', default='X', type=str)
    parser.add_argument('--config_path', default='./linear_SCVI_config.gin', type=str)
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    run_linear_SCVI(mdata, batch_key=args.batch_key, labels_key=args.labels_key, 
                    prog_key=args.prog_key, data_key=args.data_key, layer=args.layer, 
                    config_path=args.config_path, n_jobs=args.n_jobs, inplace=args.output)