import os
import argparse

import mudata

import numpy as np
import pandas as pd

import decoupler as dc

from joblib import Parallel, delayed
from tqdm.auto import tqdm


ENRICHMENT_METHODS = {
    "aucell": dc.run_aucell,
    "gsea": dc.run_gsea,
    "gsva": dc.run_gsva,
    "mdt": dc.run_mdt,
    "mlm": dc.run_mlm,
    "ora": dc.run_ora,
    "udt": dc.run_udt,
    "ulm": dc.run_ulm,
    "viper": dc.run_viper,
    "wmean": dc.run_wmean,
    "wsum": dc.run_wsum,
    "consensus": dc.run_consensus,
}


def compute_tf_enrichment(
    mdata: mudata.MuData,
    grn_key: str = "grn",
    method: str = "aucell",
    prog_key: str = "tf_enrichment",
    celltype: str = None,
    data_key: str = "rna",
    layer: str = None,
    min_n: int = 5,
    inplace: bool = True,
    seed: int = 1234,
):
    # Extract GRN
    grn = mdata.uns[grn_key].copy()
    if celltype is not None:
        grn = grn.query("celltype == @celltype")

    # Extract data
    adata = mdata[data_key].copy()
    if layer is not None:
        adata.X = adata.layers[layer].copy()

    # Run enrichment
    ENRICHMENT_METHODS[method](
        mat=adata,
        net=grn,
        source="source",
        target="target",
        min_n=min_n,
        use_raw=False,
        seed=seed
    )

    # Get activities
    acts = dc.get_acts(adata, obsm_key=f'{method}_estimate')

    # Create adjacency matrix
    adj = grn.pivot_table(index="target", columns="source", values="weight", fill_value=0)
    adj = adj[acts.var_names]
    adj = adj.reindex(adata.var_names, fill_value=0)
    acts.varm['loadings'] = adj.T.values

    # Add to mdata
    mdata[prog_key] = acts
    if not inplace:
        mdata = mudata.MuData({prog_key: mdata[prog_key].copy(),
                               data_key: mdata[data_key].copy()})
        return mdata
    