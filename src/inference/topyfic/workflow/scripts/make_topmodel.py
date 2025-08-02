import os
import sys
import Topyfic
import pandas as pd
import scanpy as sc
import numpy as np
import time


def find_best_n_topic(n_topics, names, topmodel_outputs):
    df = pd.DataFrame(columns=['Name', 'k', 'N'])
    n = len(names)

    for i in range(n):
        for n_topic in n_topics:
            name = f"{topmodel_outputs[i]}/topModel_{names[i]}_{n_topic}.p"
            top_model = Topyfic.read_topModel(name)

            tmp = pd.DataFrame([[names[i], n_topic, top_model.N]], columns=['Name', 'k', 'N'])
            df = pd.concat([df, tmp])

    res = dict()
    for i in range(n):
        tmp = df[df.Name == names[i]]
        tmp['diff'] = tmp['N'] - tmp['K']
        tmp['diff'].idxmin()
        res[names[i]] = tmp.loc[tmp['diff'].idxmin(), 'K']
    return df, res


def make_top_model(trains, adata_paths, n_top_genes=50, resolution=1, max_iter_harmony=10, min_cell_participation=None,
                   topmodel_output=None):

    # Check if the file exists
    if not os.path.exists(topmodel_output):
        os.makedirs(topmodel_output, mode=0o777)

    adata = None
    print(adata_paths)
    if isinstance(adata_paths, str):
        adata = sc.read_h5ad(adata_paths)
    else:
        for adata_path in adata_paths:
            tmp = sc.read_h5ad(adata_path)
            if adata is None:
                adata = sc.read_h5ad(adata_path)
            else:
                adata = adata.concatenate(tmp)

    top_model, clustering, adata_topmodel = Topyfic.calculate_leiden_clustering(trains=trains,
                                                                                data=adata,
                                                                                n_top_genes=n_top_genes,
                                                                                resolution=resolution,
                                                                                max_iter_harmony=max_iter_harmony,
                                                                                min_cell_participation=min_cell_participation)

    print(top_model.name, topmodel_output)
    top_model.save_topModel(save_path=topmodel_output)

    adata_topmodel.write_h5ad(f"{topmodel_output}/topic_weight_umap.h5ad")
    clustering.to_csv(f"{topmodel_output}/topic_cluster_mapping.csv")

    analysis_top_model = Topyfic.Analysis(Top_model=top_model)
    analysis_top_model.calculate_cell_participation(data=adata)
    analysis_top_model.save_analysis(save_path=topmodel_output)
