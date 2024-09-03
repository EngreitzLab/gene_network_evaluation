import os
import sys
import Topyfic
import scanpy as sc


def make_top_model(trains, adata_path, n_top_genes=50, resolution=1, max_iter_harmony=10, min_cell_participation=None,
                   topmodel_output=None):
    if not os.path.isfile(adata_path):
        sys.exit("Anndata does not exist!")

    for train in trains:
        if not os.path.isfile(train):
            sys.exit(f"{train} does not exist!")

    # Check if the file exists
    if not os.path.exists(topmodel_output):
        os.makedirs(topmodel_output)

    adata = sc.read_h5ad(adata_path)

    top_model, clustering, adata_topmodel = Topyfic.calculate_leiden_clustering(trains=trains,
                                                                                data=adata,
                                                                                n_top_genes=n_top_genes,
                                                                                resolution=resolution,
                                                                                max_iter_harmony=max_iter_harmony,
                                                                                min_cell_participation=min_cell_participation)

    top_model.save_topModel(save_path=topmodel_output)

    adata_topmodel.write_h5ad(f"{topmodel_output}/topic_weight_umap.h5ad")
    clustering.to_csv(f"{topmodel_output}/topic_cluster_mapping.csv")
