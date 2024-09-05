import sys
import Topyfic
import scanpy as sc
import os


def make_single_train_model(name, adata_path, k, random_state, train_output):
    """
    Train model for only one random state

    :param name: name of the dataset
    :type name: str
    :param adata_path: path of the input anndata
    :type adata_path: str
    :param k: number of topic
    :type k: int
    :param random_state: random state used to train
    :type random_state: int
    :param train_output: path of train object
    :type train_output: str
    """

    if not os.path.isfile(adata_path):
        sys.exit("Anndata does not exist!")

    # Check if the file exists
    if not os.path.exists(train_output):
        os.makedirs(train_output, mode=0o777)

    adata = sc.read_h5ad(adata_path)

    train = Topyfic.Train(name=f"{name}_{k}_{random_state}",
                          k=k,
                          n_runs=1,
                          random_state_range=[random_state])
    train.run_LDA_models(adata, batch_size=128, max_iter=5, n_jobs=1, n_thread=1)

    train.save_train(save_path=train_output)


def make_train_model(name, adata_path, k, n_runs, random_state, train_output):
    """
    Combine all Train models that were built using different random seed

    :param name: name of the dataset
    :type name: str
    :param adata_path: path of the input anndata
    :type adata_path: str
    :param k: number of topic
    :type k: int
    :param n_runs: number of single LDA models
    :type n_runs: int
    :param random_state: list of random seeds that were used to train LDA
    :type random_state: list
    :param train_output: path of train object
    :type train_output: str
    """

    if not os.path.isfile(adata_path):
        sys.exit("Anndata does not exist!")

    # Check if the file exists
    if not os.path.exists(train_output):
        os.makedirs(train_output, mode=0o777)

    adata = sc.read_h5ad(adata_path)

    main_train = Topyfic.Train(name=f"{name}_{k}",
                               k=k,
                               n_runs=n_runs,
                               random_state_range=random_state)
    trains = []
    for i in random_state:
        print(f"{train_output}train_{name}_{k}_{i}.p")
        train = Topyfic.read_train(f"{train_output}train_{name}_{k}_{i}.p")
        print(train.random_state_range)
        trains.append(train)

    print(main_train.random_state_range)
    main_train.combine_LDA_models(adata, single_trains=trains)
    print(f"{train_output}{main_train.name}")
    main_train.save_train(save_path=train_output)

