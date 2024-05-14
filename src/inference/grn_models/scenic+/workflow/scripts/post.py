import pandas as pd
import numpy as np
import mudata as mu
import scanpy as sc
import argparse

# Init args
parser = argparse.ArgumentParser()

# Parse args
path_eReg_direct = "/cellar/users/aklie/data/datasets/neurips2021_small/analysis/scenic+/2024_05_08/eRegulon_direct.tsv"
path_data = "/cellar/users/aklie/data/datasets/neurips2021_small/small.h5mu"
path_r2g = "/cellar/users/aklie/data/datasets/neurips2021_small/analysis/scenic+/2024_05_08/r2g.csv"
path_grn = "/cellar/users/aklie/data/datasets/neurips2021_small/analysis/scenic+/2024_05_08/grn.csv"
path_out = "/cellar/users/aklie/data/datasets/neurips2021_small/analysis/scenic+/2024_05_08/scenic+.h5mu"

# Save data
data.write_h5mu(path_mdata)
