import pandas as pd
import argparse, os, sys

"""
4th script for dictys pipeline. 
Assume that there exist a data directory and individual sub-folder for each cluster. There are TF-gene linking file for each individual cell clusters.
This script and many of the following scripts need to be ran for each cluster individually.

This script computes weights for the directed graph from TF to target genes using interaction mask inferred in the previous step.
"""


parser = argparse.ArgumentParser(description="")
parser.add_argument('--subdir', type=str, help="Individual cell cluster's directory")
parser.add_argument('--device', type=str, help="Whether using CPU (denoted as cpu) or GPU (denoted as cuda:N)")
parser.add_argument('--threads', type=int)


args = parser.parse_args()
d = args.subdir
dev = args.device
threads = args.threads


os.system(f"python3 -m dictys network reconstruct --device {dev} --nth {threads} {d}/expression.tsv.gz {d}/binlinking.tsv.gz {d}/net_weight.tsv.gz {d}/net_meanvar.tsv.gz {d}/net_covfactor.tsv.gz {d}/net_loss.tsv.gz {d}/net_stats.tsv.gz")
os.system(f"python3 -m dictys network normalize --nth {threads} {d}/net_weight.tsv.gz {d}/net_meanvar.tsv.gz {d}/net_covfactor.tsv.gz {d}/net_nweight.tsv.gz")
os.system(f"python3 -m dictys network indirect --nth {threads} --fi_meanvar {d}/net_meanvar.tsv.gz {d}/net_weight.tsv.gz {d}/net_covfactor.tsv.gz {d}/net_iweight.tsv.gz")

