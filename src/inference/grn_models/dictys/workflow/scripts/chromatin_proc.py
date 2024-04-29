import pandas as pd
import argparse, os, sys

"""
3rd script for dictys pipeline. 
Assume that there exist a data directory and individual sub-folder for each cluster. There are bam file and peaks for each cluster of cells.
This script and many of the following scripts need to be ran for each cluster individually.

This script computes TF-region link (wellington and homer), Region-gene link (TSS distance), and finally TF-gene link (binlinking).
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument('--subdir', type=str, help="Individual cell cluster's directory")
parser.add_argument('--motif', type=str, help="")
parser.add_argument('--genome', type=str, help="")
parser.add_argument('--gene_annotation', type=str, help="")
parser.add_argument('--threads', type=int)
args = parser.parse_args()
d = args.subdir
motif = args.motif
genome = args.genome
annot = args.gene_annotation
threads = args.threads

os.system(f'python3 -m dictys chromatin wellington {d}/reads.bam {d}/reads.bai {d}/peaks.bed {d}/footprints.bed --nth {threads}')
os.system(f'python3 -m dictys chromatin homer {d}/footprints.bed {motif} {genome} {d}/expression.tsv.gz {d}/motifs.bed {d}/wellington.tsv.gz {d}/homer.tsv.gz')
    
os.system(f'python3 -m dictys chromatin binding {d}/wellington.tsv.gz {d}/homer.tsv.gz {d}/binding.tsv.gz')
os.system(f'python3 -m dictys chromatin tssdist {d}/expression.tsv.gz {d}/wellington.tsv.gz {annot} {d}/tssdist.tsv.gz')
os.system(f'python3 -m dictys chromatin linking {d}/binding.tsv.gz {d}/tssdist.tsv.gz {d}/linking.tsv.gz')
os.system(f'python3 -m dictys chromatin binlinking {d}/linking.tsv.gz {d}/binlinking.tsv.gz 20')
