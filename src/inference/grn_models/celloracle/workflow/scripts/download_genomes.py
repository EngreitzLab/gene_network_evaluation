from genomepy import install_genome
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--gdir', required=True)
parser.add_argument('-o','--organism', required=True)
args = vars(parser.parse_args())

# Get dir
gdir = args['gdir']
organism = args['organism']

# Install genomes
if organism == 'human':
    install_genome(name='hg38', genomes_dir=gdir, provider="UCSC")
elif organism == 'mouse':
    install_genome(name='mm10', genomes_dir=gdir, provider="UCSC")
