from genomepy import install_genome
import argparse

# Init args
parser = argparse.ArgumentParser(
    prog="python download_genomes.py",
    description="Download reference genomes for genomepy."
)
parser.add_argument('-d','--genome_dir', required=True, type=str, help='where the reference genome genomepy uses will be saved')
parser.add_argument('-g','--genome', required=True, type=str, help='version of the genome to download, e.g. hg38, mm10')
args = vars(parser.parse_args())

# Get dir
genome_dir = args['genome_dir']
genome = args['genome']

# Install genomes
install_genome(name=genome, genomes_dir=genome_dir, provider="UCSC")
