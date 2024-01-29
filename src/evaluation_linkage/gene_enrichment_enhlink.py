import pandas as pd

import igraph

import numpy as np

import re

from collections import Counter

from statsmodels.stats.multitest import fdrcorrection

from collections import defaultdict

from argparse import ArgumentParser
from sys import stdout

from os.path import isfile


ARGPARSER = ArgumentParser(
    description='Compute enrichment score for gene set based on epigenomic linkage (enhancer-promoter)',
    prefix_chars='-'
)

ARGPARSER.add_argument('-b', '--bedpe',
                       required=True,
                       help="""Enhlink bedpe file containing enhancer-promoter linkage.
                       first 6 columns are enhancer-promoter coordinates.
                       Next three columns are p-value, score, and Adj. score.
                       The last column is the gene symbol of the promoter""",
                       type=str,
                       metavar='file')

ARGPARSER.add_argument('-g', '--genes',
                       required=True,
                       help='''Set(s) of genes to test for enrichment.
                       if a file is used as input, two columns tab separated are required,
                       with the first column as the gene symbol and the second the set name.
                       Otherwise, a list of genes is required''',
                       nargs="+",
                       type=str,
                       metavar='str/file')

ARGPARSER.add_argument('-m', '--motif',
                       required=False,
                       help='Motiff results file for each enhancer (FIMO file)',
                       type=str,
                       default=None,
                       metavar='file')

ARGPARSER.add_argument('-a', '--analysis',
                       required=False,
                       help='''Type of analysis to perform.
                       Default is "connectivity" that determines how connected are genes within the set.
                       Alternative is "motif" that determines the highest motif enrichment score
                       based on the enhancers' motifs''',
                       type=str,
                       default='connectivity',
                       metavar='{"connectivity", "motif"}')

ARGPARSER.add_argument('-bl', '--baseline',
                       required=True,
                       help='''list of baseline genes.
                       List of all the genes used in the analysis ''',
                       type=str,
                       default='',
                       metavar='file')

ARGPARSER.add_argument('-t', '--trials',
                       required=False,
                       help='Random trials to perform from the baseline genes to compute the p-value',
                       type=int,
                       default=1000,
                       metavar='int')

ARGPARSER.add_argument('-o', '--out',
                       required=True,
                       help='Output file',
                       type=str,
                       metavar='file')

ARGS = ARGPARSER.parse_args()


ANALYSIS = {"connectivity", "motif"}

assert(ARGS.analysis in ANALYSIS)

if ANALYSIS == "motif":
    assert(ARGS.motif is not None)


def main():
    """
    """
    baseline = load_baseline(ARGS.baseline)
    genes = load_gene_sets(ARGS.genes)
    linkage = load_bedpe(ARGS.bedpe)

    if ARGS.analysis == "connectivity":
        results = connectivity_analysis(linkage, genes, baseline, ARGS.trials)
    elif ARGS.analysis == "motif":
        results = motif_analysis(linkage, genes, baseline, ARGS.motif, ARGS.trials)

    results.to_csv(ARGS.out, sep="\t", index=None)
    print("\nfile written: {0}".format(ARGS.out))


def load_baseline(baseline_file):
    """
    """
    baseline = pd.read_csv(baseline_file, sep="\t", index_col=0, header=None)

    return baseline.index

def load_gene_sets(gene_input):
    """
    """
    gene_dict = {}

    if isinstance(gene_input, list) and len(gene_input) == 1:
        gene_input = gene_input[0]

    if isfile(gene_input):
        dfr = pd.read_csv(gene_input, sep="\t", header=None)
        dfr.columns = ["gene", "group", "p-value"]

        for group in set(dfr["group"]):
            index = dfr["group"] == group
            gene_dict[group] = dfr.loc[index, "gene"]
    else:
        if not isinstance(gene_input, list):
            gene_input = [gene_input]

        gene_dict["group_0"] = gene_input

    return gene_dict


def motif_analysis(linkage, genes, baseline, motif_file, trials=200):
    """
    """
    graph, nodes = create_linkage_graph(linkage)
    motifs = load_motif_file(motif_file)
    results = defaultdict(list)

    for group, gene_set in genes.items():
        stdout.write("\rProcessing group: {0}...".format(group))
        stdout.flush()

        subgraph, gene_set = get_subgraph(gene_set, nodes, graph)
        enhancers = set([v["name"] for v in subgraph.vs if v["type"] == "enhancer"])

        counter = Counter()

        for enh in enhancers:
            counter.update(motifs[enh])

        count_random = get_motif_count_for_random_genes(
        trials, baseline, gene_set, graph, nodes, motifs)

        pvals = Counter()

        for motif, count in counter.items():
            count_array = count_random[motif]
            count_array.append(count)
            pvals[motif] = 1.0 - np.argsort(count_array)[-1] / trials

        create_motif_results(pvals, gene_set, group, results)

    results = pd.DataFrame(results)
    results = results.sort_values(by="nb.motifs", ascending=False)

    return results


def create_motif_results(pvals, gene_set, group, results):
    """
    sort and filter motifs based on FDR pvals and test if TF is inside gene list
    """
    motif_arr, pval_arr = [], []

    for motif, pval in pvals.items():
        motif_arr.append(motif)
        pval_arr.append(pval)

    motif_arr = np.asarray(motif_arr)
    pval_arr = np.asarray(pval_arr)

    index = np.argsort(pval_arr)

    motif_arr = motif_arr[index]
    pval_arr = pval_arr[index]

    _, adj_pvals = fdrcorrection(pval_arr)

    top_pval = pval_arr[0]
    top_pval_adj = adj_pvals[0]

    sig_motifs = motif_arr[adj_pvals < 0.05]

    tfs = [x.split('_')[0].lower() for x in sig_motifs]
    tfs = [x for x in tfs if x in gene_set]

    results["group"].append(group)
    results["nb.motifs"].append(len(sig_motifs))
    results["min.pval"].append(top_pval)
    results["min.FDR"].append(top_pval_adj)
    results["motifs"].append(";".join(sig_motifs))
    results["TFs"].append(";".join(tfs))

    return results


def get_motif_count_for_random_genes(
        trials, baseline, gene_set, graph, nodes, motifs):
    """
    for each trial, sample random genes from baseline and count motifs within
    """
    count_random = defaultdict(list)

    print("\n")
    for trial in range(trials):
        stdout.write("\rtrial: {0}...".format(trial))
        stdout.flush()
        random_genes = np.random.choice(
            baseline, replace=False, size=len(gene_set))

        subgraph, random_genes = get_subgraph(random_genes, nodes, graph)
        enhancers = set([v["name"] for v in subgraph.vs if v["type"] == "enhancer"])

        counter_rand = Counter()

        for enh in enhancers:
            counter_rand.update(motifs[enh])

        for key, count in counter_rand.items():
            count_random[key].append(count)

    return count_random


def load_motif_file(motif_file):
    """
    """
    motifs = defaultdict(lambda:Counter())

    print("loading FIMO motif results file...")
    dfr = pd.read_csv(motif_file, sep="\s+", header=None, comment="#", index_col=None)
    dfr.columns = ["motifID", "enhancer", "start", "stop",
                   "strand", "score", "pval", "pval.adj", "motif"]

    for row in dfr.loc[:, ["motifID", "enhancer"]].itertuples():
        motifs[row.enhancer][row.motifID] += 1

    return motifs

def create_linkage_graph(linkage):
    """
    """
    print("creating linkage graph...")
    graph = igraph.Graph()
    nodes = set(list(linkage["gene"]) + list(linkage["enhancer"]))
    nodes = list(nodes)

    regex = re.compile("chr\w+:")
    node_type = ["enhancer" if regex.search(x) else "gene" for x in nodes]

    edges = linkage.loc[:, ["gene", "enhancer"]]

    graph.add_vertices(nodes)
    graph.vs["type"] = node_type
    graph.add_edges(edges.itertuples(index=False))

    return graph, set(nodes)

def connectivity_analysis(linkage, genes, baseline, trials=200):
    """
    """
    graph, nodes = create_linkage_graph(linkage)
    results = defaultdict(list)

    for group, gene_set in genes.items():
        stdout.write("\rProcessing group: {0}...".format(group))
        stdout.flush()
        subgraph, gene_set = get_subgraph(gene_set, nodes, graph)
        score, score_adj, linked_genes = get_number_of_connected_genes(subgraph, gene_set)

        scores = []
        # we compare the score inferred with an array of scores obtained
        # using sets of random genes
        for _ in range(trials):
            random_genes = np.random.choice(
                baseline, replace=False, size=len(gene_set))

            subgraph, random_genes = get_subgraph(random_genes, nodes, graph)
            rand_score, _, _ = get_number_of_connected_genes(subgraph, random_genes)
            scores.append(rand_score)

        if score == 0:
            pval = 1.0
        else:
            scores.append(score)
            scores = np.asarray(scores)
            pval = 1.0 - scores.argsort()[-1] / trials

            results["group"].append(group)
            results["p-value"].append(pval)
            results["score"].append(score)
            results["score_adj"].append(score_adj)
            results["linked_genes"].append(linked_genes)

    results = pd.DataFrame(results)
    results = results.sort_values(by="score_adj", ascending=False)

    return results


def get_subgraph(gene_set, nodes, graph):
    """
    """
    gene_set = [a.lower() for a in gene_set]
    gene_set = nodes.intersection(gene_set)
    neighbors = graph.neighborhood(gene_set)
    neighbors = set([b for a in neighbors for b in a])
    neighbors = [a["name"] for a in graph.vs[neighbors]]

    vs_set = gene_set.union(neighbors)

    subgraph = graph.subgraph(vs_set)

    return subgraph, gene_set


def get_number_of_connected_genes(subgraph, gene_set):
    """
    """
    if not gene_set:
        return 0.0, 0.0, ''

    gene_list = list(gene_set)
    mat = np.array(subgraph.distances(gene_list, gene_list))
    mat[mat == np.inf] = 0
    mat[mat > 1] = True
    xnz, ynz = np.nonzero(mat)

    linked_genes = np.array(gene_list)[list(set(xnz))]

    score = mat.sum() / 2
    # score adj. is the score normalized according to the number of elements
    # in the diagonal of the square matrix: nb_genes x nb_genes
    score_adj = score / (mat.shape[0] * mat.shape[0] / 2 - mat.shape[0])
    return score, score_adj, ';'.join(linked_genes)


def load_bedpe(bedpe):
    """
    """
    print("loading Enhlink linkage file...")
    dfr = pd.read_csv(bedpe, sep="\t", header=None, comment="#")

    dfr.columns = ["chrE", "startE", "stopE", "chrP", "startP", "stopP", "pval", "score", "adj.score", "gene"]
    dfr["enhancer"] = dfr["chrE"] + ":" + dfr["startE"].astype("str") + "-" + dfr["stopE"].astype("str")
    dfr["promoter"] = dfr["chrP"] + ":" + dfr["startP"].astype("str") + "-" + dfr["stopP"].astype("str")
    dfr["link"] = dfr["enhancer"] + "|" + dfr["promoter"]

    return dfr


if __name__ == '__main__':
    main()
