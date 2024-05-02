import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

import seaborn as sns
from matplotlib import pyplot as plt


def convertDatTraits(data):
    """
    get data trait module base on samples information

    :return: a dataframe contains information in suitable format for plotting module trait relationship heatmap
    :rtype: pandas dataframe
    """
    datTraits = pd.DataFrame(index=data.index)
    for i in range(data.shape[1]):
        data.iloc[:, i] = data.iloc[:, i].astype(str)
        if len(np.unique(data.iloc[:, i])) == 2:
            datTraits[data.columns[i]] = data.iloc[:, i]
            org = np.unique(data.iloc[:, i]).tolist()
            rep = list(range(len(org)))
            datTraits.replace(to_replace=org, value=rep,
                              inplace=True)
        elif len(np.unique(data.iloc[:, i])) > 2:
            for name in np.unique(data.iloc[:, i]):
                datTraits[name] = data.iloc[:, i]
                org = np.unique(data.iloc[:, i])
                rep = np.repeat(0, len(org))
                rep[np.where(org == name)] = 1
                org = org.tolist()
                rep = rep.tolist()
                datTraits.replace(to_replace=org, value=rep, inplace=True)

    return datTraits


def plot_topic_trait_relationship_heatmap(cell_topic_participation,
                                          metaData,
                                          annotation=False,
                                          save=True,
                                          show=True,
                                          file_format="pdf",
                                          file_name='topic-traitRelationships'):
    """
    plot topic-trait relationship heatmap

    :param cell_topic_participation: Anndata contains cell topic participation where cell information is obs and topic information is var
    :type cell_topic_participation: Anndata
    :param metaData: traits you would like to see the relationship with topics (must be column name of cell_participation.obs)
    :type metaData: list
    :param annotation: indicate if you want to add correlation and p_values as a text in each square (default:False)
    :type annotation: bool
    :param save: indicate if you want to save the plot or not (default: True)
    :type save: bool
    :param show: indicate if you want to show the plot or not (default: True)
    :type show: bool
    :param file_format: indicate the format of plot (default: pdf)
    :type file_format: str
    :param file_name: name and path of the plot use for save (default: topic-traitRelationships)
    :type file_name: str
    """
    datTraits = convertDatTraits(cell_topic_participation.obs[metaData])

    topicsTraitCor = pd.DataFrame(index=cell_topic_participation.var.index,
                                  columns=datTraits.columns,
                                  dtype="float")
    topicsTraitPvalue = pd.DataFrame(index=cell_topic_participation.var.index,
                                     columns=datTraits.columns,
                                     dtype="float")

    min_cell_participation = cell_topic_participation.to_df().min().min()
    for i in cell_topic_participation.var.index:
        for j in datTraits.columns:
            tmp = cell_topic_participation.to_df()[
                ~np.isclose(cell_topic_participation.to_df()[i],
                            min_cell_participation, atol=min_cell_participation)]

            tmp = stats.spearmanr(tmp[i], datTraits.loc[tmp.index, j], alternative='greater')
            topicsTraitCor.loc[i, j] = tmp[0]
            topicsTraitPvalue.loc[i, j] = tmp[1]

    topicsTraitCor.fillna(0.0, inplace=True)
    topicsTraitPvalue.fillna(1.0, inplace=True)

    for i in range(topicsTraitPvalue.shape[0]):
        rejected, tmp = fdrcorrection(topicsTraitPvalue.iloc[i, :])
        if not rejected.all():
            topicsTraitPvalue.iloc[i, :] = tmp

    xlabels = cell_topic_participation.to_df().columns
    ylabels = datTraits.columns

    if annotation:
        fig, ax = plt.subplots(figsize=(topicsTraitPvalue.shape[0] * 1.5,
                                        topicsTraitPvalue.shape[1] * 1.5), facecolor='white')

        # Loop over data dimensions and create text annotations.
        tmp_cor = topicsTraitCor.T.round(decimals=3)
        tmp_pvalue = topicsTraitPvalue.T.round(decimals=3)
        labels = (np.asarray(["{0}\n({1})".format(cor, pvalue)
                              for cor, pvalue in zip(tmp_cor.values.flatten(),
                                                     tmp_pvalue.values.flatten())])) \
            .reshape(topicsTraitCor.T.shape)

        sns.set(font_scale=1.5)
        res = sns.heatmap(topicsTraitCor.T, annot=labels, fmt="", cmap='RdBu_r',
                          vmin=-1, vmax=1, ax=ax, annot_kws={'size': 20, "weight": "bold"},
                          xticklabels=xlabels, yticklabels=ylabels)

    else:
        fig, ax = plt.subplots(figsize=(topicsTraitPvalue.shape[0],
                                        topicsTraitPvalue.shape[1]), facecolor='white')

        sns.set(font_scale=1.5)
        res = sns.heatmap(topicsTraitCor.T, cmap='RdBu_r',
                          vmin=-1, vmax=1, ax=ax, annot_kws={'size': 20, "weight": "bold"},
                          xticklabels=xlabels, yticklabels=ylabels)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize=20, fontweight="bold", rotation=90)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize=20, fontweight="bold")
    plt.yticks(rotation=0)
    ax.set_title(f"Topic-trait Relationships heatmap",
                 fontsize=30, fontweight="bold")
    ax.set_facecolor('white')

    if save:
        fig.savefig(f"{file_name}.{file_format}", bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()
