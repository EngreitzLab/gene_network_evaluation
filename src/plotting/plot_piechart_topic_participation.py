import pandas as pd
import random

import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt


def plot_piechart_topic_participation(cell_topic_participation,
                                      level,
                                      colors_topics=None,
                                      category=None,
                                      ascending=None,
                                      n=5,
                                      save=True,
                                      show=True,
                                      figsize=None,
                                      file_format="pdf",
                                      file_name="piechart_topicAvgCell"):
    """
    plot pie charts that shows contribution of each topics to each category (i.e cell type)

    :param cell_topic_participation: Anndata contains cell topic participation where cell information is obs and topic information is var
    :type cell_topic_participation: Anndata
    :param level: name of the column from cell_participation.obs
    :type level: str
    :param colors_topics: table contains color of each topic as a column called colors
    :type colors_topics: pandas dataframe
    :param category: list of items you want to plot pie charts which are subsets of cell_participation.obs[level](default: all the unique items in cell_participation.obs[level])
    :type category: list of str
    :param ascending: for each pie chart on which order you want to sort your data (default is descending for all pie charts)
    :type ascending: list of bool
    :param n: number of topics you want to annotate in pie charts (default: 5)
    :type n: int
    :param save: indicate if you want to save the plot or not (default: True)
    :type save: bool
    :param show: indicate if you want to show the plot or not (default: True)
    :type show: bool
    :param figsize: indicate the size of plot (default: (10 * (len(category) + 1), 10))
    :type figsize: tuple of int
    :param file_format: indicate the format of plot (default: pdf)
    :type file_format: str
    :param file_name: name and path of the plot use for save (default: piechart_topicAvgCell)
    :type file_name: str
    """

    if category is None:
        category = cell_topic_participation.obs[level].unique().tolist()

    if figsize is None:
        figsize = (10 * (len(category) + 1), 10)

    if cell_topic_participation.shape[1] <= n:
        n = max(0, cell_topic_participation.shape[1] - 2)

    fig, axs = plt.subplots(ncols=len(category) + 1,
                            figsize=figsize,
                            facecolor='white')

    if colors_topics is None:
        colors = sns.color_palette("turbo", cell_topic_participation.shape[1]).as_hex()

        def myfunction():
            return 0.1

        random.shuffle(colors, myfunction)
        index = cell_topic_participation.var.index.tolist()
        colors_topics = pd.DataFrame({'colors': colors}, index=index)

    colors = colors_topics

    if ascending is None:
        ascending = [False] * len(category)

    for i in range(len(category)):
        tissue = cell_topic_participation.obs[cell_topic_participation.obs[level] == category[i]]
        tmp = cell_topic_participation.to_df().loc[tissue.index, :]
        order = tmp.mean().sort_values(ascending=ascending[i]).index.tolist()
        index = tmp[order].sort_values(by=order, ascending=ascending[i]).index.tolist()
        tmp = tmp.reindex(columns=order)
        tmp = tmp.reindex(index)
        colors = colors.reindex(order)
        labels = tmp.columns.tolist()
        labels[n:] = ['' for i in range(len(labels) - n)]

        def make_autopct(values):
            def my_autopct(pct):
                if pct > values[n] * 100:
                    return '{p:.0f}%'.format(p=pct)
                else:
                    return ''

            return my_autopct

        axs[i].pie(tmp.mean(),
                   labels=labels,
                   colors=colors.colors.tolist(),
                   autopct=make_autopct(tmp.mean()),
                   wedgeprops={'linewidth': 0},
                   # labeldistance=0.8,
                   textprops={"fontsize": 25})

        axs[i].set_title(category[i], fontsize=40)

    colors = colors_topics
    handles = []
    for n in range(colors.shape[0]):
        patch = mpatches.Patch(color=colors.colors[n], label=colors.index[n])
        handles.append(patch)
        axs[len(category)].legend(loc='center left',
                                  title='Topic',
                                  ncol=3,
                                  handles=handles)
        axs[len(category)].axis('off')

    if save:
        fig.savefig(f"{file_name}.{file_format}")
    if show:
        plt.show()
    else:
        plt.close()
