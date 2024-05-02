import numpy as np
import pandas as pd
import random
from scipy.cluster.hierarchy import ward, leaves_list

import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt


def plot_structure(cell_topic_participation,
                   level,
                   colors_topics=None,
                   category=None,
                   topic_order=None,
                   ascending=None,
                   metaData=None,
                   metaData_palette=None,
                   width=None,
                   n=2,
                   order_cells=['hierarchy'],
                   save=True,
                   show=True,
                   figsize=None,
                   file_format="pdf",
                   file_name="structure_topicAvgCell"):
    """
    plot structure which shows contribution of each topics for each cells in given categories

    :param cell_topic_participation: Anndata contains cell topic participation where cell information is obs and topic information is var
    :type cell_topic_participation: Anndata
    :param colors_topics: table contains color of each topic as a column called colors
    :type colors_topics: pandas dataframe
    :param level: name of the column from cell_participation.obs
    :type level: str
    :param category: list of items you want to plot which are subsets of cell_participation.obs[level](default: all the unique items in cell_participation.obs[level])
    :type category: list of str
    :param topic_order: indicate if you want to have a specific order of topics which it should be name of topics. if None, it's gonna sort by cell participation
    :type topic_order: list of str
    :param ascending: for each structure plot on which order you want to sort your data (default is descending for all structure plot)
    :type ascending: list of bool
    :param metaData: if you want to add annotation for each cell add column name of that information (make sure you have that inforamtion in your cell_participation.obs)
    :type metaData: list
    :param metaData_palette: color palette for each metaData you add
    :type metaData_palette: dict
    :param width: width ratios of each category (default is based on the number of the cells we have in each category)
    :type width: list of int
    :param n: number of topics you want to sum if you used order_cell == 'sum' (default: 2)
    :type n: int
    :param order_cells: determine which kind of sorting options you want to use ('sum', 'hierarchy', sort by metaData); sum: sort cells by sum of top n topics; hierarchy: sort data by doing hierarchical clustring; metaData sort by metaData (default: ['hierarchy'])
    :type order_cells: list
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

    if colors_topics is None:
        colors = sns.color_palette("turbo", cell_topic_participation.shape[1]).as_hex()

        def myfunction():
            return 0.1

        random.shuffle(colors, myfunction)
        index = cell_topic_participation.var.index.tolist()
        colors_topics = pd.DataFrame({'colors': colors}, index=index)

    check = ["hierarchy", "sum"]
    if metaData is not None:
        check = check + metaData
    for order_cell in order_cells:
        if order_cell not in check:
            return "order cell was not valid"

    a = []
    for i in range(len(category)):
        a.append(cell_topic_participation.obs[cell_topic_participation.obs[level] == category[i]].shape[0])
    a.append(min(a) / 2)
    if width is None:
        width = a
    if metaData is None:
        fig, axs = plt.subplots(nrows=1,
                                ncols=len(category) + 1,
                                figsize=figsize,
                                gridspec_kw={'width_ratios': width},
                                facecolor='white')
    else:
        b = [8] + [1] * len(metaData)
        fig, axs = plt.subplots(nrows=len(metaData) + 1,
                                ncols=len(category) + 1,
                                figsize=figsize,
                                gridspec_kw={'width_ratios': width,
                                             'height_ratios': b},
                                facecolor='white')

    colors = colors_topics

    if ascending is None:
        ascending = [True] * len(category)

    for i in range(len(category)):
        tissue = cell_topic_participation.obs[cell_topic_participation.obs[level] == category[i]]
        tmp = cell_topic_participation.to_df().loc[tissue.index, :]
        if topic_order is None:
            order = tmp.mean().sort_values(ascending=False).index.tolist()
        else:
            order = topic_order
        # tmp['non_major'] = 1 - tmp.sum(axis=1)
        # tmp.non_major[tmp['non_major'] < 0] = 0
        index = tmp[order].sort_values(by=order, ascending=False).index.tolist()
        tmp = tmp.reindex(index)
        if len(order_cells) == 1:
            order_cell = order_cells[0]
            if order_cell == 'hierarchy':
                Z = ward(tmp)
                index = leaves_list(Z).tolist()
                tmp = tmp.iloc[index, :]
            elif order_cell == 'sum':
                index = tmp[order[:n]].sum(axis=1).sort_values(ascending=ascending[i]).index.tolist()
                tmp = tmp.reindex(index)
            elif order_cell in metaData:
                index = tissue.sort_values(by=[order_cell], ascending=ascending[i]).index.tolist()
                tmp = tmp.reindex(index)
        else:
            order_cell = order_cells[:-1]
            tissue.sort_values(by=order_cell, ascending=ascending[i], inplace=True)
            tmp = tmp.reindex(tissue.index)
            groups = tissue[order_cell].value_counts(sort=False).reset_index()
            count = 0
            index = []
            for j in range(groups.shape[0]):
                if groups.loc[j, 'count'] != 0:
                    sub = tissue.iloc[count:count + groups.loc[j, 'count'], :]
                    sub = tmp.loc[sub.index, :]
                    if sub.shape[0] > 1:
                        Z = ward(sub)
                        sub = leaves_list(Z).tolist()
                        sub = [x + count for x in sub]
                        index = index + sub
                    else:
                        index = index + [count]
                    count = count + groups.loc[j, 'count']

            tmp = tmp.iloc[index, :]

        # order.append('non_major')
        tmp = tmp.reindex(columns=order)

        colors = colors.reindex(order)
        if metaData is None:
            axs[i].stackplot(tmp.index.tolist(),
                             tmp.T.to_numpy(),
                             labels=tmp.columns.tolist(),
                             colors=colors.colors.tolist(),
                             linewidths=0)

            axs[i].xaxis.set_ticks([])
            axs[i].set_title(category[i], fontsize=40)
            axs[i].set_ylim(0, 1)
            axs[i].set_xlim(0, a[i])
            axs[0].set_ylabel("Topic proportion", fontsize=25)
        else:
            axs[0, i].stackplot(tmp.index.tolist(),
                                tmp.T.to_numpy(),
                                labels=tmp.columns.tolist(),
                                colors=colors.colors.tolist(),
                                linewidths=0)

            axs[0, i].xaxis.set_ticks([])
            axs[0, i].set_title(category[i], fontsize=40)
            axs[0, i].set_ylim(0, 1)
            axs[0, i].set_xlim(0, a[i])
            axs[0, 0].set_ylabel("Topic proportion", fontsize=30)

            tissue = tissue[metaData]
            tissue = tissue.reindex(tmp.index.tolist())
            for j in range(len(metaData)):
                if type(metaData_palette[metaData[j]]) == dict:
                    tissue.replace(metaData_palette[metaData[j]], inplace=True)

            x = [i for i in range(tmp.shape[0])]
            y = np.repeat(3000, len(x))
            for j in range(len(metaData)):
                color = tissue[metaData[j]].values
                if type(metaData_palette[metaData[j]]) == dict:
                    axs[j + 1, i].scatter(x, y,
                                          label=metaData_palette[metaData[j]],
                                          c=color,
                                          s=1000,
                                          marker="|",
                                          alpha=1,
                                          edgecolor='none')
                else:
                    axs[j + 1, i].scatter(x, y,
                                          label=metaData_palette[metaData[j]],
                                          c=color,
                                          cmap=metaData_palette[metaData[j]].get_cmap(),
                                          s=1000,
                                          marker="|",
                                          alpha=1,
                                          edgecolor='none')

                axs[j + 1, i].axis('off')
                axs[j + 1, i].set_xlim(0, a[i])

    colors = colors_topics

    handles = []
    for n in range(colors.shape[0]):
        patch = mpatches.Patch(color=colors.colors[n], label=colors.index[n])
        handles.append(patch)
        if metaData is None:
            axs[len(category)].legend(loc='center left',
                                      title='Topic',
                                      ncol=3,
                                      handles=handles)
            axs[len(category)].axis('off')
        else:
            axs[0, len(category)].legend(loc='center left',
                                         title='Topic',
                                         ncol=4,
                                         handles=handles)
            axs[0, len(category)].axis('off')

    if metaData is not None:
        for j in range(len(metaData)):
            handles = []
            if type(metaData_palette[metaData[j]]) == dict:
                for met in metaData_palette[metaData[j]].keys():
                    patch = mpatches.Patch(color=metaData_palette[metaData[j]][met], label=met)
                    handles.append(patch)
                    axs[j + 1, len(category)].legend(loc='center left',
                                                     title=metaData[j].capitalize(),
                                                     ncol=4,
                                                     handles=handles)
                    axs[j + 1, len(category)].axis('off')
            else:
                clb = fig.colorbar(mappable=metaData_palette[metaData[j]],
                                   ax=axs[j + 1, len(category)],
                                   orientation='horizontal',
                                   fraction=0.9)
                clb.ax.set_title(metaData[j].capitalize())
                axs[j + 1, len(category)].axis('off')

    if save:
        fig.savefig(f"{file_name}.{file_format}")
    if show:
        plt.show()
    else:
        plt.close()
