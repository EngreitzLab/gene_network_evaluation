import numpy as np
import pandas as pd 

import seaborn as sns
from matplotlib import pyplot as plt

# Generic volcano plot
def plot_volcano(pvals, scores, labels=None, title=None, alpha=.05,
                 score_thresh=None, figsize=(9,9), fontsize=10, 
                 ax=None):
                     
    if labels is not None:
        plot_df = pd.DataFrame({'label': labels, 'pval': pvals, 'scores': scores})
    else:
        plot_df = pd.DataFrame({'pval': pvals, 'scores': scores})
        
    if score_thresh is not None:
        plot_df['is_sig'] = ((plot_df.pval < alpha) & (plot_df.scores > score_thresh)).astype(int)
    else:
        plot_df['is_sig'] = (plot_df.pval < alpha).astype(int)
    plot_df.is_sig = pd.Categorical(plot_df.is_sig, categories=[0,1])
    
    plot_df['log10p'] = plot_df.pval.apply(lambda x: -np.log10(x))
    
    if ax is None: fig, ax = plt.subplots(figsize=figsize)
    
    sns.scatterplot(x='scores', y='log10p', hue='is_sig', 
                    data=plot_df, ax=ax, linewidth=0 , 
                    legend=False, palette=['lightgray', 'white'])
    if 'label' in plot_df.columns:
        for idx in plot_df.index.values:
            ax.annotate(plot_df.loc[idx, 'label'], 
                        (plot_df.loc[idx, 'scores'],
                         plot_df.loc[idx, 'log10p']))
    
    
    ax.axhline(-np.log10(alpha), c='tomato', linestyle='--')
    if score_thresh is not None:
        ax.axvline(score_thresh, c='k', linestyle='--')
    
    ax.set_xlabel('statistic')
    ax.set_ylabel('-log10(p_value)')
    if title:
        ax.set_title(title)
        


        
    
    

  