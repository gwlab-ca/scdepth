#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import numpy as np
import matplotlib.pyplot as plt
import numpy.typing as npt
from scdepth.pl.misc import figax
import numpy.typing as npt
import pandas as pd

mid_color = "#C9C9C9"
neg_color = "#2166AC"
pos_color = "#B2182B"
outline_gray = "#8A8A8A"


nmap = {'row':'Rows', 'col':'Columns'}

def _plot_scatter(ax, idx, y, sel, s, color, alpha, outline_color=None, label=None, lw=None, edgecolor=None):
    if outline_color is not None:
        ax.scatter(idx[sel], y[sel], s=s, lw=1.5, rasterized=True, color='none', edgecolor=outline_gray, alpha=1.0)
    ax.scatter(idx[sel], y[sel], s=s, rasterized=True, color=color, alpha=alpha, label=label, lw=lw, edgecolor=edgecolor)

def plot_effects(df : pd.DataFrame, sdf : pd.DataFrame, div : int, total_rows : int, total_cols : int):
    fig, axs = figax(2, 2, w=6, h=4)
    fig.subplots_adjust(wspace=0.2)

    vmax = 0
    for i, (k, tot) in enumerate(zip(('row','col'), (total_rows, total_cols))):
        d = df[(df['axis'] == k) & (df['included'] == 1)]
        sd = sdf[(sdf['axis'] == k)].iloc[0]
        idx = d['position'].values
        for j, (col, ax) in enumerate(zip(('ridge_rpm_dev_effect', 'sum_rpm_dev_effect'), axs[:,i])):
            y = d[col].values * 100.0
            _plot_scatter(ax, idx, y, (d['ridge_envelope'] == 0), s=25, alpha=0.7,
                        color=mid_color, outline_color=outline_gray, label = ('Within Envelope' if j == 0 else '_'))
            _plot_scatter(ax, idx, y, (d['ridge_envelope'] == -1), s=40, alpha=1, lw=1, edgecolor='k',
                        color=neg_color, outline_color=None, label = ('Pos. Envelope Exceeding' if j == 0 else '_'))
            _plot_scatter(ax, idx, y, (d['ridge_envelope'] == 1), s=40, alpha=1, lw=1, edgecolor='k',
                        color=pos_color, outline_color=None, label = ('Neg. Envelope Exceeding' if j == 0 else '_'))

            vmax = max(vmax, np.nanmax(np.abs(y)))
        dd = tot * 0.01
        for ax in axs[:,i]:
            ax.set_xlim(-dd, tot + dd)
        if k == 'row':
            axs[0,i].set_title(f'1 x {div} bins', fontsize=20)
        else:
            axs[0,i].set_title(f'{div} x 1 bins', fontsize=20)
        vv = d['ridge_threshold'].values
        vv = vv[np.isfinite(vv)][0] * 100
        axs[0,i].axhline(vv, ls='--', color='k', label='Familywise Null Envelope')
        axs[0,i].axhline(-vv, ls='--', color='k')

    for r in axs:
        for ax in r:
            ax.set_ylim(-vmax * 1.3, vmax * 1.3)
    for ax in axs[0]:
        ax.set_xticklabels([])
    axs[0,0].set_ylabel('Ridge\nRelative RPM (%)', fontsize=18)
    axs[1,0].set_ylabel('Direct\nRelative RPM (%)', fontsize=18)

    for i, k in enumerate(('row','col')):
        axs[-1,i].set_xlabel(f'{nmap[k]}', fontsize=20)

    lgds = []
    lgds.append(axs[0,-1].legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=16))
    return fig, axs, lgds
