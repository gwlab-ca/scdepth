#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

from scdepth.fit import NBLibFit
import numpy as np
import numpy.typing as npt
import matplotlib.offsetbox as offsetbox

def plot_ztnb(ax, ztf : NBLibFit, hist : npt.NDArray[np.uint32], r, 
              fs : float = 9,
              s : float =2.5, lw : float =1.0, skip_legend : bool = False):

    bins = np.arange(1, len(hist) + 1)
    maxx = len(hist)

    ax.scatter(bins, hist/hist.sum(), s=40, edgecolor='gray', lw=1, label='Observed', color='#3b5b7a', zorder=6)
    #ax.bar(bins, hist/hist.sum(), width=1, edgecolor='gray', lw=1, label='Truth', color='#3b5b7a')
    ax.plot(bins, ztf.predict_pmf(len(hist)), color='k', label='ZT-NB Fit', zorder=4) #, 
            #lw=lw, marker='o', markersize=s, markeredgecolor='w', mew=0.6)


    mu = (ztf.rhat * (1 - ztf.phat)) / ztf.phat
    KS = ztf.KS if hasattr(ztf, 'KS') and np.isfinite(ztf.KS) else r.ztnb_KS
    text = (
        f'KS = {KS:.5f}\n'
        f'$\\hat{{r}}$ = {ztf.rhat:.3f}\n'
        f'$\\hat{{p}}$ = {ztf.phat:.3f}\n'
        f'$\\hat{{\\mu}}$ = {mu:.3f}\n'
        f'$\\hat{{L}}$ = {ztf.L / 10**6:.2f}M'
    )

    ob = offsetbox.AnchoredText(text, loc='upper right', frameon=True, prop={'size': fs})
    ob.patch.set(boxstyle='round', facecolor='white', alpha=0.8) 
    ax.add_artist(ob)
    ax.set_ylabel('PMF', fontsize=18)
    ax.set_xlabel('Reads per Molecule', fontsize=18)
    if not skip_legend:
        ax.legend(loc = 'lower right', fontsize=12, bbox_to_anchor=(1.0, 0.25))
    if maxx is not None:
        ax.set_xlim(0.5, maxx + 0.5)
    else:
        maxx = len(hist + 0.5)
    step = 2 if maxx < 16 else 5
    ax.set_xticks(np.arange(1, maxx + 1, step))
    ax.set_xticks(np.arange(1, maxx + 1, 1), minor=True)

def plot_fits(ax, fits, labels, curve_stats, ylim=None):
    ax.scatter(curve_stats.reads, curve_stats.saturation, s=35, lw=1.5, edgecolor='k', label='Truth', zorder=6)
    preds = []
    for l, f in zip(labels, fits):
        pred = f.predict(curve_stats.reads.values, truth=curve_stats.molecules.values)
        preds.append(pred)
        MAE = pred['saturation_err'].mean()
        xp = np.logspace(np.log10(curve_stats.reads.min()), np.log10(curve_stats.reads.max()), 500)
        yp = f.predict(xp)
        ax.plot(xp, yp['saturation'], label=f'{l} MAE = {MAE:.2f}pp', alpha=0.7, lw=1.5)
    ax.set_xscale('log') 
    ax.legend(loc = 'upper center', fontsize=10)
    ax.set_ylabel('Saturation (%)', fontsize=14)
    ax.set_xlabel('Countable Reads', fontsize=14)
    if ylim is not None:
        ax.set_ylim(*ylim)
