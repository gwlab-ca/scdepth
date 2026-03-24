#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter


def figax(r, c, s = None, w = None, h = None, dpi=150, **kwargs):
    if s is not None:
        kwargs['figsize'] = (c * s, r * s)
    elif w is not None and h is not None:
        kwargs['figsize'] = (c * w, r * h)
    fig, axs = plt.subplots(r, c, dpi=dpi, **kwargs)
    return fig, axs

def rand_jitter(N, max_diff, stdev = 0.120, seed=42):
    rng = np.random.default_rng(seed)
    ret = rng.normal(loc=0, size=N, scale=stdev)
    ret[ret > max_diff] = max_diff
    ret[ret < -max_diff] = -max_diff
    return ret

def jitter_box(ax, pos, pts, max_diff=0.4, scale=1, bwidth=0.95, color=None, s=3, lw=0, edgecolor='none', 
               label=None, cmap=None, c=None, alpha=1, vmin=None, vmax=None, marker='o',
               blw=1., mlw=1., wlw=1., plw=1.):

    boxprops = dict(linewidth=blw, color='black')
    medianprops = dict(linewidth=mlw, color='black')
    whiskerprops = dict(linewidth=wlw, color='black')
    capprops = dict(linewidth=plw, color='black')
    ax.scatter(rand_jitter(len(pts), max_diff) * scale + pos, pts, lw=lw, rasterized=True, 
               marker=marker, color=color, s=s, edgecolor=edgecolor, label=label, cmap=cmap, c=c,
               alpha=alpha, vmin=vmin, vmax=vmax)
    return ax.boxplot(pts, positions=[pos], widths=[bwidth], showfliers=False, 
                    boxprops=boxprops, medianprops=medianprops, 
                    whiskerprops=whiskerprops, capprops=capprops)

def jitter_boxh(ax, pos, pts, max_diff=0.4, scale=1, bwidth=0.95, color=None, s=3, lw=0, edgecolor='none', 
               label=None, cmap=None, c=None, alpha=1, vmin=None, vmax=None):
    ax.scatter(pts, rand_jitter(len(pts), max_diff) * scale + pos, lw=lw, rasterized=True, 
               marker='o', color=color, s=s, edgecolor=edgecolor, label=label, cmap=cmap, c=c,
               alpha=alpha, vmin=vmin, vmax=vmax)
    return ax.boxplot(pts, positions=[pos], widths=[bwidth], showfliers=False, 
                      medianprops={'color':'k'}, orientation='horizontal')

def add_contours(ax, df, xcol, ycol, xbins = 80, ybins=120, q_lo=0.60, q_hi=0.99, 
                n_levels=5, color='0.3', lw=1.2, alpha=0.75, xlog=True):
    if xlog:
        x_bins = np.logspace(np.log10(df[xcol].min()), np.log10(df[xcol].max()), xbins)
    else:
        x_bins = np.linspace(df[xcol].min(), df[xcol].max(), xbins)
    y_bins = np.linspace(df[ycol].min(), df[ycol].max(), ybins)
    H, xedges, yedges = np.histogram2d(df[xcol], df[ycol], bins=[x_bins, y_bins])
    H = gaussian_filter(H, sigma=(1.6, 1.2))
    Hpos = H[H > 0]

    lo = np.quantile(Hpos, q_lo)
    hi = np.quantile(Hpos, q_hi)
    levels = np.linspace(lo, hi, n_levels)

    xcenters = np.sqrt(xedges[:-1] * xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])

    #X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    X, Y = np.meshgrid(xcenters, ycenters)
    ax.contour(X, Y, H.T, levels=levels, colors=color, 
               linewidths=lw, alpha=alpha, antialiased=True)
