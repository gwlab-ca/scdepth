#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

from scdepth.bindings import Downsampler
from scdepth import fn, pl, fit, libraries
from scdepth.cmds import common

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from types import SimpleNamespace
import os, warnings

use_filter = False
use_exc_genes = True
use_vis_hd = False

def custom_args(parser):
    g_init = parser.add_argument_group('Parameters')
    g_init.add_argument('-r', '--repeats', type=int, default=10,
                        help='Number of curves to generate')
    g_init.add_argument('--sat-start', type=float, default=10, 
                        help='Minimum saturation limit')
    g_init.add_argument('--plot-sats', type=str, default='10,30,50,60,70,80',
                        help='Saturation Limits to plot')
    g_init.add_argument('--sat-step', type=float, default=5,
                        help='Saturation step for limits')

    g_low = parser.add_argument_group('Minimum depth analysis')
    g_low.add_argument('--n-cells', type=int,
                        help='Number of cells to use for low depth analyses')
    g_low.add_argument('--use-scrna', action='store_true',
                        help='Use the filtered emptydrops to get the number of cells')
    g_low.add_argument('--tissue-frac', type=float,
                        help='Fraction of visium/visium HD slide covered with tissue (0-1.0)')
    g_low.add_argument('--tissue-scale', type=float,
                        help='Override scale factor for Visium/Visium HD')
    g_low.add_argument('--cell-scale', type=float,
                        help='Add an additional scaling factor for cells')
    g_low.add_argument('--base-start', type=int,
                        help='Override pilot depth start exponent (10^start) defaults: scRNA=2, Visium=6, Visium HD=6')
    g_low.add_argument('--base-end', type=int,
                        help='Override pilot depth end exponent (10^end) defaults: scRNA=5, Visium=8, Visium HD=8')
    g_low.add_argument('--capture-format', type=str, choices=['6.5mm', '11mm'],
                        help='Visium capture format: 11mm = 2.87, 6.5mm = 1.0 everything else 1.0')

    return parser

CAPTURE_SCALE = {
    '11mm':2.87,
    '6.5mm':1.0
}

def build_parser(parser):
    common.prepare_args(parser, use_filter=use_filter, 
                      use_exc_genes=use_exc_genes, use_vis_hd=use_vis_hd,
                      command_args=custom_args)
    return parser

def resolve_args(args) -> dict:
    return {'repeats':args.repeats, 'sat_start':args.sat_start,
            'plot_sats':args.plot_sats, 'sat_step':args.sat_step,
            'n_cells':args.n_cells, 'tissue_frac':args.tissue_frac,
            'tissue_scale':args.tissue_scale, 'base_start':args.base_start, 
            'base_end':args.base_end, 'capture_format':args.capture_format,
            'cell_scale':args.cell_scale, 'use_scrna':args.use_scrna}

def _make_grid(emin : int, emax : int) -> np.ndarray:
    vals = []
    for e in range(emin, emax + 1):
        b = 10 ** e
        vals.extend((b, 2*b, 5*b))
    return np.array(vals, dtype=int)

def build_pilot(args, ci : int, ds : Downsampler, full_summary : SimpleNamespace, seed : int, 
                truth : pd.DataFrame):
    if args.n_cells is None and args.tissue_frac is None:
        raise ValueError("One of --n-cells or --tissue-frac is required.")

    ureads, greads = None, None
    tscale = 1.0
    base_start = args.base_start
    base_end = args.base_end

    if base_start is None:
        base_start = 2 if args.n_cells is not None else 6

    if base_end is None:
        base_end = 5 if args.n_cells is not None else 8
    ureads = _make_grid(base_start, base_end)
    if args.n_cells is not None:
        cell_scale = getattr(args, "cell_scale")
        if cell_scale is None:
            cell_scale = 1.0
        if cell_scale <= 0:
            raise ValueError("cell_scale must be > 0")
        greads = np.round(ureads * float(args.n_cells) * float(cell_scale))
        tscale = cell_scale
    elif args.tissue_frac is not None:
        if args.tissue_scale is not None:
            tscale = args.tissue_scale
        elif args.capture_format in CAPTURE_SCALE:
            tscale =  CAPTURE_SCALE[args.capture_format]
        if tscale <= 0:
            raise ValueError("tissue scale must be > 0")
        greads = np.round(ureads * float(args.tissue_frac) * float(tscale))

    greads = greads
    fracs = greads / full_summary.countable_reads
    keep = fracs <= 1.0
    ureads = ureads[keep]
    greads = greads[keep]
    fracs = fracs[keep]
    if len(fracs) == 0:
        print('[error] All the requested low depths are less than the total depth')
        return None


    ds.downsample(fracs, umi_len=full_summary.umi_length, 
                seed=seed, threads=args.threads, aggregate_only=True)

    curve_stats = fn.bulk_stats(ds, full_summary)
    curve_stats.insert(loc=0, column='seed', value=seed)
    curve_stats.insert(loc=0, column='curve_index', value=ci)
    meta = {
        'seed':seed,
        'curve_index':ci,
        'base_start':base_start, 
        'base_end':base_end, 
        'unit_min':ureads.min(),
        'unit_max':ureads.max(),
        'reads_min':greads.min(),
        'reads_max':greads.max(),
        'frac_min':fracs.min(),
        'frac_max':fracs.max(),
        'scale':tscale, 
        'mode': 'scrna' if args.n_cells else 'spatial'
    }

    if args.n_cells is not None:
        meta['n_cells'] = args.n_cells
    elif args.capture_format is not None:
        meta['capture_format'] = args.capture_format

    ds_summary = []

    hists = {}
    p_curves = {}
    p_curves[f'curve_{ci}_reads'] = truth['reads'].values
    p_curves[f'curve_{ci}_saturation'] = truth['saturation'].values

    for i, (l, g, r) in enumerate(zip(ureads, greads, curve_stats.itertuples())):
        nbl = fit.NBLibFit()
        hist = fn.get_rpm_hist(ds=ds)[i]
        nbl.fit(hist, reads=r.reads, molecules=r.molecules)
        hists[f'rpm_{ci}_{i}'] = hist
        pred1 = nbl.predict(truth['reads'], truth=truth['molecules'])
        row = [ci, l, g, nbl.KS, nbl.nll_per_mol, nbl.tail_mass, pred1['saturation_err'].mean() , 
                r.saturation, r.reads, r.molecules, r.fraction, nbl.rhat, nbl.phat, nbl.L]


        p_curves[f'curve_{ci}_{i}_nb_lib'] = pred1['saturation']

        ds_summary.append(row)

    col_names = ['curve','unit_limit','limit', 'ztnb_ks', 'ztnb_nll', 'ztnb_tail_mass', 'nb_lib_err', 
                'saturation', 'reads', 'molecules', 'fraction', 'ztnb_rhat', 'ztnb_phat', 'nb_lib_L']

    ds_summary = pd.DataFrame(ds_summary, columns=col_names)
    return curve_stats, ds_summary, p_curves, hists, meta

def process_pilot(args, oprefix, pilot_stats, pilot_summary, pilot_meta):
    stats_df  = pd.concat(pilot_stats)
    summary_df  = pd.concat(pilot_summary)
    meta_df = pd.DataFrame(pilot_meta)

    stats_df.insert(loc=0, column='sample', value=args.sample)
    summary_df.insert(loc=0, column='sample', value=args.sample)
    meta_df.insert(loc=0, column='sample', value=args.sample)

    stats_df.to_csv(f'{oprefix}_pilot_stats.txt', sep='\t', index=False)
    summary_df.to_csv(f'{oprefix}_pilot_summary.txt', sep='\t', index=False)
    meta_df.to_csv(f'{oprefix}_pilot_curve_meta.txt', sep='\t', index=False)

    fig, axs = pl.figax(1, 4, w=6, h=4)
    gmed = summary_df[['unit_limit', 'saturation', 'reads', 'ztnb_ks', 'nb_lib_err']].groupby('unit_limit').median()
    axs[0].scatter(summary_df['unit_limit'], summary_df['ztnb_ks'], label='Seeds')
    axs[0].plot(gmed.index, gmed['ztnb_ks'], lw=2, color='r', label='Median KS')
    axs[0].set_ylabel('ZT-NB KS')
    axs[0].set_xscale('log')

    axs[1].scatter(summary_df['unit_limit'], summary_df['nb_lib_err'], label='Seeds')
    axs[1].plot(gmed.index, gmed['nb_lib_err'], lw=2, color='r', label='Median MAE (pp)')
    axs[1].set_xscale('log')
    ymax = summary_df['nb_lib_err'].max() * 1.1
    axs[1].set_ylim(0, ymax)
    axs[1].set_ylabel('ZT-NB KS')
    axs[1].set_ylabel('NB-Lib MAE (pp)')

    axs[2].scatter(summary_df['unit_limit'], summary_df['saturation'], label='Seeds')
    axs[2].plot(gmed.index, gmed['saturation'], lw=2, color='r', label='Median Saturation (%)')
    axs[2].set_xscale('log')
    ymax = summary_df['saturation'].max() * 1.1
    axs[2].set_ylim(0, ymax)
    axs[2].set_ylabel('Saturation (%)')

    axs[3].scatter(summary_df['unit_limit'], summary_df['reads'], label='Reads')
    axs[3].plot(gmed.index, gmed['reads'], lw=2, color='r', label='Median Countable Reads')
    axs[3].set_xscale('log')
    axs[3].set_yscale('log')
    axs[3].set_ylabel('Countable Reads')

    for ax in axs:
        ax.set_xlabel('Effective Reads')
        ax.legend()

    axs[1].set_ylabel('NB-Lib MAE (pp)')
    fig.savefig(f'{oprefix}_pilot_summary.svg', bbox_inches='tight')


def build_limits(args, full_summary):
    oprefix = None
    oprefix = f'{args.prefix}_limit'

    keep_limits = set(map(int, args.plot_sats.split(',')))
    ldata = libraries.library2ns(args.library)
    calc_sau = (ldata.probe_based == 0)

    ds = Downsampler()
    ds.init(args.prefix, max_hist=args.max_hist, build_matrices=False, exclude_file=args.exclude_file,
            calc_sau=calc_sau)

    if not os.path.isfile(args.prefix + '_fit_baseline.txt'):
        print('The baseline fit does not exist please run scdepth fit first')
        return

    nbl, bstats, full_stats = fit.baseline_fitter(args.prefix)

    ds.downsample([bstats.fraction], umi_len=full_summary.umi_length, 
            seed=args.seed, threads=args.threads, aggregate_only=True)

    if args.n_cells is None and args.use_scrna:
        barcodes = fn.read_barcodes_meta(ds, args.prefix)
        #df = fn.barcode_df(ds=ds, df=barcodes.copy(), step=0)
        args.n_cells = barcodes['passed'].sum()

    limits = np.arange(args.sat_start, full_stats.saturation + 1e-12, args.sat_step)
    lfracs = fit.saturation_fracs(nbl, full_summary, limits) #, full=full_stats)
    lfracs = lfracs[lfracs['fraction'] <= 1.0].reset_index(drop=True)

    rng = np.random.default_rng(args.seed)
    seeds = rng.choice(2**32, size=args.repeats, replace=False)
    limits = lfracs['saturation'].values
    ds_curves = {'limit':[l for l in limits]}
    ds_summary = []
    ds_hists = []
    all_stats = []

    pilot_stats = []
    pilot_summary = []
    pilot_meta = []
    pilot_hists = []
    pilot_curves = []


    warnings.filterwarnings(
        "ignore",
        message="invalid value encountered in subtract",
        category=RuntimeWarning,
        module="scipy.optimize._numdiff"
    )

    for ci, s in enumerate(seeds):
        ds.downsample(lfracs['fraction'].values, umi_len=full_summary.umi_length, 
                    seed=s, threads=args.threads, aggregate_only=True)
        curve_stats = fn.bulk_stats(ds, full_summary)
        curve_stats.insert(loc=0, column='seed', value=s)
        curve_stats.insert(loc=0, column='curve_index', value=ci)

        all_stats.append(curve_stats)

        ds_curves[f'curve_{ci}_reads'] = curve_stats.reads
        ds_curves[f'curve_{ci}_saturation'] = curve_stats.saturation
        hists = {}


        for i, (l, r) in enumerate(zip(limits, curve_stats.itertuples())):
            nbl = fit.NBLibFit()
            hist = fn.get_rpm_hist(ds=ds)[i]
            hists[f'total_rpm_{ci}_{i}'] = hist

            nbl.fit(hist, reads=r.reads, molecules=r.molecules)


            pred1 = nbl.predict(curve_stats['reads'], truth=curve_stats['molecules'])
            ds_curves[f'curve_{ci}_{i}_nb_lib'] = pred1['saturation']
            hc = hist.sum()
            row = [ci, l, nbl.KS, nbl.nll_per_mol, nbl.tail_mass, pred1['saturation_err'].mean() , 
                   r.saturation, r.reads, r.molecules, r.fraction, nbl.rhat, nbl.phat, nbl.L]

            ds_summary.append(row)
        ds_hists.append(hists)

        if args.tissue_frac is not None or args.n_cells is not None:
            pstats, psum, pcurve, phist, pmeta = build_pilot(args, ci, ds, full_summary, s, curve_stats)
            pilot_hists.append(phist)
            pilot_stats.append(pstats)
            pilot_summary.append(psum)
            pilot_meta.append(pmeta)
            pilot_curves.append(pcurve)

    col_names = ['curve','limit','ztnb_ks', 'ztnb_nll', 'ztnb_tail_mass', 'nb_lib_err', 'saturation', 
                'reads', 'molecules', 'fraction', 'ztnb_rhat', 'ztnb_phat', 'nb_lib_L']
    ds_summary = pd.DataFrame(ds_summary, columns=col_names)
    ds_curves  = pd.DataFrame(ds_curves)
    all_stats  = pd.concat(all_stats)

    hh = {'bins':np.arange(1, args.max_hist)}
    for h in ds_hists:
        hh.update(h)
    hdf = pd.DataFrame(hh)
    hdf.to_csv(f'{oprefix}_rpm.txt.gz', sep='\t', index=False)

    ds_summary.insert(loc=0, column='sample', value=args.sample)
    ds_curves.insert(loc=0, column='sample', value=args.sample)
    all_stats.insert(loc=0, column='sample', value=args.sample)

    ds_summary.to_csv(f'{oprefix}_summary.txt', sep='\t', index=False)
    ds_curves.to_csv(f'{oprefix}_curves.txt', sep='\t', index=False)
    all_stats.to_csv(f'{oprefix}_curve_stats.txt', sep='\t', index=False)

    fig, axs = pl.figax(1, 2, w=6, h=4)
    gmed = ds_summary[['limit', 'ztnb_ks', 'nb_lib_err']].groupby('limit').median()
    axs[0].scatter(ds_summary['limit'], ds_summary['ztnb_ks'], label='Seeds')
    axs[0].plot(gmed.index, gmed['ztnb_ks'], lw=2, color='r', label='Median KS')
    axs[0].set_ylabel('ZT-NB KS')
    axs[1].scatter(ds_summary['limit'], ds_summary['nb_lib_err'], label='Seeds')
    axs[1].plot(gmed.index, gmed['nb_lib_err'], lw=2, color='r', label='Median MAE (pp)')

    ymax = ds_summary['nb_lib_err'].max() * 1.1
    axs[1].set_ylim(0, ymax)

    for ax in axs:
        ax.set_xlabel('Saturation Limit')
        ax.legend()

    axs[1].set_ylabel('NB-Lib MAE (pp)')
    axs[1].set_title(f'{args.sample} Full Saturation = {full_stats.saturation:.2f}% Limit = {limits[-1]:.2f}%')
    fig.savefig(f'{oprefix}_summary.svg', bbox_inches='tight')

    fig, ax = pl.figax(1, 1, w=6, h=3)
    ci = 0
    idx = 0
    for i, l in enumerate(limits):
        if i < (len(limits) - 1) and (l not in keep_limits):
            continue
        ax.plot(ds_curves[f'curve_{ci}_reads'], ds_curves[f'curve_{ci}_{i}_nb_lib'], color=plt.cm.tab10(idx), alpha=1, label=f'{l:.1f}% Limit', lw=2)
        idx += 1
    ax.set_ylim(0, 100)
    ax.set_xscale('log')
    ax.scatter(ds_curves[f'curve_{ci}_reads'], ds_curves[f'curve_{ci}_saturation'], color='0.7', edgecolor='k', lw = 1,
            s=25, label='Truth', zorder=5)
    ax.legend(fontsize=8)
    ax.set_ylabel('NB-Lib\nSaturation (%)', fontsize=14)
    fig.savefig(f'{oprefix}_curve_0.svg', bbox_inches='tight')

    if len(pilot_stats) > 0:
        process_pilot(args, oprefix, pilot_stats, pilot_summary, pilot_meta)
        hh = {'bins':np.arange(1, args.max_hist)}
        for h in pilot_hists:
            hh.update(h)
        hdf = pd.DataFrame(hh)
        hdf.to_csv(f'{oprefix}_pilot_rpm.txt.gz', sep='\t', index=False)

        pc = {}
        for c in pilot_curves:
            pc.update(c)
        ppdf = pd.DataFrame(pc)
        ppdf.to_csv(f'{oprefix}_pilot_curves.txt', sep='\t', index=False)



def main(parser, args):
    args, summary = common.resolve_args(parser, args, resolve_args=resolve_args, 
                                    use_filter=use_filter, use_exc_genes=use_exc_genes, 
                                    use_vis_hd=use_vis_hd)

    build_limits(args, summary)
