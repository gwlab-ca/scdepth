#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

from scdepth.bindings import Downsampler, aggregate_visium_bins, calculate_gene_overlaps
from scdepth import fn, pl, fit, filt, libraries
from scdepth.cmds import common
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import os

use_filter = True
use_exc_genes = True
use_vis_hd = True

def custom_args(parser):
    g_init = parser.add_argument_group('Parameters')
    g_init.add_argument('-r', '--repeats', type=int, default=10,
                        help='Number of curves to generate')
    g_init.add_argument('--recoveries', type=str,
                        help='Recovery targets for downsampling (default: 30,40,50,60,70)')
    g_init.add_argument('--write', action='store_true',
                        help='Write barcode overlap data')
    g_init.add_argument('--min_obs', type=int, default=10,
                        help='Minimum number of genes with molecular count to keep for analyses')

    return parser

def build_parser(parser):
    common.prepare_args(parser, use_filter=use_filter, 
                      use_exc_genes=use_exc_genes, use_vis_hd=use_vis_hd,
                      command_args=custom_args)
    return parser

def resolve_args(args) -> dict:
    return {'repeats':args.repeats, 'recoveries':args.recoveries,
            'write':args.write, 'min_obs':args.min_obs}

def build_overlaps(args, full_summary):
    oprefix = f'{args.prefix}_genes'
    tprefix = f'{args.prefix}_genes_tmp'
    ldata = libraries.library2ns(args.library)
    is_hd = 'visium_hd' in ldata.library

    ds = Downsampler()
    ds.init(args.prefix, max_hist=args.max_hist, build_matrices=True, exclude_file=args.exclude_file)
    barcodes = fn.read_barcodes_meta(ds, args.prefix)

    if is_hd:
        ds.init_visium(rows=barcodes['array_row'].to_numpy(), cols=barcodes['array_col'].to_numpy(),
                    in_tissue=barcodes['in_tissue'].to_numpy(), countable=barcodes['countable'].to_numpy(), 
                    total=barcodes['total'].to_numpy(), total_rows=args.total_rows, total_cols=args.total_cols)
    else:
        args.bin_div = 0


    nbl, bline, _ = fit.baseline_fitter(args.prefix)

    if args.recoveries is None:
        args.recoveries = '30,40,50,60,70'


    recoveries = sorted(set(map(int, args.recoveries.split(','))))
    #saturations = sorted(set(map(int, args.saturations.split(','))))

    ds.downsample([1.0], umi_len=full_summary.umi_length, 
                seed=args.seed, threads=args.threads, aggregate_only=False)
    fstats = fn.bulk_stats(ds, full_summary)

    if is_hd:
        cdf, _ = aggregate_visium_bins(ds, 0, args.bin_div, args.bin_div)
        cdf = pd.DataFrame(cdf)
        fdf, _ = filt.filter_barcodes(cdf, **args)
    else:
        fdf = barcodes[barcodes['passed'] == 1].copy()

    idx = np.array(fdf.index.values, dtype=np.uint32)

    full_file = f'{tprefix}_full.gz'
    if not ds.write_gene_baseline(full_file, idx=idx, bin_div=args.bin_div):
        print('[error] Writing baseline file')
        return
    rng = np.random.default_rng(args.seed)
    seeds = rng.choice(2**32, size=args.repeats, replace=False)


    fracs = nbl.reads_for_recovery(recoveries) / full_summary.countable_reads
    #sat_fracs = nbl.reads_for_saturation(saturations) / full_summary.countable_reads
    fracs = fracs[fracs <= 1.0]
    slab = 'recovery'
    targs = recoveries
    #sat_fracs = sat_fracs[sat_fracs <= 1.0]



    seed_files = []
    all_stats  = []

    ratios = []

    for ci, s in enumerate(seeds):
        ds.downsample(fracs, umi_len=full_summary.umi_length,     
                    seed=s, threads=args.threads, aggregate_only=False)

        curve_stats = fn.bulk_stats(ds, full_summary)
        curve_stats.insert(loc=0, column=f'{slab}_target', value=targs)
        curve_stats.insert(loc=0, column='seed', value=s)
        curve_stats.insert(loc=0, column='curve_index', value=ci)
        curve_stats['read_ratio'] = fstats['reads'] / curve_stats['reads']
        ratios.append(fstats['reads'].values / curve_stats['reads'].values)
        all_stats.append(curve_stats)

        seed_file = f'{tprefix}_seed_{ci}.gz'
        if not ds.write_gene_mats(seed_file, idx=idx, bin_div=args.bin_div):
            print('[error] writing gene matrices')
            for s in seed_files:
                os.remove(s)
            os.remove(full_file)
            return
        seed_files.append(seed_file)

    all_stats  = pd.concat(all_stats)
    all_stats.insert(loc=0, column='sample', value=args.sample)
    all_stats.to_csv(f'{oprefix}_curve_stats.txt', sep='\t', index=False)

    gdf = pd.DataFrame(calculate_gene_overlaps(full_file, seed_files))
    if args.write:
        if is_hd:
            gdf.to_parquet(f'{oprefix}_barcode_summary.parquet')
        else:
            gdf.to_csv(f'{oprefix}_barcode_summary.txt.gz', sep='\t', float_format='%.5f')

    for s in seed_files:
        os.remove(s)

    color_jaccard = "#D55E00"
    color_stability = "#0072B2"

    fig, axs = pl.figax(1, 4, s=4)
    fig.subplots_adjust(hspace=0.4)
    axs = axs.flatten()
    sprefs = ['stability', 'stability_1', 'stability_2', 'stability_3p']
    cprefs = ['counts', 'counts_1', 'counts_2', 'counts_3p']
    fprefs = ['fracs', 'fracs_1', 'fracs_2', 'fracs_3p']
    gprefs = ['genes', 'genes_1', 'genes_2', 'genes_3p']
    fgprefs = ['', 'fgenes_1', 'fgenes_2', 'fgenes_3p']

    min_obs = args.min_obs
    kpercs = []


    rows = {slab:targs, 'ratio':np.array(ratios).mean(axis=0)}
    for r in cprefs + fprefs:
        rows[r] = []

    for r in sprefs + gprefs + fgprefs:
        if r != '':
            rows[r] = []

    for s, c, f, gp, fg, ax in zip(sprefs, cprefs, fprefs, gprefs, fgprefs, axs):
        g = gdf[gdf[c] >= min_obs]
        if c == 'counts':
            for i in range(len(targs)):
                rows[gp].append(np.median(gdf['counts'].values))
        else:
            for i in range(len(targs)):
                rows[gp].append(np.median(gdf[c].values))
                rows[fg].append(100 * np.median(gdf[c].values / gdf['counts'].values))


        g = gdf[gdf[c] >= min_obs]
        kpercs.append(100.0 * len(g) / len(gdf))
        for i in range(len(targs)):
            rows[c].append(len(g))
            rows[f].append(100.0 * len(g) / len(gdf))
            v = g[f'{s}_{i}'].values
            rows[s].append(np.median(v))
            pl.plot_violin(ax, pos=i, v=v, color=color_stability, alpha=0.8)

    odf = pd.DataFrame(rows)
    odf.to_csv(f'{oprefix}_summary.txt', sep='\t', index=False)

    for ax in axs:
        ax.set_ylim(0, 1)
        ax.set_xticks(np.arange(len(targs)))
        ax.set_xticklabels(targs)
        ax.tick_params(axis='both', labelsize=14)
        ax.tick_params(axis='both', which='both', length=4)
        ax.set_xlabel(f'{slab.title()} (%)', fontsize=24)
        ax.grid(axis='y', color='lightgray', ls='--', lw=1, which='both')
        ax.set_axisbelow(True)

    axs[0].set_ylabel('Mean Stability\nPer Barcode', fontsize=24, labelpad=15)

    axs[0].set_title(f'Global [{kpercs[0]:.2f}% Barcodes]', fontsize=24, pad=15)
    axs[1].set_title(f'Molecules = 1 [{kpercs[1]:.2f}% Barcodes]', fontsize=24, pad=15)
    axs[2].set_title(f'Molecules = 2 [{kpercs[2]:.2f}% Barcodes]', fontsize=24, pad=15)
    axs[3].set_title(f'Molecules = 3+ [{kpercs[3]:.2f}% Barcodes]', fontsize=24, pad=15)

    legend_elements = [
        Line2D([0], [0],
            color='black',
            marker='o',
            linestyle='None',
            markersize=12,
            label='Median')
    ]

    axs[-1].legend(handles=legend_elements, frameon=True, loc='lower right', fontsize=16)

    fig.savefig(f'{oprefix}_barcode_summary.svg', bbox_inches='tight')

    os.remove(full_file)

def main(parser, args):
    args, summary = common.resolve_args(parser, args, resolve_args=resolve_args, 
                                    use_filter=use_filter, 
                                    use_exc_genes=use_exc_genes, use_vis_hd=use_vis_hd)

    build_overlaps(args, summary)
