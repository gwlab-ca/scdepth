#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import argparse
from scdepth.bindings import Downsampler
from scdepth.cmds import common
from scdepth import fn, pl, fit, libraries
import numpy as np
import pandas as pd

use_filter = False
use_exc_genes = True
use_vis_hd = False

def custom_args(parser : argparse.ArgumentParser) -> None:
    g_init = parser.add_argument_group('Parameters')
    g_init.add_argument('-b', '--baseline-sat', type=float, default=50.0,
                        help='Baseline saturation for curves libraries over max_sat they will be downsampled')
    g_init.add_argument('-m', '--max-sat', type=float, default=55.0,
                        help='Libraries with saturation greater than this will be downsampled to the baseline for MM-NB estimation')
    g_init.add_argument('-p', '--points', type=int, default=8,
                        help = 'Number of points for curve building')
    g_init.add_argument('--sat-end', type=float, default=90,
                        help='Curve points maximum saturation')
    g_init.add_argument('--sat-start', type=float, default=10,
                        help='Curve points minimum saturation')
    g_init.add_argument('--sat-step', type=float, default=1,
                        help='Step to output read costs')

def resolve_args(args) -> dict:
    return {
        'baseline_sat':args.baseline_sat,
        'max_sat':args.max_sat,
        'points':args.points,
        'sat_start':args.sat_start,
        'sat_end':args.sat_end,
        'sat_step':args.sat_step,
        #'plot_sats': args.plot_sats
    }

def build_parser(parser):
    common.prepare_args(parser, use_filter=use_filter, 
                      use_exc_genes=use_exc_genes, use_vis_hd=use_vis_hd,
                      command_args=custom_args)
    return parser

def build_curve(args, full_summary):
    oprefix = None
    print(args)
    oprefix = f'{args.prefix}_fit'


    ldata = libraries.library2ns(args.library)
    calc_sau = (ldata.probe_based == 0)
    ds = Downsampler()
    ds.init(args.prefix, max_hist=args.max_hist, build_matrices=False, exclude_file=args.exclude_file,
            calc_sau=calc_sau)

    full_stats, bfrac = fit.find_target_saturation(ds, full_summary, max_sat=args.max_sat, 
                        target_sat=args.baseline_sat, threads=args.threads, seed=args.seed) 

    ds.downsample([bfrac], umi_len=full_summary.umi_length, 
            seed=args.seed, threads=args.threads, 
            aggregate_only=True)

    bstats = fn.bulk_stats(ds, full_summary)
    bstats.insert(loc=0, column='sample', value=args.sample)
    r = bstats.iloc[0]

    nbl = fit.NBLibFit()
    nbl.fit(fn.get_rpm_hist(ds=ds)[0], reads=r.reads, molecules=r.molecules)
    full_fracs = fit.full_fracs(nbl, full_summary, args.sat_start, full_stats.loc[0, 'saturation'], True, steps=args.points)

    bhist = fn.get_rpm_hist(ds)[0].copy()
    hists = {
        'bin':np.arange(1, len(bhist) + 1),
        'total':bhist,
    }

    sau_params = {}
    if calc_sau:
        for k in fn.READ_CLASSES[:-1]:
            kr = r[f'{k}_reads']
            km = r[f'{k}_molecules']
            khist = fn.get_rpm_hist(ds, key=k)[0].copy()
            hists[k] = khist
            nbc = fit.NBLibFit()
            nbc.fit(khist, reads=kr, molecules=km)
            for kk, v in nbc.serialize().items():
                sau_params[f'{k}_{kk}'] = v


    ds.downsample(list(full_fracs['fraction'].to_numpy(float)), umi_len=full_summary.umi_length, 
                seed=args.seed, threads=args.threads, aggregate_only=True)

    if calc_sau:
        for k in fn.READ_CLASSES[:-1]:
            khist = fn.get_rpm_hist(ds, key=k)[-1].copy()
            hists[f'full_{k}'] = khist

    hists['full_total'] = fn.get_rpm_hist(ds)[-1].copy()

    curve_stats = fn.bulk_stats(ds, full_summary)
    curve_stats['target_saturation'] = full_fracs['saturation']
    curve_stats.insert(loc=0, column='sample', value=args.sample)

    fig, (ax, fax) = pl.figax(1, 2, w=6, h=4)
    fig.subplots_adjust(hspace=0.3, wspace=0.15)

    lab = 'NB-Lib'
    preds = nbl.predict(curve_stats['reads'].to_numpy(int), truth=curve_stats['molecules'].to_numpy(int))

    pl.fits.plot_ztnb(ax, nbl, bhist, r, fs=10, lw=2, s=6)
    pl.fits.plot_fits(fax, [nbl], [lab], curve_stats)

    fig.savefig(f'{oprefix}_fits.svg', bbox_inches='tight')

    fit_params = nbl.serialize()
    fit_params['nb_lib_saturation_MAE'] = preds['saturation_err'].mean()
    fit_params['nb_lib_recovery'] = nbl.predict_recovery(r.reads),

    fit_params['countable_reads'] = full_summary.countable_reads
    fit_params['total_reads'] = full_summary.total_reads
    fit_params['full_saturation'] = full_stats.iloc[0].saturation
    fit_params['full_reads'] = full_stats.iloc[0].reads
    fit_params['full_molecules'] = full_stats.iloc[0].molecules

    hdf = pd.DataFrame(hists)
    hdf.to_csv(f'{oprefix}_rpm_histograms.txt', sep='\t', index=False)

    for k, v in fit_params.items():
        bstats[k] = v
    if calc_sau:
        for k, v in sau_params.items():
            bstats[k] = v
        for k in fn.READ_CLASSES[:-1]:
            bstats[f'full_{k}_reads'] = full_stats[f'{k}_reads'].values[0]
            bstats[f'full_{k}_molecules'] = full_stats[f'{k}_molecules'].values[0]
    bstats.to_csv(f'{oprefix}_baseline.txt', sep='\t', index=False)

    curve_stats['nb_lib_molecules'] = preds['molecules']
    curve_stats['nb_lib_recovery'] = nbl.predict_recovery(curve_stats['reads'].values)
    curve_stats.to_csv(f'{oprefix}_curve.txt', sep='\t', index=False)

    sats = np.arange(args.sat_start, args.sat_end + args.sat_step, args.sat_step)

    costs = pd.DataFrame(
        {
            'saturation':sats,
            'reads':nbl.reads_for_saturation(sats),
        }, 
    )
    arate = full_summary.countable_reads / full_summary.total_reads 
    costs['molecules'] = np.round(nbl.predict_molecules(costs['reads'])).astype(int)
    costs['sequenced_reads'] = np.round(costs['reads'] / arate).astype(int)
    costs['recovery'] = nbl.predict_recovery(costs['reads'].values)
    costs.to_csv(f'{oprefix}_saturation_costs.txt', sep='\t', index=False)

    costs = pd.DataFrame(
        {
            'recovery':sats,
            'reads':nbl.reads_for_recovery(sats),
        }, 
    )
    costs['molecules'] = np.round(nbl.predict_molecules(costs['reads'])).astype(int)
    costs['sequenced_reads'] = np.round(costs['reads'] / arate).astype(int)
    costs['saturation'] = nbl.predict_saturation(costs['reads'].values)
    costs.to_csv(f'{oprefix}_recovery_costs.txt', sep='\t', index=False)

def main(parser, args):
    args, summary = common.resolve_args(parser, args, resolve_args=resolve_args, 
                                        use_filter=use_filter, 
                                        use_exc_genes=use_exc_genes, use_vis_hd=use_vis_hd)
    build_curve(args, summary)
