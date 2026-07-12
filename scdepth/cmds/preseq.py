#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

from scdepth.cmds import common
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os, argparse
import subprocess

def run_preseq_bound_pop(preseq : str, hist_path: str,
                out_path):

    if os.path.isfile(out_path):
        os.remove(out_path)

    cmd = [
        preseq, 'bound_pop',
        '-H', hist_path,
        '-o', out_path,
    ]

    subprocess.run(cmd, check=True) #, stdout=subprocess.DEVNULL)
    lines = open(out_path, 'r').readlines()
    return float(lines[1].split()[0])


def preseq_complexity(preseq : str, hist_path : str, out_path : str):
    defect = 0
    try:
        res = run_preseq_bound_pop(preseq, hist_path, out_path)
    except:
        defect = 1
        return np.nan, defect
    return res, defect


def run_preseq_lc_extra(preseq : str, hist_path: str,
                out_path: str, e_max_reads: int, step, use_D):
    if os.path.isfile(out_path):
        os.remove(out_path)
    if use_D:
        cmd = [
            preseq, 'lc_extrap',
            '-D', '-Q',
            '-H', hist_path,
            '-e', str(int(e_max_reads)),
            '-s', str(int(step)),
            '-o', out_path,
        ]
    else:
        cmd = [
            preseq, 'lc_extrap',
            '-Q', 
            '-H', hist_path,
            '-e', str(int(e_max_reads)),
            '-s', str(int(step)),
            '-o', out_path,
        ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def parse_preseq(out_path: str):
    R_list, M_list = [], []
    with open(out_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                r = float(parts[0])
                m = float(parts[1])
            except ValueError:
                continue
            if r <= 0:
                continue
            R_list.append(r)
            M_list.append(m)

    R = np.asarray(R_list, dtype=float)
    M = np.asarray(M_list, dtype=float)

    order = np.argsort(R)
    R, M = R[order], M[order]

    if np.any(np.diff(R) == 0):
        _, idx = np.unique(R, return_index=True)
        R, M = R[idx], M[idx]

    return R, M

def preseq_interp(preseq : str, hist_path: str, reads_target: np.ndarray,
                  out_path: str, e_max_reads: int | None = None, step: int = 100_000):
    reads_target = np.asarray(reads_target, dtype=float)
    if e_max_reads is None:
        e_max_reads = int(np.max(reads_target))

    defect = 0
    try:
        run_preseq_lc_extra(preseq, hist_path, out_path, e_max_reads=e_max_reads, step=step, use_D=False)
    except:
        defect = 1
        try:
            run_preseq_lc_extra(preseq, hist_path, out_path, e_max_reads=e_max_reads, step=step, use_D=True)
        except:
            return np.full(len(reads_target), np.nan), defect
    try:
        R, M = parse_preseq(out_path)
        x = np.log10(R)
        xt = np.log10(reads_target)
        M_interp = np.interp(xt, x, M)
        sat = 100.0 * (1.0 - (M_interp / reads_target))
    except:
        return np.full(len(reads_target), np.nan), 1


    return (sat, defect)

def build_parser(parser) -> argparse.ArgumentParser:
    g_init = parser.add_argument_group('Parameters')
    g_init.add_argument('--preseq', type=str, required=True,
                        help='Path to preseq executable')
    g_init.add_argument('-s', '--sample', type=str, default='sample',
                        help='Sample name to include in output')

    parser.add_argument('prefix', help='scdepth cache prefix (also used as output prefix)')
    return parser

def main(parser, args):
    oprefix, prefix = None, None
    oprefix = f'{args.prefix}_preseq'
    prefix = args.prefix
    preseq = args.preseq

    if not os.path.isfile(preseq):
        print('Preseq command not found')
        return

    rpms = pd.read_csv(prefix + '_limit_rpm.txt.gz', sep='\t')
    curves = pd.read_csv(prefix + '_limit_curves.txt', sep='\t')
    summary = pd.read_csv(prefix + '_limit_summary.txt', sep='\t')
    bins = rpms['bins'].to_numpy(int)

    htmp = oprefix + '_tmp_hist.txt'
    ptmp = oprefix + '_tmp_preseq.txt'

    cmin = summary['reads'].min()
    cmax = summary['reads'].max()
    cmin = 10 ** int(np.floor(np.log10(cmin)))

    psats = {}
    rows = []
    pout = {'limit':curves['limit'].values}
    for c in rpms.columns[1:]:
        v = rpms[c]
        _, ci, di = c.rsplit('_', maxsplit=2)
        ci = int(ci)
        di = int(di)
        lim = curves['limit'].values[di]
        Lhat = summary.loc[(summary.curve == ci) & (summary.limit == lim)].iloc[0].nb_lib_L
        fsat = summary.loc[(summary.curve == ci) & (summary.limit == lim)].iloc[0].saturation
        reads_target= curves[f'curve_{ci}_reads'].to_numpy(float)
        tsats = curves[f'curve_{ci}_saturation'].to_numpy(float)
        sats = curves[f'curve_{ci}_{di}_nb_lib'].to_numpy(float)
        with open(htmp, 'w') as fo:
            for b,h in zip(bins, v):
                fo.write(f'{b}\t{h}\n')

        psats, defect = preseq_interp(preseq, 
            hist_path=htmp, reads_target=reads_target,
            out_path=ptmp, e_max_reads=int(cmax),
            step=min(500_000, cmin)
        )

        pL, defect_L = preseq_complexity(preseq, htmp, ptmp)

        pout[f'curve_{ci}_{di}_sat'] = tsats
        pout[f'curve_{ci}_{di}_nb_lib'] = sats
        pout[f'curve_{ci}_{di}_preseq'] = psats
        pout[f'curve_{ci}_{di}_defect'] = defect
        rows.append((ci, lim, fsat, np.abs(tsats - sats).mean(), np.abs(tsats - psats).mean(), Lhat, pL, defect_L))

    df = pd.DataFrame(rows, columns=['curve', 'limit', 'saturation', 'nb_lib_mae', 'preseq_mae', 'nb_lib_L', 'preseq_L', 'preseq_L_defect'])
    cdf = pd.DataFrame(pout)
    df.to_csv(oprefix + '_summary.txt', sep='\t', index=False)
    cdf.to_csv(oprefix + '_curves.txt', sep='\t', index=False)
    df.insert(loc=0, column='sample', value=args.sample)
    cdf.insert(loc=0, column='sample', value=args.sample)

    if not os.path.isfile(f'{prefix}_limit_pilot_rpm.txt.gz'):
        return


    rpms = pd.read_csv(prefix + '_limit_pilot_rpm.txt.gz', sep='\t')
    psummary = pd.read_csv(prefix + '_limit_pilot_summary.txt', sep='\t')
    pcurves = pd.read_csv(prefix + '_limit_pilot_curves.txt', sep='\t')
    bins = rpms['bins'].to_numpy(int)

    cu = psummary[psummary['curve'] == 0].copy()
    pout = {}
    rows = []

    for c in rpms.columns[1:]:
        v = rpms[c]
        _, ci, di = c.rsplit('_', maxsplit=2)
        ci = int(ci)
        di = int(di)
        reads_target= pcurves[f'curve_{ci}_reads'].to_numpy(float)
        tsats = pcurves[f'curve_{ci}_saturation'].to_numpy(float)
        sats = pcurves[f'curve_{ci}_{di}_nb_lib'].to_numpy(float)
        with open(htmp, 'w') as fo:
            for b,h in zip(bins, v):
                fo.write(f'{b}\t{h}\n')

        step = min(250_000, cmin)
        psats, defect = preseq_interp(preseq, 
            hist_path=htmp, reads_target=reads_target,
            out_path=ptmp, e_max_reads=int(cmax + step),
            step=min(250_000, cmin)
        )

        pL, defect_L = preseq_complexity(preseq, htmp, ptmp)

        pout[f'curve_{ci}_{di}_sat'] = tsats
        pout[f'curve_{ci}_{di}_nb_lib'] = sats
        pout[f'curve_{ci}_{di}_preseq'] = psats
        pout[f'curve_{ci}_{di}_defect'] = defect
        rows.append((ci, cu.unit_limit.values[di], cu.limit.values[di], np.abs(tsats - sats).mean(), np.abs(tsats - psats).mean(), Lhat, pL, defect_L))

    pdf = pd.DataFrame(rows, columns=['curve','unit_limit', 'limit', 'nb_lib_mae', 'preseq_mae', 'nb_lib_L', 'preseq_L', 'preseq_L_defect'])
    pdf.insert(loc=0, column='sample', value=args.sample)
    pcdf = pd.DataFrame(pcurves)
    pcdf.insert(loc=0, column='sample', value=args.sample)
    pdf.to_csv(oprefix + '_pilot_summary.txt', sep='\t', index=False)
    pcdf.to_csv(oprefix + '_pilot_curves.txt', sep='\t', index=False)

    if os.path.isfile(ptmp):
        os.remove(ptmp)
    if os.path.isfile(htmp):
        os.remove(htmp)
