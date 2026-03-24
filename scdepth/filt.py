#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation
#from scdepth import libraries

def filter_barcodes(
        df: pd.DataFrame, library: str, bin_div : int, read_mads: float | None = None,
        min_molecules: int | None = None, dry_run = False, skip_MT = False, 
        add_passed : bool = False, 

        #mt_mads: float | None = None, mt_cap: float | None = None, 
        **_) -> tuple[pd.DataFrame | None, dict]:

    tkey = 'is_cell' if 'is_cell' in df.columns else 'in_tissue'
    bco = min_molecules
    mads = read_mads

    #ldata = libraries.library2ns(library)
    mads = 3
    bco = 0
    if min_molecules is None:
        if library == '10X_visium_hd' or library == '10X_visium_hd_3p':
            if bin_div == 1:
                bco = 10
            elif bin_div < 4:
                bco = 20
            elif bin_div == 4:
                bco = 50
            else:
                bco = 100
            mads = 4
        else:
            bco = 500
    else:
        bco = min_molecules

    n_total = len(df)
    keep = np.ones(n_total, dtype=bool)
    molecules_all = df['molecules'].to_numpy()

    called_lost_mask = (df[tkey].to_numpy() <= 0)
    keep[called_lost_mask] = False
    called_lost = int(called_lost_mask.sum())

    info = {
        'total': n_total,
        'total_passed': 0,
        'called_lost': called_lost,
        'called_lost_perc': 100.0*called_lost/n_total,
        'called_key': tkey,
    }

    molecules_all = df['molecules'].to_numpy()
    baseline_lost_mask = (molecules_all < bco) & keep
    keep[baseline_lost_mask] = False
    baseline_lost = int(baseline_lost_mask.sum())

    kept_idx = np.flatnonzero(keep)
    rp = molecules_all[kept_idx].astype(float)
    rp = rp[rp > 0]
    if rp.size == 0:
        # nothing left
        if not dry_run:
            d = df.iloc[[]].copy()
            return d, info
        else:
            return None, info

    lm = np.log10(rp)
    med_lm = float(np.median(lm))
    mad_lm = float(median_abs_deviation(lm, scale=1.0))

    if mads <= 0:
        mad_cutoff = float(bco)
    else:
        if mad_lm == 0 or not np.isfinite(mad_lm):
            mad_cutoff = float(bco)
        else:
            mad_cutoff = float(max(bco, 10 ** (med_lm - mads * mad_lm)))

    mad_lost_mask = (molecules_all < mad_cutoff) & keep
    keep[mad_lost_mask] = False
    mad_lost = int(mad_lost_mask.sum())
    info['molecules_baseline_cutoff'] = float(bco)
    info['molecules_baseline_lost_perc'] = float(100*baseline_lost / n_total)
    info['molecules_mad_lost_perc'] = float(100*mad_lost / n_total)
    info['molecules_mad_cutoff'] = float(mad_cutoff)

    #if library != 'visium_hd' and not skip_MT and mt_cap is not None and mt_mads is not None:
    #    kept_idx = np.flatnonzero(keep)
    #    mol = df['molecules'].to_numpy(dtype=float)[kept_idx]
    #    mtm = df['mt_molecules'].to_numpy(dtype=float)[kept_idx]

    #    # avoid divide-by-zero / NaNs
    #    ok = np.isfinite(mol) & np.isfinite(mtm) & (mol > 0)
    #    mt_frac = np.full_like(mol, np.nan)
    #    mt_frac[ok] = mtm[ok] / mol[ok]
    #    mt_valid = mt_frac[ok]

    #    med = float(np.median(mt_valid))
    #    mad = float(median_abs_deviation(mt_valid, scale=1.0))
    #    if mad <= 0 or not np.isfinite(mad):
    #        mt_cutoff = float(mt_cap)
    #    else:
    #        mt_cutoff = float(min(mt_cap, med + mt_mads * mad))
    #    mt_lost_local = np.isfinite(mt_frac) & (mt_frac > mt_cutoff)
    #    # map back to global indices
    #    keep[kept_idx[mt_lost_local]] = False
    #    mt_lost = int(mt_lost_local.sum())
    #    info['mt_cutoff'] = float(mt_cutoff)
    #    info['mt_lost'] = int(mt_lost)
    #    info['mt_lost_perc'] = float(100 * mt_lost / n_total)

    tp = int(keep.sum())
    info['total_passed'] = tp
    info['total_passed_perc'] = float(tp / n_total * 100)
    if not dry_run:
        if add_passed:
            df['passed'] = keep.astype(int)
            return (df, info)
        d = df[keep].copy()
        return (d, info)
    return None, info

