#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import numpy as np
import pandas as pd
from scipy.optimize import root_scalar
from scdepth.bindings import Downsampler
import sys
from numpy.typing import NDArray
from scipy import optimize
from scipy.stats import nbinom
from scipy.special import logit
from types import SimpleNamespace
import scdepth.fn as fn
#from .ztnbfast import _ztnb_nll_numba_wrapper

def log1mexp(a):
    if a > -0.6931471805599453:
        return np.log(-np.expm1(a))
    else:
        return np.log1p(-np.exp(a))

def _ztnb_nll(params, k_vals, freq):
    log_r, log_mu = params
    r = np.exp(log_r)
    mu = np.exp(log_mu)
    p = r / (r + mu)

    log_pmf_k = nbinom.logpmf(k_vals, r, p)
    log_p0 = nbinom.logpmf(0, r, p)
    log_trunc = log1mexp(log_p0)
    return -np.sum(freq * (log_pmf_k - log_trunc))

class NBLibFit(object):
    rhat : float = np.nan
    phat : float = np.nan
    L : float = np.nan
    L_init : float = np.nan

    KS : float = np.nan
    error_flag : bool = True
    nll_per_mol : float = np.nan
    tail_mass : float = np.nan

    def __init__(self, rhat : float = np.nan, L : float = np.nan, 
                 phat : float = np.nan):
        self.rhat = rhat
        self.phat = phat
        self.L = L

    def fit(self, rpm_hist : NDArray[np.uint32], reads : int, molecules : int) -> bool:
        R = float(reads)
        M = float(molecules)
        self._fit_ztnb(rpm_hist)
        if self.error_flag:
            return False

        self.L_init = self._estimate_L(R)
        self.L = self.L_init

        return True

    def serialize(self) -> dict:
        d = {
            'ztnb_rhat':float(self.rhat),
            'ztnb_phat':float(self.phat),
            'ztnb_KS':float(self.KS),
            'ztnb_nll_per_mol':float(self.nll_per_mol),
            'ztnb_tail_mass':float(self.tail_mass),
            'ztnb_error_flag':int(self.error_flag),
        }
        d['nb_lib_L'] = self.L_init
        return d


    def unserialize(self, d: dict):
        self.rhat = float(d['ztnb_rhat'])
        self.phat = float(d['ztnb_phat'])
        self.L_init = float(d['nb_lib_L'])
        self.L = self.L_init

        self.KS = float(d['ztnb_KS'])
        self.nll_per_mol = float(d['ztnb_nll_per_mol'])
        self.tail_mass = float(d['ztnb_tail_mass'])
        self.error_flag = int(d['ztnb_error_flag']) > 0

    def predict(self, reads, truth = None) -> pd.DataFrame:
        reads = np.asarray(reads)
        pred = self.predict_molecules(reads)
        psat = (100 * (1 - pred/reads))
        prec = 100 * pred / self.L
        df = pd.DataFrame({
            'reads':reads, 
            'molecules':pred,
            'saturation':psat,
            'recovery':100 * pred / self.L
        })
        if truth is not None:
            truth = np.asarray(truth)
            tsat = 100 * (1 - truth / reads)
            trec = 100 * truth / self.L
            df[f'molecules_truth'] = truth
            df[f'molecules_pdiff'] = 100.0 * np.abs(truth - pred) / truth
            df[f'saturation_truth'] = tsat
            df[f'saturation_err'] = np.abs(psat - tsat)
            df[f'recovery_truth'] = trec
            df[f'recovery_err'] = np.abs(prec - trec)
        return df

    def predict_rpm(self, bins : int, M : int) -> NDArray[np.double]:
        K = np.arange(1, bins + 1)
        pmf = nbinom.pmf(K, self.rhat, self.phat)
        p0  = nbinom.pmf(0, self.rhat, self.phat)
        pmf_ztnb = pmf / (1.0 - p0)
        pmf_ztnb = pmf_ztnb / pmf_ztnb.sum()
        return M * pmf_ztnb

    def predict_pmf(self, bins : int) -> NDArray[np.double]:
        K = np.arange(1, bins + 1)
        pmf = nbinom.pmf(K, self.rhat, self.phat)
        p0  = nbinom.pmf(0, self.rhat, self.phat)
        pmf_ztnb = pmf / (1.0 - p0)
        pmf_ztnb = pmf_ztnb / pmf_ztnb.sum()
        return pmf_ztnb

    def predict_molecules(self, reads):
        return self._nb_model(reads)

    def predict_saturation(self, reads):
        M = self._nb_model(reads)
        return 100.0 - 100.0 * M/reads

    def predict_recovery(self, reads, squeeze : bool = True):
        if np.isscalar(reads):
            reads = [reads]
        reads = np.asarray(reads)
        ret = np.array([self._reads2recovery(s) for s in reads])
        if squeeze:
            return ret[0] if len(ret) == 1 else ret
        return ret

    def reads_for_recovery(self, recs, squeeze : bool = True):
        if np.isscalar(recs):
            recs = [recs]
        recs = np.asarray(recs) / 100.0
        ret = np.array([self._recovery2reads(s) for s in recs])
        if squeeze:
            return ret[0] if len(ret) == 1 else ret
        return ret

    def reads_for_saturation(self, sats, squeeze=True):
        if np.isscalar(sats):
            sats = [sats]
        sats = np.asarray(sats) / 100.0
        ret = np.array([self._sat2reads(s) for s in sats])
        if squeeze:
            return ret[0] if len(ret) == 1 else ret
        return ret

    def _recovery2reads(self, rec: float) -> float:
        return (self.rhat * self.L) * ((1.0 - rec)**(-1.0 / self.rhat) - 1.0)

    def _reads2recovery(self, R: float) -> float:
        #M = float(self._nb_model(R))
        p = (self.L * self.rhat) / (R + self.L * self.rhat)
        return 100.0 * (1.0 - p**self.rhat)

    def _nb_model(self, R, L : float | None = None):
        if L is None:
            L = self.L
        R = np.asarray(R, dtype=np.double)
        out = L * (1.0 - np.power(self.rhat / (self.rhat + R / L), self.rhat))
        return float(out) if out.ndim == 0 else out

    def _sat2reads(self, S_target: float, tol: float = 1e-9, max_iter: int = 200) -> float:
        y = 1.0 - S_target  # target yield = M/R
        def g(R):
            return self._nb_model(R) - y * R

        Rhi = max(self.L, 10.0)
        Rlo = 1.0
        for _ in range(100):
            if g(Rhi) <= 0:
                break
            Rhi *= 2.0
        else:
            raise RuntimeError("Failed to bracket root; check L, r, and S_target.")

        # Ensure we bracket: g(R_lo) >= 0 and g(R_hi) <= 0

        if g(Rlo) < 0:
            Rlo = 0.0
            if g(Rlo) < 0:
                raise RuntimeError("Unexpected: g(0) < 0; cannot bracket.")

        if g(Rhi) > 0:
            raise RuntimeError("Upper bracket not sufficient; increase R_hi.")

        # Bisection
        lo, hi = float(Rlo), float(Rhi)
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            val = g(mid)
            if abs(val) < tol:
                return mid
            if val > 0:
                lo = mid
            else:
                hi = mid

        return 0.5 * (lo + hi)

    def _estimate_L(self, R : float) -> float:
        R = float(R)
        return R * self.phat / (self.rhat * (1 - self.phat))

    def _fit_ztnb(self, hist, eps=1e-12, rmin=0.25):
        hist = np.asarray(hist, dtype=float)
        k = np.arange(1, len(hist) + 1, dtype=int)

        mu_min, mu_max = 1e-6, 200 
        rmax = 50
        bounds = [(np.log(rmin), np.log(rmax)), (np.log(mu_min), np.log(mu_max))]

        # Multi-start grid (small, robust)
        starts_r = np.log(np.array([0.5, 1.0, 2.0, 5.0, 10.0]))
        starts_mu = np.log(np.array([0.25, 0.5, 1.0, 2.0, 5.0, 10.0]))

        best = None
        for lr in starts_r:
            for lmu in starts_mu:
                x0 = np.array([lr, lmu], float)
                res = optimize.minimize(
                    _ztnb_nll, x0=x0, args=(k, hist),
                    method="L-BFGS-B", bounds=bounds
                )
                if best is None or (res.success and res.fun < best.fun):
                    best = res

        # fall back to single start if all failed
        if best is None or not best.success:
            x0 = np.array([np.log(1.0), np.log(1.0)])
            best = optimize.minimize(
                _ztnb_nll, x0=x0, args=(k, hist),
                method="L-BFGS-B", bounds=bounds
            )

        log_rhat, log_muhat = best.x
        rhat = float(np.exp(log_rhat))
        muhat = float(np.exp(log_muhat))
        phat = float(rhat / (rhat + muhat))
        #phat = float(1.0 / (1.0 + np.exp(-log_phat)))

        pmf_raw = nbinom.pmf(k, rhat, phat)
        p0 = nbinom.pmf(0, rhat, phat)
        pred_pmf = pmf_raw / (1 - p0)
        pred_pmf = pred_pmf / pred_pmf.sum()

        true_pmf = hist / hist.sum()
        true_cdf = np.cumsum(true_pmf)
        pred_cdf = np.cumsum(pred_pmf)
        KS = float(np.max(np.abs(true_cdf - pred_cdf)))

        loglik = float(np.sum(hist * np.log(np.clip(pred_pmf, eps, 1.0))))
        nll_per_mol = float(-loglik / hist.sum())
        Fk = nbinom.cdf(int(k[-1]), rhat, phat)
        tail_mass = float((1.0 - Fk) / (1.0 - p0))

        flag = (rhat > (rmax - 1)) | (rhat < (rmin - 0.1))
        #if flag:
        #    print(f'[warning] ZT-NB fit found an anomalous rhat = {rhat:.5f} phat = {phat:.5f}', file=sys.stderr)
        self.rhat = rhat
        self.phat = phat
        self.KS = KS
        self.nll_per_mol = nll_per_mol
        self.tail_mass = tail_mass
        self.error_flag = bool(flag)

def full_fracs(fit : NBLibFit, full_summary : SimpleNamespace, 
               min_sat : float, max_sat : float, 
               use_full : bool, steps : int = 8) -> pd.DataFrame:
    sats = np.linspace(min_sat, max_sat, steps)
    reads = fit.reads_for_saturation(sats)
    fracs = reads / full_summary.countable_reads
    if use_full:
        fracs[-1] = 1.0
    return pd.DataFrame({'reads':reads, 'fraction':fracs, 'saturation':sats})

def saturation_fracs(fit : NBLibFit, full_summary : SimpleNamespace, sats, full : pd.Series | None = None) -> pd.DataFrame:
    reads = fit.reads_for_saturation(sats)
    fracs = reads / full_summary.countable_reads
    df = pd.DataFrame({'reads':reads, 'fraction':fracs, 'saturation':sats})
    if full is not None:
        df.loc[len(df)] = [full.reads, full.fraction, full.saturation]
    return df

def find_target_saturation(ds: Downsampler, full_summary: SimpleNamespace, 
                           target_sat: float = 50.0, max_sat : float = 55, 
                           seed: int = 42, threads: int = 1) -> tuple[pd.DataFrame, float]:

    ds.downsample([1.0], umi_len=full_summary.umi_length, 
            seed=seed, threads=threads, aggregate_only=True)
    full_stats = fn.bulk_stats(ds, full_summary)
    pd.set_option('display.float_format', '{:.2f}'.format)
    pd.set_option("styler.format.thousands", ",")
    ss = full_stats.iloc[0]
    if ss.saturation < max_sat:
        return full_stats, 1.0
    nbl = NBLibFit()
    nbl.fit(fn.get_rpm_hist(ds=ds)[0], reads=ss.reads, molecules=ss.molecules)
    frac = nbl.reads_for_saturation(target_sat) / full_summary.countable_reads
    ds.downsample([frac], umi_len=full_summary.umi_length, seed=seed, 
                  threads=threads, aggregate_only=True)

    df = fn.bulk_stats(ds, full_summary)
    b = df.iloc[0]
    nbl = NBLibFit()
    nbl.fit(fn.get_rpm_hist(ds=ds)[0], reads=b.reads, molecules=b.molecules)
    frac = nbl.reads_for_saturation(target_sat) / full_summary.countable_reads
    return full_stats, frac

def baseline_fitter(prefix : str, barcode_prefix : str = ''):
    if barcode_prefix == '':
        bdf = pd.read_csv(prefix + '_fit_baseline.txt', sep='\t')
        cdf = pd.read_csv(prefix + '_fit_curve.txt', sep='\t')
    else:
        bdf = pd.read_csv(f'{prefix}_{barcode_prefix}_fit_baseline.txt', sep='\t')
        cdf = pd.read_csv(f'{prefix}_{barcode_prefix}_fit_curve.txt', sep='\t')
    cdf = cdf[cdf.fraction >= 1].copy()
    nbl = NBLibFit()
    nbl.unserialize(bdf.iloc[0].to_dict())
    return nbl, bdf.iloc[0], cdf.iloc[0]

