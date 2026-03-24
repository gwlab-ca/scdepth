import numpy as np
from scipy import optimize
from scipy.stats import nbinom
from scipy.special import expit, logit, logsumexp

def _ztnb2_pred_pmf(k, pi, r1, r2, p, eps=1e-300):
    """
    Returns truncated mixture pmf over k (k>=1), normalized over full support (k>=1).
    """
    # log f(k) = log( pi * NB(k;r1,p) + (1-pi) * NB(k;r2,p) )
    log_fk = logsumexp(
        np.vstack([
            np.log(pi) + nbinom.logpmf(k, r1, p),
            np.log1p(-pi) + nbinom.logpmf(k, r2, p),
        ]),
        axis=0
    )

    # mixture P(X=0) = pi*p^r1 + (1-pi)*p^r2
    # compute stably:
    p0 = pi * (p ** r1) + (1.0 - pi) * (p ** r2)
    p0 = np.clip(p0, 0.0, 1.0 - 1e-16)

    # zero-truncate: log P(k|>0) = log f(k) - log(1-p0)
    log_trunc = log_fk - np.log1p(-p0)

    pmf = np.exp(log_trunc)
    pmf = np.clip(pmf, eps, np.inf)

    # In theory this already sums to 1 over k>=1.
    # Over a finite k grid (1..Kmax), renormalize like your current code.
    pmf = pmf / pmf.sum()
    return pmf

def _ztnb2_nll(theta, k, hist, rmin, rmax, pmin, pmax, eps=1e-12):
    """
    theta = [logit(pi), log(r1), log(delta), logit(p)]
    with r2 = r1 + exp(delta) (enforces r2 > r1)
    """
    lpi, lr1, ldelta, lp = theta
    pi = expit(lpi)
    r1 = np.exp(lr1)
    r2 = r1 + np.exp(ldelta)
    p  = expit(lp)

    # hard clamps to avoid wandering into numerically awful regions
    if not (rmin <= r1 <= rmax and rmin <= r2 <= rmax and pmin <= p <= pmax):
        return np.inf
    if pi < 1e-6 or pi > 1 - 1e-6:
        return np.inf

    pred = _ztnb2_pred_pmf(k, pi, r1, r2, p, eps=1e-300)
    true = hist / hist.sum()

    nll = -np.sum(hist * np.log(np.clip(pred, eps, 1.0)))
    return float(nll)

def fit_ztnb2(hist, rmin=0.25, rmax=50.0, pmin=1e-6, pmax=1-1e-6, eps=1e-12):
    hist = np.asarray(hist, dtype=float)
    k = np.arange(1, len(hist) + 1, dtype=int)

    # bounds in unconstrained space
    # lpi, lr1, ldelta, lp
    bounds = [
        (logit(1e-4), logit(1 - 1e-4)),           # pi
        (np.log(rmin), np.log(rmax)),             # r1
        (np.log(1e-6), np.log(rmax - rmin + 1)),  # delta (controls r2-r1)
        (logit(pmin), logit(pmax)),               # p
    ]

    # multi-start grid (tune as needed)
    starts_pi = [0.2, 0.5, 0.8]
    starts_r1 = [0.5, 1.0, 2.0, 5.0]
    starts_delta = [0.25, 1.0, 3.0]  # r2 roughly r1+delta
    starts_p = [0.1, 0.3, 0.5, 0.7]

    best = None
    for spi in starts_pi:
        for sr1 in starts_r1:
            for sdel in starts_delta:
                for sp in starts_p:
                    x0 = np.array([
                        logit(spi),
                        np.log(sr1),
                        np.log(max(sdel, 1e-6)),
                        logit(sp),
                    ], float)

                    res = optimize.minimize(
                        _ztnb2_nll, x0=x0,
                        args=(k, hist, rmin, rmax, pmin, pmax, eps),
                        method="L-BFGS-B",
                        bounds=bounds,
                    )
                    if best is None or (res.success and res.fun < best.fun):
                        best = res

    if best is None or not best.success:
        raise RuntimeError("2-component ZT-NB fit failed (all starts unsuccessful).")

    lpi, lr1, ldelta, lp = best.x
    pi = float(expit(lpi))
    r1 = float(np.exp(lr1))
    r2 = float(r1 + np.exp(ldelta))
    p  = float(expit(lp))

    pred_pmf = _ztnb2_pred_pmf(k, pi, r1, r2, p, eps=1e-300)

    # some diagnostics similar to your single-component version
    true_pmf = hist / hist.sum()
    true_cdf = np.cumsum(true_pmf)
    pred_cdf = np.cumsum(pred_pmf)
    KS = float(np.max(np.abs(true_cdf - pred_cdf)))
    nll_per_mol = float(best.fun / hist.sum())

    # mixture tail mass estimate over finite grid (like yours)
    # compute F(kmax | >0) and report remaining mass
    kmax = int(k[-1])
    # untruncated CDF at kmax:
    F1 = nbinom.cdf(kmax, r1, p)
    F2 = nbinom.cdf(kmax, r2, p)
    # mixture zero:
    p0 = pi * (p ** r1) + (1.0 - pi) * (p ** r2)
    # mixture truncated CDF at kmax:
    F_trunc = (pi * F1 + (1.0 - pi) * F2 - p0) / (1.0 - p0)
    tail_mass = float(max(0.0, 1.0 - F_trunc))

    return {
        "pi": pi,
        "r1": r1,
        "r2": r2,
        "p": p,
        "KS": KS,
        "nll_per_mol": nll_per_mol,
        "tail_mass": tail_mass,
        "opt": best,
        "pred_pmf": pred_pmf,
    }

def ztnb_component_contributions(k, pi, r1, r2, p):
    """
    Returns:
      total_trunc_pmf : P(k | X>0) for the mixture
      c1, c2          : component contributions in the truncated space
                        (they sum to total_trunc_pmf on the k-grid)
    """
    k = np.asarray(k, int)

    # raw NB pmfs
    nb1 = nbinom.pmf(k, r1, p)
    nb2 = nbinom.pmf(k, r2, p)

    # component zero masses
    p01 = p**r1
    p02 = p**r2

    # ZT component pmfs
    zt1 = nb1 / (1.0 - p01)
    zt2 = nb2 / (1.0 - p02)

    # mixture zero mass + observed mixture weight
    p0 = pi * p01 + (1.0 - pi) * p02
    pi_obs = (pi * (1.0 - p01)) / (1.0 - p0)

    # contributions in observed/truncated space
    c1 = pi_obs * zt1
    c2 = (1.0 - pi_obs) * zt2
    total = c1 + c2

    s = total.sum()
    c1, c2, total = c1/s, c2/s, total/s

    return total, c1, c2, pi_obs
