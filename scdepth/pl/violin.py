import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe


def _density_on_grid(x: np.ndarray, grid: np.ndarray, bandwidth: float | None = None):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.zeros_like(grid, dtype=float), np.nan

    grid = np.asarray(grid, dtype=float)
    if grid.ndim != 1 or grid.size < 2:
        raise ValueError('grid must be a 1D array with at least 2 points')
    if not np.all(np.diff(grid) > 0):
        raise ValueError('grid must be strictly increasing')

    dy = float(np.median(np.diff(grid)))
    if dy <= 0:
        raise ValueError('grid spacing must be positive')

    # Build bin edges from grid centers
    edges = np.concatenate(([grid[0] - dy / 2], grid + dy / 2))
    counts, _ = np.histogram(x, bins=edges)

    # Bandwidth in data units
    if bandwidth is None:
        n = x.size
        s = float(np.std(x, ddof=1)) if n > 1 else dy
        bandwidth = s * (n ** (-1 / 5)) if (n > 1 and s > 0) else dy

    sigma_bins = max(float(bandwidth) / dy, 1e-6)

    # Smooth histogram counts (mode controls boundary handling)
    smooth = gaussian_filter1d(counts.astype(float), sigma=sigma_bins, mode='nearest')

    # Normalize to density (area under curve = 1)
    density = smooth / (smooth.sum() * dy + 1e-300)
    return density, float(bandwidth)

def plot_violin(ax, v, pos=0.0, width=0.9, color='C0', alpha=0.9,
    edgecolor='k', linewidth=1.0, nbins=400, ylims=None, bandwidth=None, max_points_kde=200_000, median_size=36,
    zorder=2, median_color='k'):
    nbins = int(nbins)

    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]

    def _maybe_subsample(x):
        if x.size > max_points_kde:
            rng = np.random.default_rng(0)
            return x[rng.choice(x.size, max_points_kde, replace=False)]
        return x

    vs = _maybe_subsample(v)

    if ylims is None:
        ymins, ymaxs = float(np.min(v)), float(np.max(v))
        pad = 0.05 * (ymaxs - ymins + 1e-12)
        ymins -= pad
        ymaxs += pad
    else:
        ymins, ymaxs = ylims

    y = np.linspace(ymins, ymaxs, nbins)

    d1, _ = _density_on_grid(vs, y, bandwidth=bandwidth)

    scale = (width / 2) / max(d1.max(), 1e-12)
    xl = pos - d1 * scale
    xr = pos + d1 * scale

    ax.fill_betweenx(y, xl, xr, facecolor=color, alpha=alpha, edgecolor=edgecolor, 
                     linewidth=linewidth, zorder=zorder)

    med = np.percentile(v, 50)
    ax.scatter([pos], [med], s=median_size, marker='o',
        color=median_color, edgecolors='w', linewidths=1.2,
        zorder=zorder + 10,
    )
