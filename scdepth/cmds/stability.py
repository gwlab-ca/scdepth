#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

from scdepth.bindings import Downsampler, aggregate_visium_bins
from scdepth import fn, pl, fit, filt, libraries
from scdepth.cmds import common
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import warnings
from sklearn.metrics import adjusted_rand_score
import matplotlib.pyplot as plt
import json
from itertools import combinations

use_filter = True
use_exc_genes = True
use_vis_hd = True

def custom_args(parser):
    g_init = parser.add_argument_group('Parameters')
    g_init.add_argument('-r', '--repeats', type=int, default=10,
                        help='Number of curves to generate')
    g_init.add_argument('--recoveries', type=str,
                        help='Recovery targets for downsampling (default: 30,40,50,60,70)')
    g_init.add_argument('--resolution', type=str, default='0.6,0.8,1.0',
                        help='Leiden resolution to use for clustering')
    g_init.add_argument('--min-cells', type=int, default=10,
                        help='scanpy min_cell argument for filter_genes')
    g_init.add_argument('--hvgs', type=int, default=2000,
                        help='Number of highly variable genes')
    g_init.add_argument('--scale', type=int, default=10000,
                        help='Normalization scale')
    g_init.add_argument('--pcs', type=int, default=30,
                        help='Number of PCs')
    g_init.add_argument('--neighbours', type=int, default=30,
                        help='Number of neighbours')

    return parser

def build_parser(parser):
    common.prepare_args(parser, use_filter=use_filter, 
                      use_exc_genes=use_exc_genes, use_vis_hd=use_vis_hd,
                      command_args=custom_args)
    return parser

def resolve_args(args) -> dict:
    return {'repeats':args.repeats, 'recoveries':args.recoveries,
            'resolution':args.resolution, 'min_cells':args.min_cells, 
            'hvgs':args.hvgs, 'scale':args.scale, 'pcs':args.pcs,
            'neighbours':args.neighbours}

def build_adata(ds, step, kidx, obs, genes, bin_div = 0):
    if bin_div > 0:
        X_ref = ds.total_csr_bin(step, bin_div)
    else:
        X_ref = ds.total_csr(step)

    X_ref = X_ref[kidx]
    var = genes.copy()
    var.index = var['gene_id'].astype(str)
    return ad.AnnData(X=X_ref, obs=obs, var=var)
 
def run_seed_pipeline(adata, resolutions, min_cells, hvgs, scale, pcs, neighbours, **_):
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.highly_variable_genes(
        adata,
        flavor='seurat_v3',
        n_top_genes=hvgs,
        subset=False
    )

    hvg_mask = adata.var['highly_variable'].values

    sc.pp.normalize_total(adata, target_sum=scale)
    sc.pp.log1p(adata)
    adata = adata[:, hvg_mask].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=pcs)
    sc.pp.neighbors(adata, n_neighbors=neighbours, n_pcs=pcs)
    labels = []
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added='leiden', flavor='igraph', n_iterations=2)
        labels.append(adata.obs['leiden'].to_numpy())
    return adata, labels

def get_hvg_set(adata):
    return set(adata.var_names.tolist())

def get_knn_indices(adata, k=None):
    D = adata.obsp['distances'].tocsr()
    if k is None:
        k = int(adata.uns['neighbors']['params']['n_neighbors'])

    knn = []
    for i in range(D.shape[0]):
        cols = D.indices[D.indptr[i]:D.indptr[i+1]]
        vals = D.data[D.indptr[i]:D.indptr[i+1]]
        if cols.size == 0:
            knn.append(np.array([], dtype=np.int32))
            continue
        order = np.argsort(vals)[:k]   # smaller distance = closer
        knn.append(cols[order])
    return knn

def scalar_hvg_stability(hvgs_j):
    vals = hvgs_j.values
    iu = np.triu_indices_from(vals, k=1)
    return np.median(vals[iu])

def scalar_knn_stability(knn_pair_summary):
    return np.median(knn_pair_summary['median'].values)    


def jaccard_set(a, b):
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0

def compare_hvgs(hvg_sets, seed_names=None):
    S = len(hvg_sets)
    if seed_names is None:
        seed_names = [f'seed{i}' for i in range(S)]
    mat = np.zeros((S, S), dtype=float)
    for i in range(S):
        mat[i, i] = 1.0
        for j in range(i+1, S):
            mat[i, j] = mat[j, i] = jaccard_set(hvg_sets[i], hvg_sets[j])
    return pd.DataFrame(mat, index=seed_names, columns=seed_names)

def compare_knn_jaccard(knn_lists, seed_names=None):
    S = len(knn_lists)
    n_cells = len(knn_lists[0])
    if seed_names is None:
        seed_names = [f'seed{i}' for i in range(S)]

    # sanity: all have same cell count
    for s in range(1, S):
        assert len(knn_lists[s]) == n_cells, 'Cell count mismatch across seeds.'

    pair_cols = []
    per_cell_data = {}

    for (a, b) in combinations(range(S), 2):
        col = f'{seed_names[a]}__vs__{seed_names[b]}'
        pair_cols.append(col)

        ja = np.zeros(n_cells, dtype=float)
        for i in range(n_cells):
            Na = set(knn_lists[a][i])
            Nb = set(knn_lists[b][i])
            if not Na and not Nb:
                ja[i] = 1.0
            else:
                inter = len(Na & Nb)
                union = len(Na | Nb)
                ja[i] = inter / union if union else 0.0

        per_cell_data[col] = ja

    per_cell = pd.DataFrame(per_cell_data)
    summary = pd.DataFrame({
        'pair': pair_cols,
        'median': [np.median(per_cell[c].values) for c in pair_cols],
        'mean': [np.mean(per_cell[c].values) for c in pair_cols],
    }).set_index('pair')

    return per_cell, summary


def summarize_ari(labels_by_seed, res):
    faris = []
    lcounts = []

    NG = len(res)
    N = len(labels_by_seed)
    for z in range(NG):
        aris = []
        lcount = []
        for i in range(N):
            si = labels_by_seed[i][z]

            lcount.append(len(set(si)))
            for j in range(i+1, N):
                sj = labels_by_seed[j][z]
                ari = adjusted_rand_score(
                    si,
                    sj,
                )
                aris.append(ari)
        faris.append(aris)
        lcounts.append(lcount)

    faris = np.array(faris)
    lcounts = np.array(lcounts)

    return pd.DataFrame({
        'resolution':res,
        'ari_median':np.median(faris, axis=1),
        'ari_mean':np.mean(faris, axis=1),
        'labs_median':np.median(lcounts, axis=1),
        'labs_mean':np.mean(lcounts, axis=1),
    })

def summarize_upper_triangle(mat_df):
    vals = mat_df.values
    iu = np.triu_indices_from(vals, k=1)
    v = vals[iu]
    return dict(median=float(np.median(v)),
                mean=float(np.mean(v)))

def summarize_pair_df(pair_df, col='median'):
    v = pair_df[col].to_numpy()
    return dict(median=float(np.median(v)),
                mean=float(np.mean(v)))


def run_stability(args, full_summary):
    if args.barcode_prefix is None or args.barcode_prefix == '':
        oprefix = f'{args.prefix}_stability'
    else:
        oprefix = f'{args.prefix}_{args.barcode_prefix}_stability'

    ldata = libraries.library2ns(args.library)
    is_hd = 'visium_hd' in ldata.library

    ds = Downsampler()
    ds.init(args.prefix, max_hist=args.max_hist, build_matrices=True, exclude_file=args.exclude_file)
    barcodes = fn.read_barcodes_meta(ds, args.prefix, barcode_prefix=args.barcode_prefix)

    if is_hd:
        ds.init_visium(rows=barcodes['array_row'].to_numpy(), cols=barcodes['array_col'].to_numpy(),
                    in_tissue=barcodes['in_tissue'].to_numpy(), 
                    countable=barcodes['countable'].to_numpy(), 
                    total=barcodes['total'].to_numpy(), 
                    total_rows=args.total_rows, total_cols=args.total_cols)
    else:
        args.bin_div = 0


    nbl, _, _ = fit.baseline_fitter(args.prefix, barcode_prefix = args.barcode_prefix)

    if args.recoveries is None:
        args.recoveries = '30,40,50,60,70'

    recoveries = sorted(set(map(int, args.recoveries.split(','))))
    resolutions = sorted(set(map(float, args.resolution.split(','))))
    genes = pd.read_csv(f'{args.prefix}_genes.txt.gz', sep='\t')


    ds.downsample([1.0], umi_len=full_summary.umi_length, 
                seed=args.seed, threads=args.threads, aggregate_only=False)

    if is_hd:
        cdf, _ = aggregate_visium_bins(ds, 0, args.bin_div, args.bin_div)
        cdf = pd.DataFrame(cdf)
        fdf, meta = filt.filter_barcodes(cdf, **args)
        meta['bin_div'] = args.bin_div
        with open(oprefix + '_bin_filtering.txt', 'w') as f:
            json.dump(meta, f, indent=4)

    else:
        fdf = barcodes[barcodes['passed'] == 1].copy()
    if fdf is None:
        return

    kidx = np.array(fdf.index.values, dtype=np.uint32)

    rng = np.random.default_rng(args.seed)
    seeds = rng.choice(2**32, size=args.repeats, replace=False)
    fracs = nbl.reads_for_recovery(recoveries) / full_summary.countable_reads
    fracs = fracs[fracs <= 1.0]

    curves = []
    hvg_sets = [[] for _ in range(len(recoveries))]
    knn_lists = [[] for _ in range(len(recoveries))]
    label_lists = [[] for _ in range(len(recoveries))]

    warnings.filterwarnings(
        'ignore',
        message='zero-centering a sparse array/matrix densifies it.'
    )

    df = barcodes.loc[kidx, ['barcode']].copy()
    df['barcode'] = df['barcode'].astype(str)
    df.set_index('barcode', inplace=True, drop=True)

    for ci, s in enumerate(seeds):
        ds.downsample(fracs, umi_len=full_summary.umi_length,
                    seed=s, threads=args.threads, aggregate_only=False)
        curve_stats = fn.bulk_stats(ds, full_summary)
        curve_stats.insert(loc=0, column=f'recovery_target', value=recoveries)
        curve_stats.insert(loc=0, column='seed', value=s)
        curve_stats.insert(loc=0, column='curve_index', value=ci)
        curves.append(curve_stats)
        for j in range(len(fracs)):
            if is_hd:
                adata = build_adata(ds, j, kidx, df, genes, bin_div=args.bin_div)
            else:
                adata = build_adata(ds, j, kidx, df, genes)
            adata, labels = run_seed_pipeline(adata, resolutions, **args)
            hvg_sets[j].append(set(adata.var_names))
            knn_lists[j].append(get_knn_indices(adata))
            label_lists[j].append(labels)
            del adata

    x = []
    y = [[] for _ in resolutions]
    rows = []
    for i, r in enumerate(recoveries):
        hvgs_jaccard = compare_hvgs(hvg_sets[i], seed_names=seeds)
        _, knn_pair_summary = compare_knn_jaccard(knn_lists[i], seed_names=seeds)
        ari = summarize_ari(label_lists[i], resolutions)
        b = scalar_knn_stability(knn_pair_summary)

        hvg_sum = summarize_upper_triangle(hvgs_jaccard)
        knn_sum = summarize_pair_df(knn_pair_summary, col='median')

        row = {
            'sample': args.sample,
            'recovery': float(r),
            'hvg_jaccard_median': hvg_sum['median'],
            'hvg_jaccard_mean': hvg_sum['mean'],
            'knn_jaccard_median': knn_sum['median'],
            'knn_jaccard_mean': knn_sum['mean'],
        }

        for j, t in enumerate(ari.itertuples()):
            y[j].append(t.ari_median)    
            res = str(t.resolution).replace('.', 'p')  # 0.6 -> 0p6
            row[f'ari_median_res{res}'] = float(t.ari_median)
            row[f'n_clusters_median_res{res}'] = float(t.labs_median)
        x.append(b)
        rows.append(row)

    all_stats  = pd.concat(curves)
    all_stats.insert(loc=0, column='sample', value=args.sample)
    all_stats.to_csv(f'{oprefix}_curve_stats.txt', sep='\t', index=False)
    summary = pd.DataFrame(rows).sort_values(['sample', 'recovery'])
    summary.to_csv(f'{oprefix}_summary.txt', sep='\t', index=False)

def main(parser, args):
    args, summary = common.resolve_args(parser, args, resolve_args=resolve_args, 
                                    use_filter=use_filter, use_exc_genes=use_exc_genes, 
                                    use_vis_hd=use_vis_hd)

    run_stability(args, summary)
