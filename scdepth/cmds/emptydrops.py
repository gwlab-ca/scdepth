#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import argparse
from scdepth.bindings import Downsampler
from scdepth.cmds import common
from scdepth import fn, filt
from typing import Optional, Sequence
import pandas as pd
import numpy as np
import scipy.sparse as sp

# Globals for lazy-loaded rpy2 stuff
ro = None
vectors = None
importr = None
localconverter = None
default_converter = None
pandas2ri = None
Matrix = None
DropletUtils = None

_RPY2_LOADED = False
_RPY2_ERROR = None

use_filter = True
use_exc_genes = True
use_vis_hd = False

def _load_rpy2():
    """
    Lazily import rpy2 and required R packages.
    Sets module-level globals so other functions can use them.
    Raises the original exception if something is missing.
    """
    global _RPY2_LOADED, _RPY2_ERROR
    global ro, vectors, importr, localconverter, default_converter, pandas2ri
    global Matrix, DropletUtils

    if _RPY2_LOADED:
        # We've already tried once; re-raise the error if there was one
        if _RPY2_ERROR is not None:
            raise _RPY2_ERROR
        return

    try:
        import rpy2.robjects as ro_
        from rpy2.robjects import vectors as vectors_
        from rpy2.robjects.packages import importr as importr_
        from rpy2.robjects.conversion import localconverter as localconverter_
        from rpy2.robjects import default_converter as default_converter_, pandas2ri as pandas2ri_

        Matrix_ = importr_("Matrix")
        DropletUtils_ = importr_("DropletUtils")

        # Assign to globals
        ro = ro_
        vectors = vectors_
        importr = importr_
        localconverter = localconverter_
        default_converter = default_converter_
        pandas2ri = pandas2ri_
        Matrix = Matrix_
        DropletUtils = DropletUtils_

        _RPY2_LOADED = True
        _RPY2_ERROR = None

    except Exception as e:
        _RPY2_LOADED = True
        _RPY2_ERROR = e
        raise

def custom_args(parser):
    g_init = parser.add_argument_group('Parameters')
    g_init.add_argument('-F', '--FDR', type=float, default=0.01,
                        help='Emptydrops FDR cutoff used for significant cell filtering')
    return parser

def resolve_args(args) -> dict:
    return {'FDR':args.FDR}

def build_parser(parser) -> argparse.ArgumentParser:
    common.prepare_args(parser, use_filter=use_filter, 
                      use_exc_genes=use_exc_genes, use_vis_hd=use_vis_hd,
                      command_args=custom_args)
    return parser

def run_emptydrops_csr(
    X: sp.csr_matrix,
    exclude_file : str,
    barcodes: Optional[Sequence[str]] = None,
    genes: Optional[Sequence[str]] = None,
    lower: int = 100,
    niters: int = 10000,
    retain: Optional[int] = None,
    test_ambient: bool = False,
    seed: Optional[int] = None,
) -> pd.DataFrame:

    try:
        _load_rpy2()
    except Exception as e:
        raise RuntimeError("EmptyDrops R dependencies missing: " + str(e))

    if not sp.isspmatrix_csr(X):
        raise TypeError("X must be a scipy.sparse.csr_matrix (barcodes x genes).")

    n_barcodes, n_genes = X.shape

    # Safety checks for names
    if barcodes is not None and len(barcodes) != n_barcodes:
        raise ValueError("Length of 'barcodes' does not match number of rows in X.")
    if genes is not None and len(genes) != n_genes:
        raise ValueError("Length of 'genes' does not match number of columns in X.")

    # Ensure integer counts
    X = X.astype(np.int64, copy=False)

    # Convert CSR -> COO for triplet representation
    X_coo = X.tocoo()

    # Build R vectors (NOTE: R is 1-based indexing)
    # Rows = barcodes; Cols = genes
    i = vectors.IntVector((X_coo.row + 1).astype(np.int32))
    j = vectors.IntVector((X_coo.col + 1).astype(np.int32))
    x = vectors.FloatVector(X_coo.data.astype(float))
    dims = vectors.IntVector([n_barcodes, n_genes])

    # Construct an R sparse matrix: barcodes x genes
    r_mat = Matrix.sparseMatrix(i=i, j=j, x=x, dims=dims)

    # Transpose to genes x barcodes for EmptyDrops
    r_counts = ro.r["t"](r_mat)

    # Optionally set dimnames (rownames = genes, colnames = barcodes)
    if genes is not None or barcodes is not None:
        # Fill missing names with generic strings if one side is None
        if genes is None:
            genes_r = [f"gene_{k}" for k in range(n_genes)]
        else:
            genes_r = list(map(str, genes))

        if barcodes is None:
            barcodes_r = [f"barcode_{k}" for k in range(n_barcodes)]
        else:
            barcodes_r = list(map(str, barcodes))

        r_gene_names = vectors.StrVector(genes_r)
        r_barcode_names = vectors.StrVector(barcodes_r)

        dimnames = ro.ListVector({"rownames": r_gene_names, "colnames": r_barcode_names})
        # r_counts is an S4 object; assign the Dimnames slot
        r_counts.do_slot_assign("Dimnames", dimnames)

    # Optional: set seed in R for reproducibility
    if seed is not None:
        ro.r["set.seed"](seed)

    # Build arguments for emptyDrops call
    edrops_kwargs = {
        "lower": lower,
        "niters": niters,
        "test.ambient": test_ambient,
    }
    if retain is not None:
        edrops_kwargs["retain"] = retain

    # Run emptyDrops on the genes x barcodes matrix
    r_res = DropletUtils.emptyDrops(r_counts, **edrops_kwargs)
    r_as_df = ro.r["as.data.frame"](r_res)


    with localconverter(default_converter + pandas2ri.converter):
        res_df: pd.DataFrame = ro.conversion.rpy2py(r_as_df)

    if barcodes is not None:
        res_df.index = pd.Index(barcodes, name="barcode")

    return res_df


def main(parser, args):
    try:
        import rpy2
    except ModuleNotFoundError:
        raise RuntimeError(
            "rpy2 is not installed. "
            "Please install scdepth with the 'emptydrops' extra:\n"
            "    pip install scdepth[emptydrops]\n"
            "and ensure R is installed with packages 'Matrix' and 'DropletUtils'."
        )

    args, _ = common.resolve_args(parser, args, resolve_args=resolve_args, 
                                        use_filter=use_filter, use_exc_genes=use_exc_genes, use_vis_hd=use_vis_hd)


    df = pd.read_csv(args.prefix + '_summary.txt', sep='\t')
    umi_length = int(df.loc[df['key'] == 'umi_length', 'value'].values[0])
    gdf = pd.read_csv(args.prefix + '_genes.txt.gz', sep='\t')

    ds = Downsampler()
    ds.init(args.prefix, max_hist=40, build_matrices=True) #, calc_usa=True)
    ds.downsample([1.0], umi_len=umi_length, seed=args.seed, threads=args.threads, aggregate_only=False,
            primer_mode=args.primer_mode)

    barcodes = np.array([b for b in ds.barcodes['barcode']])
    tot = ds.total_csr(0)
    gene_ids = gdf['gene_id'].values

    if args.samples is not None:
        print('FIX ME AND MAKE THIS AUTOMATIC')
        sample_cells = []
        for s in args.samples:
            idx = np.array([i for i, b in enumerate(barcodes) if b.startswith(f'{sample}_')])
            if len(idx) == 0:
                print(f'Warning: no barcodes found for {s}: {len(idx):,} barcodes')
                continue
            bcs = barcodes[idx]

            total_mols = np.asarray(tot.sum(axis=1)).ravel()
            stot = tot[idx]
            print(f'Calling emptydrops for {s}')
            cells = run_emptydrops_csr(X=stot, barcodes=bcs, genes=gene_ids, exclude_file=args.exclude_file)
            cells["total_mols"] = total_mols[idx]
            cells["is_cell"] = (
                cells["FDR"].notna()
                & (cells["FDR"] <= args.FDR)
            )

            sample_cells.append(cells)

        if not sample_cells:
            raise ValueError("No barcodes matched any values in args.samples")

        # Produces one DataFrame equivalent to single-sample mode.
        cells = pd.concat(sample_cells, axis=0, ignore_index=False)
    else:
        total_mols = np.asarray(tot.sum(axis=1)).ravel()

        print('Calling emptydrops')
        cells = run_emptydrops_csr(
            X=tot,
            barcodes=barcodes,
            genes=gene_ids,
            exclude_file=args.exclude_file,
        )

        cells["total_mols"] = total_mols
        cells["is_cell"] = (
            cells["FDR"].notna()
            & (cells["FDR"] <= args.FDR)
        )

    cells.to_csv(args.prefix + '_emptydrops_raw.txt.gz', sep='\t')

    barcodes = fn.read_barcodes_meta(ds, args.prefix, meta=args.prefix + '_emptydrops_raw.txt.gz')
    df = fn.barcode_df(ds=ds, df=barcodes.copy(), step=0)
    df, info = filt.filter_barcodes(df, add_passed=True, **args)

    df = df[['barcode','is_cell', 'passed']]

    ed = df['is_cell'].sum()
    pp = df['passed'].sum()
    print(f'Emptydrops Passed = {ed:,} [{ed/len(df) * 100.0:.2f}%]')
    print(f'Filtering Passed  = {pp:,} [{pp/len(df) * 100.0:.2f}%]')
    print(f'Total Barcodes    = {len(df):,}')
    #print(info)
    if args.barcode_prefix is None or args.barcode_prefix == '':
        df.to_csv(args.prefix + '_emptydrops.txt.gz', sep='\t', index=False)
    else:
        df.to_csv(f'{args.prefix}_{args.barcode_prefix}_emptydrops.txt.gz', sep='\t', index=False)
