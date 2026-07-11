#rSPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

from scdepth.bindings import Downsampler
from numpy.typing import NDArray
import numpy as np
import pandas as pd
from types import SimpleNamespace
import os

READ_CLASSES = ['spliced', 'ambiguous', 'unspliced', 'total']

def get_rpm_hist(ds : Downsampler, key: str = 'total', include_tail = False) -> NDArray[np.uint32]:
    ret = np.array([], dtype=np.uint32)
    if key not in READ_CLASSES:
        print(f'Invalid read class {key}')
        return ret
    if include_tail:
        return getattr(ds, f'{key}_mhist').copy()
    return getattr(ds, f'{key}_mhist')[:,:-1].copy()

def read_barcodes(prefix : str):
    barcodes = pd.read_csv(prefix + '_barcode_index.txt.gz', sep='\t')
    barcodes.set_index('barcode', inplace=True, drop=True)
    return barcodes

def read_barcodes_meta(ds : Downsampler, prefix : str, meta : str | None = None, barcode_prefix : str = '') -> pd.DataFrame:
    barcodes = pd.DataFrame(ds.barcodes)
    barcodes.set_index('barcode', inplace=True, drop=True)

    barcodes_pos = None
    if meta is not None:
        meta_file = meta
        barcodes_pos = pd.read_csv(meta_file, sep='\t')[['barcode','is_cell']]
        barcodes_pos['is_cell'] = barcodes_pos['is_cell'].astype(int)
    elif os.path.isfile(prefix + '_positions.parquet'):
        meta_file = prefix + '_positions.parquet'
        barcodes_pos = pd.read_parquet(meta_file)
    elif os.path.isfile(prefix + '_positions.csv'):
        meta_file = prefix + '_positions.csv'
        barcodes_pos = pd.read_csv(meta_file)
    elif barcode_prefix != '' and os.path.isfile(f'{prefix}_{barcode_prefix}_emptydrops.txt.gz'):
        meta_file = f'{prefix}_{barcode_prefix}_emptydrops.txt.gz'
        barcodes_pos = pd.read_csv(meta_file, sep='\t')[['barcode','is_cell', 'passed']]
        barcodes_pos['is_cell'] = barcodes_pos['is_cell'].astype(int)
    elif os.path.isfile(prefix + '_emptydrops.txt.gz'):
        meta_file = prefix + '_emptydrops.txt.gz'
        barcodes_pos = pd.read_csv(meta_file, sep='\t')[['barcode','is_cell', 'passed']]
        barcodes_pos['is_cell'] = barcodes_pos['is_cell'].astype(int)
    else:
        print('Could not find a barcode annotation file (scdepth_positions.parquet, scdepth_emptydrops.txt.gz or scdepth_positions.csv')
        return pd.DataFrame()

    barcodes_pos.set_index('barcode', inplace=True, drop=True)
    barcodes = barcodes.merge(barcodes_pos, left_index=True, right_index=True, how='left')
    barcodes.insert(0, 'barcode', barcodes.index, allow_duplicates=False)
    barcodes.reset_index(inplace=True, drop=True)
    del barcodes_pos
    return barcodes

def _get_bc_count(ds : Downsampler, rtype : str, key: str, 
                  squeeze : bool = True) -> NDArray[np.uint32]:

    ret = np.array([], dtype=np.uint32)
    if key not in READ_CLASSES:
        print(f'Invalid read class {key}')
        return ret

    ret = getattr(ds, f'{key}_bc_{rtype}').copy()
    if squeeze and ret.ndim == 2 and ret.shape[0] == 1:
        ret = ret.squeeze(axis=0)
    return ret

def _get_bc_basic(ds : Downsampler, key : str, squeeze : bool = True) -> NDArray[np.uint32]:
    ret = getattr(ds, f'bc_{key}')
    if squeeze and ret.ndim == 2 and ret.shape[0] == 1:
        ret = ret.squeeze(axis=0)
    return ret

def get_bc_reads(ds : Downsampler, key: str, 
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_count(ds, 'reads', key, squeeze=squeeze)

def get_bc_mols(ds : Downsampler, key: str,
                squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_count(ds, 'mols', key, squeeze=squeeze)

def get_bc_MT(ds : Downsampler, key: str,
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_count(ds, 'mt', key, squeeze=squeeze)

def get_bc_mod(ds : Downsampler, key: str,
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_count(ds, 'mod', key, squeeze=squeeze)

def get_bc_discarded_reads(ds : Downsampler,
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_basic(ds, 'discarded_reads', squeeze=squeeze)

def get_bc_excluded_reads(ds : Downsampler,
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_basic(ds, 'excluded_reads', squeeze=squeeze)

def get_bc_discarded_molecules(ds : Downsampler,
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_basic(ds, 'discarded_reads', squeeze=squeeze)

def get_bc_discarded_ambiguous(ds : Downsampler,
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_basic(ds, 'discarded_ambiguous', squeeze=squeeze)

def get_bc_merged_umis(ds : Downsampler,
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_basic(ds, 'merged_umis', squeeze=squeeze)

def get_bc_extra_umis(ds : Downsampler,
                 squeeze : bool = True) -> NDArray[np.uint32]:

    return _get_bc_basic(ds, 'extra_umis', squeeze=squeeze)

def get_bc_genes(ds : Downsampler, key: str,
                 squeeze : bool = True) -> NDArray[np.uint32]:
    return _get_bc_count(ds, 'genes', key, squeeze=squeeze)


def value_convert(x):
    try: 
        return int(x)
    except:
        try:
            return float(x)
        except:
            return x

def parse_summary(prefix : str) -> SimpleNamespace:
    df = pd.read_csv(prefix + '_summary.txt', sep='\t')
    d = {k:value_convert(v) for k, v in zip(df['key'].values, df['value'].values)}
    return SimpleNamespace(**d)

def bulk_stats(ds : Downsampler, summary : SimpleNamespace) -> pd.DataFrame:
    dd = {
        'fraction':ds.fracs.copy(), 
        'reads':ds.total_reads.copy(), 
        'molecules':ds.total_molecules.copy(), 
        'reads_discarded':ds.reads_discarded.copy(), 
        'reads_excluded':ds.reads_excluded.copy(), 
    }

    if ds.has_sau:
        dd['spliced_reads'] = ds.spliced_reads.copy()
        dd['spliced_molecules'] = ds.spliced_molecules.copy()
        dd['ambiguous_reads'] = ds.ambiguous_reads.copy()
        dd['ambiguous_molecules'] = ds.ambiguous_molecules.copy()
        dd['unspliced_reads'] = ds.unspliced_reads.copy()
        dd['unspliced_molecules'] = ds.unspliced_molecules.copy()
        dd['total_reads'] = ds.total_reads.copy()
        dd['total_molecules'] = ds.total_molecules.copy()

    df = pd.DataFrame(dd)
    disc = df['reads_discarded'] + df['reads_excluded']
    df['downsampled_frac'] = (disc + df['reads']) / summary.countable_reads
    df['saturation'] = 100.0 - 100 * df['molecules'] / df['reads']
    return df

def barcode_df(ds : Downsampler, df : pd.DataFrame, step: int,
               key : str = 'total') -> pd.DataFrame:

    df['reads'] = get_bc_reads(ds, key, squeeze=False)[step]
    df['molecules'] = get_bc_mols(ds, key, squeeze=False)[step]
    df['genes'] = get_bc_genes(ds, key, squeeze=False)[step]
    if ds.has_MT:
        df['mt_molecules'] = get_bc_MT(ds, key, squeeze=False)[step]
    if ds.has_module:
        df['mod_molecules'] = get_bc_mod(ds, key, squeeze=False)[step]
    if key == 'total':
        if ds.has_excluded:
            df['excluded_reads'] = get_bc_excluded_reads(ds, squeeze=False)[step]
        df['discarded_reads'] = get_bc_discarded_reads(ds, squeeze=False)[step]

    return df

