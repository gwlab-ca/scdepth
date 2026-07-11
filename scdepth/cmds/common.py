#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import argparse, os
from typing import Any
from scdepth import fn, libraries
from pathlib import Path
from collections.abc import Callable
from collections.abc import MutableMapping

class Args(MutableMapping):
    def __init__(self, **data):
        self.__dict__.update(data)

    # ---- dict interface ----
    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __delitem__(self, key):
        del self.__dict__[key]

    def __iter__(self):
        return iter(self.__dict__)

    def __len__(self):
        return len(self.__dict__)

    # ---- nice repr ----
    def __repr__(self):
        items = ", ".join(f"{k}={v!r}" for k, v in self.items())
        return f"Args({items})"

AddArgsFn = Callable[[argparse.ArgumentParser], None]
ResolveArgsFn = Callable[[argparse.ArgumentParser], dict]

#def add_mt_gene_args(parser: argparse.ArgumentParser):
#    grp = parser.add_argument_group('Barcode filtering')
#    grp.add_argument(
#            '-m', '--mt-prefix', 
#            type=str, 
#            default='',
#            help='Mitochondria gene name prefix'
#    )
#
#    grp.add_argument(
#            '-M', '--mt-file', 
#            type=str, 
#            default='',
#            help='Mitochondria gene ID file (one gene per line no header)'
#    )
#
#    grp.add_argument(
#            '-X', '--mod-file', 
#            type=str, 
#            default='',
#            help='Additional module gene ID file (one gene per line no header)'
#    )

#def resolve_mt_gene_args(args) -> dict:
#    return {
#        'mt_file': args.mt_file,
#        'mt_prefix': args.mt_prefix,
#        'mod_file': args.mod_file,
#    }

def add_barcode_filter_args(parser: argparse.ArgumentParser):
    grp = parser.add_argument_group('Barcode filtering')

    grp.add_argument(
        '--min-molecules',
        type=int,
        default=None,
        help=(
            'Baseline minimum molecules before MAD filtering. '
            'If not set, a technology-specific default is used '
            '(scRNA/Visium: 500, Visium HD = 10 for 1x1, 20 for 2x2/3x3, 50 for 4x4, otherwise 100'
        ),
    )

    grp.add_argument(
        '--molecule-mads',
        type=float,
        help=(
            'Number of MADs below the median (on log10(molecules)) '
            'used to define the adaptive minimum read cutoff. Default 3.0 for scRNA/Visium 4.0 for Visium HD'
        ),
    )

    #grp.add_argument(
    #    '--mt-mads',
    #    type=float,
    #    default=4.0,
    #    help=(
    #        'Upper-tail MAD multiplier for mitochondrial fraction filtering '
    #        '(median + k*MAD). Set to 0 to disable MAD-based MT filtering.'
    #        ' Not used for visium HD'
    #    ),
    #)

    #grp.add_argument(
    #    '--mt-cap',
    #    type=float,
    #    default=0.30,
    #    help=(
    #        'Absolute maximum allowed mitochondrial fraction. '
    #        'Final MT cutoff is min(mt_cap, median + k*MAD).'
    #        ' Not used for visium HD'
    #    ),
    #)

def resolve_filter_args(args) -> dict:
    # --- MT parameters ---
    #if args.mt_file == '' and args.mt_prefix == '':
    #    mt_mads = None
    #    mt_cap = None
    #else:
    #    mt_mads = args.mt_mads
    #    mt_cap = args.mt_cap

    return {
        'min_molecules': args.min_molecules,
        'molecule_mads': args.molecule_mads,
        #'mt_mads': mt_mads,
        #'mt_cap': mt_cap,
        #'mt_file': args.mt_file,
        #'mt_prefix': args.mt_prefix,
        'total_rows':args.vis_rows,
        'total_cols':args.vis_cols,
        'bin_div':args.vis_bins,
    }

def add_visium_hd_args(parser: argparse.ArgumentParser):
    grp = parser.add_argument_group('Barcode filtering / efficiency stability')
    grp.add_argument(
        '--vis-rows',
        type=int,
        default=3350,
        help=(
            'Number of rows in the visium grid'
        ),
    )

    grp.add_argument(
        '--vis-cols',
        type=int,
        default=3350,
        help=(
            'Number of columns in the visium grid'
        ),
    )

    grp.add_argument(
        '--vis-bins',
        type=int,
        default=8,
        help=(
            'Bin divisor for rows/cols for standard analyses such as efficiency and plotting (ie BxB)'
        ),
    )

    return grp

def resolve_vis_hd_args(args) -> dict:
    return {
        'total_rows':args.vis_rows,
        'total_cols':args.vis_cols,
        'bin_div':args.vis_bins,
    }

def add_exc_gene_args(parser: argparse.ArgumentParser):
    grp = parser.add_argument_group('Gene filtering')
    grp.add_argument(
        '--no-probes',
        action='store_true',
        help=(
            'If a scdepth_probes.csv file exists do not use it and do not warn if its missing for the library'
        ),
    )

    grp.add_argument(
        '--exclude-genes',
        type=str,
        help=(
            'Exclude gene file to use instead of the probe file. Important. It is up to the user to'
            ' include the filtered probe based gene_ids'
        ),
    )

    return grp

def resolve_exc_gene_args(args, has_probes: bool):
    if args.exclude_genes:
        return {'exclude_file': args.exclude_genes}

    if not has_probes:
        return {'exclude_file':''}

    if args.no_probes:
        print('[warning] Library is probe-based but --no-probes was set and no --exclude-genes was provided')
        return {'exclude_file':''}

    ff = args.prefix + '_exclude.txt'
    if not os.path.isfile(ff):
        print(f'[warning] Library is probe-based but scdepth_exclude.txt file not found: {ff}. Provide --exclude-genes or use `scdepth probes`.')
        return {'exclude_file':''}

    return {'exclude_file': ff}

def prepare_args(parser: argparse.ArgumentParser, use_exc_genes : bool, use_filter : bool,
                 use_vis_hd : bool, command_args: AddArgsFn | None = None):
    g_common = parser.add_argument_group('Common')
    g_common.add_argument('-s', '--sample', type=str, default='sample',
                        help='Sample name to include in output')
    g_common.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of downsampling threads')
    g_common.add_argument('--max-hist', type=int, default=50,
                        help='Maximum bin for the reads per molecule histogram')
    g_common.add_argument('--custom-libs', type=Path, help='Path to JSON file defining custom libraries '
                          '(must be specified if a custom library was used when caching tags)')
    g_common.add_argument('--barcode-prefix', type=str, help='Only process scRNA library barcodes starting with this prefix, '
                          'useful to filter specific samples from multiplexed libraries such as Parse WT data',
                          default='')
    g_common.add_argument('--primer-type', type=str, help='For scRNA libaries such as Parse WT this permits selecting a primer type or merging both: '
                          'options are merge/polyA/random_hex', choices=['merge','polyA','random_hex'], default='merge')
    g_common.add_argument('-S', '--seed', type=int, default=42,
                        help='Seed used for downsampling')

    if command_args is not None:
        command_args(parser)
    #if use_MT:
    #    add_mt_gene_args(parser)
    if use_exc_genes:
        add_exc_gene_args(parser)
    if use_filter:
        add_barcode_filter_args(parser)
    if use_vis_hd:
        add_visium_hd_args(parser)

    parser.add_argument('prefix', help='scdepth cache prefix (also used as output prefix)')

def resolve_args(parser, args, use_exc_genes : bool, use_filter : bool,
                 use_vis_hd : bool, resolve_args: ResolveArgsFn | None = None) -> tuple[Any, Any]:
    params = {
        'prefix':args.prefix,
        'sample':args.sample,
        'threads':args.threads,
        'max_hist':args.max_hist,
        'seed':args.seed,
    }
    summary = fn.parse_summary(args.prefix)
    if args.custom_libs:
        if not args.custom_libs.is_file():
            raise SystemExit(f'Custom libs file not found: {args.custom_libs}')
        libraries.add_custom_libraries(args.custom_libs)
    if not libraries.is_library(summary.library_string):
        parser.error(
            f'Sample summary library "{summary.library_string}" was not found '
            f'in default or custom libraries. If it is custom, pass --custom-libs.'
        )
    params['library'] = summary.library_string

    #if use_MT:
    #    params.update(resolve_mt_gene_args(args))
    if use_exc_genes:
        params.update(resolve_exc_gene_args(args, libraries.is_probe_based(summary.library_string)))
    if use_filter:
        params.update(resolve_filter_args(args))
    if use_vis_hd:
        params.update(resolve_vis_hd_args(args))

    if resolve_args is not None:
        params.update(resolve_args(args))

    return Args(**params), summary

def interactive_args(
                  prefix: str,
                  overrides: dict,
                  use_exc_genes = False,
                  use_filter = False, use_vis_hd = False):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    prepare_args(parser, use_exc_genes=use_exc_genes, use_filter=use_filter, use_vis_hd=use_vis_hd)
    ns = parser.parse_args([prefix])

    for k, v in overrides.items():
        setattr(ns, k, v)

    return resolve_args(parser, ns, use_exc_genes=use_exc_genes, 
                        use_filter=use_filter, use_vis_hd=use_vis_hd)
