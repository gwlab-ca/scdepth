#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import argparse
import pandas as pd

def build_parser(parser) -> argparse.ArgumentParser:
    g_init = parser.add_argument_group('Parameters')
    g_init.add_argument('-b', '--barcodes', type=str, required=True,
                        help='Barcodes file')
    parser.add_argument('prefix', help='scdepth cache prefix (also used as output prefix)')
    return parser

def main(_, args):
    if args.barcodes.endswith('parquet'):
        bc = pd.read_parquet(args.barcodes)
        bc[['barcode','array_row','array_col','in_tissue']].to_parquet(args.prefix + '_positions.parquet')
    elif args.barcodes.endswith('csv'):
        bc = pd.read_csv(args.barcodes)
        bc[['barcode','array_row','array_col','in_tissue']].to_csv(args.prefix + '_positions.csv', index=False)
    else:
        print('Unrecognized barcode format for ', args.barcodes)

