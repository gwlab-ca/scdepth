#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import argparse
import pandas as pd

def build_parser(parser) -> argparse.ArgumentParser:
    g_init = parser.add_argument_group('Parameters')
    g_init.add_argument('-p', '--probe', type=str, required=True,
                        help='Probe set file')
    parser.add_argument('prefix', help='scdepth cache prefix (also used as output prefix)')
    return parser

def main(parser, args):
    df = pd.read_csv(args.probe, comment='#')
    #print(df.head())
    #print(df.included.sum(), len(df))
    exc = df[df.included == False]

    with open(f'{args.prefix}_exclude.txt', 'w') as fp:
        eset = set()
        for t in exc.itertuples():
            gg = t.probe_id.split('|')[0]
            if gg.startswith('DEPRECATED_'):
                eset.add(gg[len('DEPRECATED_'):])
            else:
                eset.add(gg)
        print(f'Writing {len(eset)} excluded genes')
        for g in sorted(eset):
            fp.write(f'{g}\n')

