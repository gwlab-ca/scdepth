#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import argparse
from scdepth.cmds import subs
def main():
    parser = argparse.ArgumentParser(prog='scdepth')
    subparsers = parser.add_subparsers(help='scdepth command help', dest='command')

    for n, s, h in subs:
        p = subparsers.add_parser(n, help=h, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        s.build_parser(p)

    args = parser.parse_args()
    for n, s, h in subs:
        if n == args.command:
            s.main(parser, args)
            break
