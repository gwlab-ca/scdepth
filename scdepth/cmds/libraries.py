#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import argparse
from scdepth.libraries import print_libs

def build_parser(parser) -> argparse.ArgumentParser:
    return parser

def main(*_):
    print_libs()

