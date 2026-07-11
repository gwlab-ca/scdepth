#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import scdepth.bindings
import tqdm, argparse
import pandas as pd
from scdepth import libraries
from pathlib import Path

def build_parser(parser) -> argparse.ArgumentParser:
    g_init = parser.add_argument_group('Required')
    g_init.add_argument('-g', '--gtf', required=True, type=str,
                        help='Annotation GTF file')
    g_init.add_argument('--custom-libs', type=Path, help='Path to JSON file defining custom libraries '
                          ' must be specified if a custom library is used')
    g_init.add_argument('-l', '--lib-type', required=True,
                        help='Library identifier can from scdepth libraries for a list or a custom name')
    g_init.add_argument('--sample-tag', required=False, default="",
                        help='Sample identifier tag for multiplex libraries such as parse')
    g_init.add_argument('--samples', required=True,
                        help='Comma separated list of valid samples')
    #g_init.add_argument('--benchmark', action='store_true',
    #                    help='Calculate timing / memory usage')
    g_init.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress output such as tqdm')

    g_filter = parser.add_argument_group('Filtering (set_count_parameters)')
    g_filter.add_argument('--min-gene', type=float, default=0.95,
                          help='Minimum fraction of gene overlap')
    g_filter.add_argument('--min-gene-bases', type=int, default=40,
                          help='Minimum number of bases overlapping a gene')
    g_filter.add_argument('--min-exonic', type=float, default=0.95,
                          help='Minimum exonic fraction')
    g_filter.add_argument('--min-intronic', type=int, default=15,
                          help='Minimum intronic overlap')
    g_filter.add_argument('--min-qual', type=int, default=255,
                          help='Minimum read quality score')
    g_filter.add_argument('--discard-unknown-juncs', action='store_true',
                          help='Discard reads with unknown splice junctions')

    g_prep = parser.add_argument_group('BAM preparation (prepare_bam)')
    g_prep.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of processing threads')
    g_prep.add_argument('--max-tags', type=int, default=200000000,
                        help='Approximate maximum number of tags to keep in memory')
    g_prep.add_argument('--max-tag-frac', type=float, default=0.95,
                        help='Maximum tag fraction')
    g_prep.add_argument('--chunk', type=int, default=1000000,
                        help='Number of reads per processing chunk')


    parser.add_argument('bam', help='Input BAM file')
    parser.add_argument('out', help='Output Prefix')

    return parser

def process_bam(
    args, gtf: str, bam: str, out: str, lib_type : str, threads: int = 1, 
    max_tags: int = 200000000, max_tag_frac: float = 0.95, 
    chunk: int = 1000000, min_gene: float = 0.95, min_gene_bases: int = 40, 
    min_exonic: float = 0.95, min_intronic: int = 15, min_qual: int = 255, 
    discard_unknown_juncs: bool = False,  quiet: bool = False, samples = "", sample_tag = "",
    **_
    ):


    res = libraries.library2ns(lib_type)

          #'random_hex_re': '',
          #'random_hex_value': '',
    bc = scdepth.bindings.BarcodeCounter()
    bc.init(lib_type, fwd=(not res.five_prime_like), barcode_tag=res.CB_tag, barcode_re=res.CB_re, umi_tag=res.UR_tag, 
            barcode_length = res.CB_length, umi_length = res.UR_length, samples=samples, sample_tag=sample_tag,
            random_hex_re=res.random_hex_re, random_hex_value=res.random_hex_value)
    if not quiet and res.probe_based:
        print('Using probe based gene ids')

    bc.set_count_parameters(min_gene=min_gene, min_gene_bases=min_gene_bases, 
                            min_exonic=min_exonic, min_intronic=min_intronic, 
                            min_qual=min_qual, discard_unknown_juncs=discard_unknown_juncs, 
                            probes=res.probe_based)

    if not quiet:
        print('Parsing GTF and building indices')
    bc.prepare_bam(gtf=gtf, bam=bam, out=out, threads=threads, max_tags=max_tags, 
                    max_tag_frac=max_tag_frac);

    pbar = tqdm.tqdm(total = 0, disable=quiet, unit=' reads', desc='Overall Process')
    while not bc.done():
        processed = bc.process_reads(chunk)
        pbar.set_postfix({'Countable Rate':f'{100.0 * bc.countable_reads()/bc.total_reads():.2f}%'})
        pbar.update(processed)
    pbar.close()
    print('Merging in memory + shard files and writing the output')
    bc.finish()    

def main(parser, args):
    kwargs = vars(args)
    if args.custom_libs:
        if not args.custom_libs.is_file():
            raise SystemExit(f'Custom libs file not found: {args.custom_libs}')
        libraries.add_custom_libraries(args.custom_libs)

    if not libraries.is_library(args.lib_type):
        parser.error(
            f'Library "{args.lib_type}" was not found '
            f'in default or custom libraries. If it is custom, pass --custom-libs.')
    #if 'benchmark' in kwargs:
    #    del kwargs['benchmark']
    #    _, res = scdepth.benchmark.run(process_bam, **kwargs)
    #    with open(kwargs['out'] + '_benchmark.txt', 'w') as fo:
    #        fo.write('key\tvalue\n') 
    #        for k, v in res.items():
    #            fo.write(f'{k}\t{v}\n')
    #else:
    #    process_bam(**kwargs)
    process_bam(args, **kwargs)

