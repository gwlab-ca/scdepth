#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

from scdepth.cmds import counter, libraries, emptydrops, curve, \
        limits, probes, barcodes, genes, preseq, stability

subs = (
    ('libraries', libraries, 'List pre-defined libraries'),
    ('cache', counter, 'Cache molecule tags for downsampling etc'),
    ('barcodes', barcodes, 'Build an annotated barcode list for Visium and Visium HD samples from Space Ranger output'),
    ('probes', probes, 'Build a gene exclusion file from the 10X probe set file'),
    ('emptydrops', emptydrops, 'Run emptydrops using rpy2'),
    ('fit', curve, 'Fit the full curve for a cached sample'),
    ('limit', limits, 'Limited saturation analyses'),
    ('preseq', preseq, 'Run preseq on limited saturation analyses'),
    ('genes', genes, 'Gene overlap/stability analyses'),
    ('stability', stability, 'Barcode stability analyses'),
)

