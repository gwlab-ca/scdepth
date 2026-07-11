"""
scdepth.bindings
"""
from __future__ import annotations
import collections.abc
import numpy
import numpy.typing
import typing
__all__: list[str] = ['BarcodeCounter', 'Downsampler', 'aggregate_visium_bins', 'calculate_gene_overlaps']
class BarcodeCounter:
    """
    Caches read data from a scRNA-seq sam/bam/cram file
    """
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def countable_reads(self) -> int:
        """
        Countable Reads Processed
        """
    def done(self) -> bool:
        """
        True when processing is complete
        """
    def finish(self) -> bool:
        """
        Finalize and flush outputs
        """
    def init(self, lib_string: str, fwd: bool, barcode_tag: str, barcode_re: str, umi_tag: str, sample_tag: str = '', samples: collections.abc.Sequence[str] = [], random_hex_re: str = '', random_hex_value: str = '', barcode_length: typing.SupportsInt | typing.SupportsIndex = 0, umi_length: typing.SupportsInt | typing.SupportsIndex = 0) -> None:
        """
        Configure the scRNA-seq library properties
        """
    def prepare_bam(self, gtf: str, bam: str, out: str, threads: typing.SupportsInt | typing.SupportsIndex = 1, max_tags: typing.SupportsInt | typing.SupportsIndex = 200000000, max_tag_frac: typing.SupportsFloat | typing.SupportsIndex = 0.95) -> bool:
        """
        Prepare the bam for processing
        """
    def process_reads(self, chunk: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
        Process up to `chunk` reads and return the number processed
        """
    def set_count_parameters(self, min_gene: typing.SupportsFloat | typing.SupportsIndex = 0.95, min_gene_bases: typing.SupportsInt | typing.SupportsIndex = 40, min_exonic: typing.SupportsFloat | typing.SupportsIndex = 0.95, min_intronic: typing.SupportsInt | typing.SupportsIndex = 15, min_qual: typing.SupportsInt | typing.SupportsIndex = 255, discard_unknown_juncs: bool = False, probes: bool = False) -> None:
        """
        Configure read filtering
        """
    def total_reads(self) -> int:
        """
        Total Reads Processed
        """
class Downsampler:
    def __init__(self) -> None:
        ...
    def ambiguous_csr(self, step: typing.SupportsInt | typing.SupportsIndex) -> typing.Any:
        ...
    def clear_output(self) -> None:
        """
        Clear output memory
        """
    def downsample(self, fracs: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], umi_len: typing.SupportsInt | typing.SupportsIndex, seed: typing.SupportsInt | typing.SupportsIndex, threads: typing.SupportsInt | typing.SupportsIndex = 1, aggregate_only: bool = False, umi_mode: str = 'directed', correct_multi_umis: bool = True, primer_mode: str = 'merge') -> bool:
        """
        Run downsampling for the given barcodes (optional) and fractions
        """
    def init(self, prefix: str, mt_prefix: str = '', mt_file: str = '', mod_file: str = '', exclude_file: str = '', barcode_prefix: str = '', max_hist: typing.SupportsInt | typing.SupportsIndex = 50, build_matrices: bool = False, calc_sau: bool = False) -> bool:
        """
        Initialize the downsampler from a cached prefix
        """
    def init_visium(self, rows: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], cols: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], in_tissue: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], countable: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], total: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], total_rows: typing.SupportsInt | typing.SupportsIndex, total_cols: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Initialize the downsampler to create binned tiles/barcodes for visium HD data
        """
    def reset_visium(self) -> bool:
        """
        Reset visium tile aggregation information
        """
    def spliced_csr(self, step: typing.SupportsInt | typing.SupportsIndex) -> typing.Any:
        ...
    def total_csr(self, step: typing.SupportsInt | typing.SupportsIndex) -> typing.Any:
        ...
    def total_csr_bin(self, step: typing.SupportsInt | typing.SupportsIndex, bin_div: typing.SupportsInt | typing.SupportsIndex) -> typing.Any:
        ...
    def unspliced_csr(self, step: typing.SupportsInt | typing.SupportsIndex) -> typing.Any:
        ...
    def write_gene_baseline(self, output: str, idx: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], step: typing.SupportsInt | typing.SupportsIndex = 0, bin_div: typing.SupportsInt | typing.SupportsIndex = 0) -> bool:
        """
        Write the gene detection full count matrix to a file for comparisons (optionally bin for visium HD data)
        """
    def write_gene_mats(self, output: str, idx: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], bin_div: typing.SupportsInt | typing.SupportsIndex = 0) -> bool:
        """
        Write the gene detection matrices to a file for comparisons (optionally bin for visium HD data)
        """
    @property
    def aggregate_only(self) -> bool:
        """
        Whether sparse matrices were built
        """
    @property
    def ambiguous_bc_genes(self) -> typing.Any:
        """
        Ambiguous genes per barcode shape (steps x barcodes)
        """
    @property
    def ambiguous_bc_mod(self) -> typing.Any:
        """
        Ambiguous module mols per barcode shape (steps x barcodes)
        """
    @property
    def ambiguous_bc_mols(self) -> typing.Any:
        """
        Ambiguous molecules per barcode shape (steps x barcodes)
        """
    @property
    def ambiguous_bc_mt(self) -> typing.Any:
        """
        Ambiguous MT mols per barcode shape (steps x barcodes)
        """
    @property
    def ambiguous_bc_reads(self) -> typing.Any:
        """
        Ambiguous reads per barcode shape (steps x barcodes)
        """
    @property
    def ambiguous_mhist(self) -> numpy.typing.NDArray[numpy.uint32]:
        """
        Histogram of ambiguous reads per MT molecule (steps x max_hist)
        """
    @property
    def ambiguous_molecules(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of ambiguous molecules
        """
    @property
    def ambiguous_reads(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of ambiguous reads
        """
    @property
    def barcode_prefix(self) -> str:
        """
        Number of downsampling fractions
        """
    @property
    def barcodes(self) -> dict:
        """
        Barcode columns as a dict-of-arrays (DataFrame-friendly)
        """
    @property
    def bc_discarded_moleculess(self) -> typing.Any:
        """
        Molecules lost from multi-gene UMIs (steps x barcodes)
        """
    @property
    def bc_discarded_reads(self) -> typing.Any:
        """
        Reads lost from multi-gene UMIs and excluded gene lists (steps x barcodes)
        """
    @property
    def bc_excluded_reads(self) -> typing.Any:
        """
        Reads lost from excluded read lists (steps x barcodes)
        """
    @property
    def build_mats(self) -> bool:
        """
        Whether sparse matrices were built
        """
    @property
    def fracs(self) -> list[float]:
        """
        Downsampling fractions actually used
        """
    @property
    def has_MT(self) -> bool:
        """
        A MT gene list was included
        """
    @property
    def has_excluded(self) -> bool:
        """
        An excluded gene list was included
        """
    @property
    def has_module(self) -> bool:
        """
        A module gene list was included
        """
    @property
    def has_sau(self) -> bool:
        """
        Has spliced/ambig/unspliced quantifications
        """
    @property
    def max_hist(self) -> int:
        """
        Maximum histogram bin (histograms have size max_hist+1)
        """
    @property
    def molecules_ambig(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of molecule subgraphs mapping to at least 2 genes.
        """
    @property
    def molecules_discarded(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of molecules lost to ambiguous UMI/gene mappings
        """
    @property
    def primer_mode(self) -> str:
        """
        Number of downsampling fractions
        """
    @property
    def reads_discarded(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of reads lost to ambiguous UMI/gene mappings
        """
    @property
    def reads_excluded(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of reads lost to excluded gene filter
        """
    @property
    def seed(self) -> int:
        """
        Downsampling seed
        """
    @property
    def spliced_bc_genes(self) -> typing.Any:
        """
        Spliced genes per barcode shape (steps x barcodes)
        """
    @property
    def spliced_bc_mod(self) -> typing.Any:
        """
        Spliced module mols per barcode shape (steps x barcodes)
        """
    @property
    def spliced_bc_mols(self) -> typing.Any:
        """
        Spliced molecules per barcode shape (steps x barcodes)
        """
    @property
    def spliced_bc_mt(self) -> typing.Any:
        """
        Spliced MT mols per barcode shape (steps x barcodes)
        """
    @property
    def spliced_bc_reads(self) -> typing.Any:
        """
        Spliced reads per barcode shape (steps x barcodes)
        """
    @property
    def spliced_mhist(self) -> numpy.typing.NDArray[numpy.uint32]:
        """
        Histogram of spliced reads per molecule (steps x max_hist)
        """
    @property
    def spliced_molecules(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of spliced molecules
        """
    @property
    def spliced_reads(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of spliced reads
        """
    @property
    def total_barcodes(self) -> int:
        """
        Number of raw barcodes in the results
        """
    @property
    def total_bc_genes(self) -> typing.Any:
        """
        Total genes per barcode shape (steps x barcodes)
        """
    @property
    def total_bc_mod(self) -> typing.Any:
        """
        Total module mols per barcode shape (steps x barcodes)
        """
    @property
    def total_bc_mols(self) -> typing.Any:
        """
        Total molecules per barcode shape (steps x barcodes)
        """
    @property
    def total_bc_mt(self) -> typing.Any:
        """
        Total MT mols per barcode shape (steps x barcodes)
        """
    @property
    def total_bc_reads(self) -> typing.Any:
        """
        Total reads per barcode shape (steps x barcodes)
        """
    @property
    def total_fracs(self) -> int:
        """
        Number of downsampling fractions
        """
    @property
    def total_genes(self) -> int:
        """
        Number of genes in the results
        """
    @property
    def total_mhist(self) -> numpy.typing.NDArray[numpy.uint32]:
        """
        Histogram of total reads per molecule (steps x max_hist)
        """
    @property
    def total_molecules(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of total molecules
        """
    @property
    def total_reads(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of total reads
        """
    @property
    def unspliced_bc_genes(self) -> typing.Any:
        """
        Unspliced genes per barcode shape (steps x barcodes)
        """
    @property
    def unspliced_bc_mod(self) -> typing.Any:
        """
        Unspliced module mols per barcode shape (steps x barcodes)
        """
    @property
    def unspliced_bc_mols(self) -> typing.Any:
        """
        Unspliced molecules per barcode shape (steps x barcodes)
        """
    @property
    def unspliced_bc_mt(self) -> typing.Any:
        """
        Unspliced MT mols per barcode shape (steps x barcodes)
        """
    @property
    def unspliced_bc_reads(self) -> typing.Any:
        """
        Unspliced reads per barcode shape (steps x barcodes)
        """
    @property
    def unspliced_mhist(self) -> numpy.typing.NDArray[numpy.uint32]:
        """
        Histogram of unspliced reads per molecule (steps x max_hist)
        """
    @property
    def unspliced_molecules(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of unspliced molecules
        """
    @property
    def unspliced_reads(self) -> numpy.typing.NDArray[numpy.uint64]:
        """
        Number of unspliced reads
        """
def aggregate_visium_bins(downsampler: Downsampler, step: typing.SupportsInt | typing.SupportsIndex, row_div: typing.SupportsInt | typing.SupportsIndex, col_div: typing.SupportsInt | typing.SupportsIndex) -> tuple:
    """
    Aggregate per-barcode outputs into Visium HD bins. Returns a dict of 1d numpy arrays.
    """
def calculate_gene_overlaps(arg0: str, arg1: collections.abc.Sequence[str]) -> dict:
    ...
