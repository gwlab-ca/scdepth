# Basic Downsampler Python API

Python bindings for the `Downsampler` class.

Note, this documentation is still a work in progress, see the example jupyter notebooks
---

## Class

```python
from scdepth.bindings import Downsampler
Downsampler()
```

Create a new `Downsampler` instance.

---

## Methods

### init

```python
init(
    prefix,
    mt_prefix="",
    mt_file="",
    mod_file="",
    exclude_file="",
    barcode_prefix="",
    max_hist=50,
    build_matrices=False,
    calc_sau=False,
)
```

Initialize the downsampler from a cached prefix.

---

### init_visium

```python
init_visium(
    rows,
    cols,
    in_tissue,
    countable,
    total,
    total_rows,
    total_cols,
)
```

Initialize the downsampler for Visium HD tiled/barcode data.

---

### downsample

```python
downsample(
    fracs,
    umi_len,
    seed,
    threads=1,
    aggregate_only=False,
    umi_mode="directed",
    correct_multi_umis=True,
)
```

Run downsampling.

---

### reset_visium

```python
reset_visium()
```

Reset Visium tile aggregation information.

---

### clear_output

```python
clear_output()
```

Clear output memory.

---

### Full list of class methods and properties

| Name | Type | Description |
|---|---|---|
| `clear_output` | method | `clear_output(self: scdepth.bindings._bindings.Downsampler) -> None` Clear output memory |
| `downsample` | method | `downsample(self: scdepth.bindings._bindings.Downsampler, fracs: collections.abc.Sequence[typing.SupportsFloat], umi_len: typing.SupportsInt, seed: typing.SupportsInt, threads: typing.SupportsInt = 1, aggregate_only: bool = False, umi_mode: str = 'directed', correct_multi_umis: bool = True) -> bool` Run downsampling for the given barcodes (optional) and fractions |
| `init` | method | `init(self: scdepth.bindings._bindings.Downsampler, prefix: str, mt_prefix: str = '', mt_file: str = '', mod_file: str = '', exclude_file: str = '', max_hist: typing.SupportsInt = 50, build_matrices: bool = False, calc_sau: bool = False) -> bool` Initialize the downsampler from a cached prefix |
| `init_visium` | method | `init_visium(self: scdepth.bindings._bindings.Downsampler, rows: collections.abc.Sequence[typing.SupportsInt], cols: collections.abc.Sequence[typing.SupportsInt], in_tissue: collections.abc.Sequence[typing.SupportsInt], countable: collections.abc.Sequence[typing.SupportsInt], total: collections.abc.Sequence[typing.SupportsInt], total_rows: typing.SupportsInt, total_cols: typing.SupportsInt) -> bool` Initialize the downsampler to create binned tiles/barcodes for visium HD data |
| `reset_visium` | method | `reset_visium(self: scdepth.bindings._bindings.Downsampler) -> bool` Reset visium tile aggregation information |
| `spliced_csr` | method | `spliced_csr(self: scdepth.bindings._bindings.Downsampler, step: typing.SupportsInt) -> object` |
| `total_csr` | method | `total_csr(self: scdepth.bindings._bindings.Downsampler, step: typing.SupportsInt) -> object` |
| `total_csr_bin` | method | `total_csr_bin(self: scdepth.bindings._bindings.Downsampler, step: typing.SupportsInt, bin_div: typing.SupportsInt) -> object` |
| `unspliced_csr` | method | `unspliced_csr(self: scdepth.bindings._bindings.Downsampler, step: typing.SupportsInt) -> object` |
| `ambiguous_csr` | method | `ambiguous_csr(self: scdepth.bindings._bindings.Downsampler, step: typing.SupportsInt) -> object` |
| `aggregate_only` | property | Whether sparse matrices were built |
| `barcode_prefix` | property | User specified primer mode |
| `priemr_mode` | property | User specified barcode prefix |
| `ambiguous_bc_genes` | property | Ambiguous genes per barcode shape (steps x barcodes) |
| `ambiguous_bc_mod` | property | Ambiguous module mols per barcode shape (steps x barcodes) |
| `ambiguous_bc_mols` | property | Ambiguous molecules per barcode shape (steps x barcodes) |
| `ambiguous_bc_mt` | property | Ambiguous MT mols per barcode shape (steps x barcodes) |
| `ambiguous_bc_reads` | property | Ambiguous reads per barcode shape (steps x barcodes) |
| `ambiguous_mhist` | property | Histogram of ambiguous reads per MT molecule (steps x max_hist) |
| `ambiguous_molecules` | property | Number of ambiguous molecules |
| `ambiguous_reads` | property | Number of ambiguous reads |
| `barcodes` | property | Barcode columns as a dict-of-arrays (DataFrame-friendly) |
| `bc_discarded_moleculess` | property | Molecules lost from multi-gene UMIs (steps x barcodes) |
| `bc_discarded_reads` | property | Reads lost from multi-gene UMIs and excluded gene lists (steps x barcodes) |
| `bc_excluded_reads` | property | Reads lost from excluded read lists (steps x barcodes) |
| `build_mats` | property | Whether sparse matrices were built |
| `fracs` | property | Downsampling fractions actually used |
| `has_MT` | property | A MT gene list was included |
| `has_excluded` | property | An excluded gene list was included |
| `has_module` | property | A module gene list was included |
| `has_sau` | property | Has spliced/ambig/unspliced quantifications |
| `max_hist` | property | Maximum histogram bin (histograms have size max_hist+1) |
| `molecules_ambig` | property | Number of molecule subgraphs mapping to at least 2 genes. |
| `molecules_discarded` | property | Number of molecules lost to ambiguous UMI/gene mappings |
| `reads_discarded` | property | Number of reads lost to ambiguous UMI/gene mappings |
| `reads_excluded` | property | Number of reads lost to excluded gene filter |
| `seed` | property | Downsampling seed |
| `spliced_bc_genes` | property | Spliced genes per barcode shape (steps x barcodes) |
| `spliced_bc_mod` | property | Spliced module mols per barcode shape (steps x barcodes) |
| `spliced_bc_mols` | property | Spliced molecules per barcode shape (steps x barcodes) |
| `spliced_bc_mt` | property | Spliced MT mols per barcode shape (steps x barcodes) |
| `spliced_bc_reads` | property | Spliced reads per barcode shape (steps x barcodes) |
| `spliced_mhist` | property | Histogram of spliced reads per molecule (steps x max_hist) |
| `spliced_molecules` | property | Number of spliced molecules |
| `spliced_reads` | property | Number of spliced reads |
| `total_barcodes` | property | Number of raw barcodes in the results |
| `total_bc_genes` | property | Total genes per barcode shape (steps x barcodes) |
| `total_bc_mod` | property | Total module mols per barcode shape (steps x barcodes) |
| `total_bc_mols` | property | Total molecules per barcode shape (steps x barcodes) |
| `total_bc_mt` | property | Total MT mols per barcode shape (steps x barcodes) |
| `total_bc_reads` | property | Total reads per barcode shape (steps x barcodes) |
| `total_fracs` | property | Number of downsampling fractions |
| `total_genes` | property | Number of genes in the results |
| `total_mhist` | property | Histogram of total reads per molecule (steps x max_hist) |
| `total_molecules` | property | Number of total molecules |
| `total_reads` | property | Number of total reads |
| `unspliced_bc_genes` | property | Unspliced genes per barcode shape (steps x barcodes) |
| `unspliced_bc_mod` | property | Unspliced module mols per barcode shape (steps x barcodes) |
| `unspliced_bc_mols` | property | Unspliced molecules per barcode shape (steps x barcodes) |
| `unspliced_bc_mt` | property | Unspliced MT mols per barcode shape (steps x barcodes) |
| `unspliced_bc_reads` | property | Unspliced reads per barcode shape (steps x barcodes) |
| `unspliced_mhist` | property | Histogram of unspliced reads per molecule (steps x max_hist) |
| `unspliced_molecules` | property | Number of unspliced molecules |
| `unspliced_reads` | property | Number of unspliced reads |
