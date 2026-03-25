# scdepth

A CLI tool for sequencing depth and saturation analysis workflows.

## Requirements
 * Python ≥ 3.8
 * A working C/C++ build environment (required by scikit-build-core)
 * For the EmptyDrops wrapper R and [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
 * To use the preseq wrapper [preseq](https://github.com/smithlabcode/preseq) must be built and available

## Installation

```
git clone --recurse-submodules https://github.com/gwlab-ca/scdepth.git
cd scdepth
pip install .
```

To install rpy2 which is necessary to run the EmptyDrops wrapper use:

```
pip install -e .[emptydrops]
```

### Verify installation

```
scdepth -h
```


## Basic Usage

```
scdepth <command> [options]
```

### Note:

All commands share a sample prefix, for example, `your_sample/scdepth `.  This prefix should be the same for each command as certain output files are automatically loaded based on previous command outputs.

All of these commands depend on [cache](#cache) which generates a binary encoded and bgzip compressed tag file. For documentation on the tag file format see [here](tag_format.md).

## API

There is a large python API associated with the scdepth downsampling framework documented [here](downsampler.md) with example notebooks in [./notebooks](./notebooks). 

## Unsupported Library Types:

The [libraries](#libraries) subcommand lists all of the currently supported libraryes.  I have also added beta support for custom libraries using a [json file](custom-libs.md). Please note that this is still in beta and we are working to improve this. The custom library json file needs to be provided to the commands using the `--custom-libs` argument. I am happy to help support adding custom libraries please create an issue if you have any trouble.

## Commands and order

* [libraries](#libraries)
* [cache](#cache)
* [barcodes](#barcodes)
* [probes](#probes)
* [emptydrops](#emptydrops)
* [fit](#fit)
* [limit](#limit)
* [preseq](#preseq)
* [genes](#genes)
* [stability](#stability)

---

## Basic Command Order

| Command       | Description |
| ------------ | ----------- |
| [cache](#cache) | Cache the bam tags   |
| [barcodes](#barcodes) | Build an annotated barcode list for Visium and Visium HD samples from Space Ranger output |
| [probes](#probes) | Build a list of excluded probes. Required for Visium/Visium HD/scRNA Flex samples  |
| [emptydrops](#emptydrops) | Run emptydrops and basic filtering. Running this step or creating a similar file is required for scRNA-seq samples | 
| [fit](#fit) | Run the baseline fit.  This is required for interactive analyses and the other steps |


## Run analyses with default parameters:

```
#libraries can be viewed with scdepth libraries
scdepth cache -t <samtools reading threads> -l 10X_3p_v31 -g <gtf path> <bam path> sample/scdepth

#for probe based libraries
#scdepth probes -p <path to probe file.csv> sample/scdepth

#perform the baseline fit. This is required all of the other analyses
scdepth fit -s <sample name> -t <processing threads> sample/scdepth

#Files needed for barcode filtering. Only needed for genes/stability commands,
#  emptydrops is used with limits for pilot analysis if cell numbers are not specified.

#for Visium
#scdepth barcodes -b <path/tissue_positions.csv> sample/scdepth
#For Visium HD
#scdepth barcodes -b <path/tissue_positions.parquet> sample/scdepth
#for scRNA
#scdepth emptydrops -t <processing threads> sample/scdepth


#perform limited saturation analyses with pilot depths
#Visium
#scdepth limit -s <sample name> -t <processing threads> --tissue-frac <slide tissue fraction> --tissue-scale <6.5 or 11> sample/scdepth
#Visium HD
#scdepth limit -s <sample name> -t <processing threads> --tissue-frac <slide tissue fraction> sample/scdepth
#scRNA-seq with scdepth empydrops
#scdepth limit -s <sample name> -t <processing threads> --use-scrna sample/scdepth
#scRNA-seq with manual cell count
#scdepth limit -s <sample name> -t <processing threads> --n-cells <cells> sample/scdepth
#The limit command can also be run without pilot depths
#scdepth limit --no-cells -s <sample name> -t <processing threads> sample/scdepth

#Run preseq scdepth limit must be executed first
scdepth preseq -s <sample name> --preseq <path_to_preseq executable> sample/scdepth

#Run gene stability analyses.  Requires emptydrops/barcodes output
scdepth genes -s <sample name> -t <processing threads> sample/scdepth
#Run basic stability analysis for HVGS/KNN/Leiden clustering
scdepth stability -s <sample name> -t <processing threads> sample/scdepth
```

---

## libraries

List pre-defined libraries.

### Usage

```
scdepth libraries
```

### Options

| Option       | Description | Default / Required |
| ------------ | ----------- | ------------------ |
| `-h, --help` | Show help   | Optional           |

---

## cache

Cache molecule tags for downstream analyses.

### Usage

```
scdepth cache [options] <bam> <out>
```

### Required Options

| Option           | Description                    | Default / Required               |
| ---------------- | ------------------------------ | -------------------------------- |
| `-g, --gtf`      | Annotation GTF file            | **Required**                     |
| `-l, --lib-type` | Library type                   | **Required**                     |
| `--custom-libs`  | JSON file for custom libraries | Required if using custom library |

### Filtering Options

| Option                    | Description                      | Default / Required |
| ------------------------- | -------------------------------- | ------------------ |
| `--min-gene`              | Minimum gene overlap fraction    | 0.95               |
| `--min-gene-bases`        | Minimum gene bases               | 40                 |
| `--min-exonic`            | Minimum exonic fraction          | 0.95               |
| `--min-intronic`          | Minimum intronic overlap         | 15                 |
| `--min-qual`              | Minimum read quality             | 255                |
| `--discard-unknown-juncs` | Discard unknown splice junctions<br>Reads that contain unknown junctions will be discarded<br>(this affects cases where a read contains a known and unknown junction) | False              |

### BAM Processing

| Option           | Description        | Default / Required |
| ---------------- | ------------------ | ------------------ |
| `-t, --threads`  | Threads            | 1                  |
| `--max-tags`     | Max tags in memory | 200000000          |
| `--max-tag-frac` | Max tag fraction   | 0.95               |
| `--chunk`        | Chunk size         | 1000000            |

### Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_barcode_index.txt.gz` | The list of all the barcodes and their tag file offsets |
| `<prefix>_summary.txt` | Basic summary of the bam file including read counts, unmapped reads, library information |
| `<prefix>_genes.txt.gz` | Gene ID mapping file extracted from the gtf. This is primarily used for matrix generation and is generated for portability |
| `<prefix>_tags.gz` | Binary encoded and bgzip compressed tag file |

---

## barcodes

Build an annotated barcode list for Visium and Visium HD samples from Space Ranger output

### Usage

```
scdepth barcodes -b <barcodes> <prefix>
```

### Options

| Option           | Description  | Default / Required |
| ---------------- | ------------ | ------------------ |
| `-b, --barcodes` | Barcode file | **Required**       |

### Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_positions.parquet` | Tissue position data for 10X Visium HD libraries |
| `<prefix>_positions.csv` | Tissue position data for 10X Visium libraries |

---

## probes

Generate gene exclusion file from probe sets.

### Usage

```
scdepth probes -p <probe> <prefix>
```

### Options

| Option        | Description    | Default / Required |
| ------------- | -------------- | ------------------ |
| `-p, --probe` | Probe set file from 10X Genomics | **Required**       |

### Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_exclude.txt` | List of genes to exclude based on 10X Genomics recommendations |

---

## emptydrops

Run EmptyDrops cell-calling.

### Usage

```
scdepth emptydrops [options] <prefix>
```

### Common Options

| Option          | Description      | Default / Required |
| --------------- | ---------------- | ------------------ |
| `-s, --sample`  | Sample name      | sample             |
| `-t, --threads` | Threads          | 1                  |
| `--max-hist`    | Histogram max    | 50                 |
| `--custom-libs` | Custom libraries | Optional           |
| `-S, --seed`    | Random seed      | 42                 |

### Parameters

| Option       | Description       | Default / Required |
| ------------ | ----------------- | ------------------ |
| `-f, --frac` | Fraction of reads | 1.0                |
| `-F, --FDR`  | FDR cutoff        | 0.01               |

### Filtering Options

| Option            | Description   | Default / Required  |
| ----------------- | ------------- | ------------------- |
| `--min-molecules` | Min molecules | 500 for scRNA-seq |
| `--molecule-mads` | MAD threshold | 3.0 for scRNA-seq |

### Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_emptydrops.txt.gz` | List of barcodes with is_cell/passed annotaitons see below |

This file is a tab deliminated file and can be created by a user. It only has three fields:

| Column            | Description   | Values  |
| ----------------- | ------------- | ------------------- |
| `barcode` | barcode | Must match the barcode format in `<prefix>_barcode_index.txt.gz` |
| `is_cell` | 0 if not a cell 1 if a cell (This column is not actually used and is only there for users) |
| `passed`       | 0 if failed filtering/emptydrops and 1 if passed        | This is the column used for cell number calculations and filtering |

Example:

```
barcode is_cell passed
CTAACCCCAACCCTAA-1      0       0
CTAACCCGTAACCCTA-1      0       0
GTTCTATCAGTAACAA-1      1       1
CCTAACCGTAACCCTA-1      0       0
TCCAGAACAAATCGTC-1      1       1
TGTTCTACAGCGTACC-1      1       1
CTACCTGGTTCAAACC-1      1       1
GTCCACTAGGTAAACT-1      1       1
TCGCAGGAGCTCCGAC-1      1       1
```

---

## fit

Fit saturation curves.

### Usage

```
scdepth fit [options] <prefix>
```

### Common Options

| Option          | Description      | Default / Required |
| --------------- | ---------------- | ------------------ |
| `-s, --sample`  | Sample name      | sample             |
| `-t, --threads` | Threads          | 1                  |
| `--max-hist`    | Histogram max    | 50                 |
| `--custom-libs` | Custom libraries | Optional           |
| `-S, --seed`    | Random seed      | 42                 |

### Parameters

| Option               | Description         | Default / Required |
| -------------------- | ------------------- | ------------------ |
| `-b, --baseline-sat` | Baseline saturation | 50.0               |
| `-m, --max-sat`      | Max saturation      | 55.0               |
| `-p, --points`       | Curve points        | 8                  |
| `--sat-start`        | Min saturation      | 10                 |
| `--sat-end`          | Max saturation      | 90                 |
| `--sat-step`         | Step                | 1                  |

### Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_fit_baseline.txt` | Tab-separated file with baseline fit and count data and full dataset count data |
| `<prefix>_fit_curve.txt` | Downsampled saturation curve with `--points`` positions starting at `--sat-start` and up to and including the maximum depth for the sample |
| `<prefix>_fit_fits.svg` | Plot of the ZT-NB fit and the fit curve |
| `<prefix>_fit_rpm_histograms.txt` | Reads-per-molecule (RPM) histograms for the available read classes at the baseline depth and the full dataset depth |
| `<prefix>_fit_recovery_costs.txt` | Predicted read costs for differing molecular recovery values starting at `--sat-start` and ending at `--sat-end` using `--sat-step` as a step |
| `<prefix>_fit_saturation_costs.txt` | PRedicted read costs for differing molecular recovery values starting at `--sat-start` and ending at `--sat-end` using `--sat-step` as a step |

---

## limit

Run limited saturation analyses along a grid of saturation points using a multiple random seeds

### Usage

```
scdepth limit [options] <prefix>
```
### Common Options

| Option          | Description      | Default / Required |
| --------------- | ---------------- | ------------------ |
| `-s, --sample`  | Sample name      | sample             |
| `-t, --threads` | Threads          | 1                  |
| `--max-hist`    | Histogram max    | 50                 |
| `--custom-libs` | Custom libraries | Optional           |
| `-S, --seed`    | Random seed      | 42                 |

### Parameters

| Option          | Description         | Default / Required |
| --------------- | ------------------- | ------------------ |
| `-r, --repeats` | Number of curves    | 10                 |
| `--sat-start`   | Min saturation      | 10                 |
| `--plot-sats`   | Saturations to plot | 10,30,50,60,70,80  |
| `--sat-step`    | Step size           | 5                  |

### Depth Analysis

| Option           | Description     | Default / Required |
| ---------------- | --------------- | ------------------ |
| `--n-cells`      | Number of cells | Number of cells to use with scRNA-seq samples for pilot depth analyses |
| `--use-scrna`    | Use scRNA       | if `--n-cells` is None the EmptyDrops output file will be used for cell numbers |
| `--tissue-frac`  | Tissue fraction | Required for Visium and Visium Hd libraries for pilot depth analyses |
| `--tissue-scale` | Scale factor    | Required for Visium pilot depth analyses, possible values are 11mm (2.87) and 6.5mm (1.0)|
| `--cell-scale`   | Cell scaling    | Optional scaling factor to apply to pilot depth reads for scRNA libraries |

#### Pilot depth analyses will only be performed if `--n-cells` or `--use-scrna` is specified for scRNA-seq and `--tissue-frac` is specified for Visium/Visium HD. `--tissue-scale` should be specified for Visium. 

### Basic Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_limit_curve_stats.txt` | Countable read data and raw curve information for each seed |
| `<prefix>_limit_curves.txt` | Saturation and fit information for each point / seed |
| `<prefix>_limit_rpm.txt.gz` | Reads-per-molecule histograms for each point / seed |
| `<prefix>_limit_summary.txt` | Summary statistics for each point / seed |
| `<prefix>_limit_curve_0.svg` | Example curve plot using the `--plot-sats`  values |
| `<prefix>_limit_summary.svg` | Summary of the ZT-NB and saturation curve fit errors |

### Pilot Depth Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_limit_pilot_curve_meta.txt` | Basic information about the pilot depth curves|
| `<prefix>_limit_pilot_curves.txt` | Raw curve data for the pilot depths |
| `<prefix>_limit_pilot_rpm.txt.gz` | Reads-per-molecule histograms for each pilot point/seed |
| `<prefix>_limit_pilot_stats.txt` | Curve data for each point/seed |
| `<prefix>_limit_pilot_summary.txt` | Summary and fit error information for the pilot curves|
| `<prefix>_limit_pilot_summary.svg` | Plots of the pilot depth fits|

---

## preseq

Run preseq on limited saturation curves. This required `scdepth limit` to have already been executed.

### Usage

```
scdepth preseq --preseq <path> <prefix>
```

### Options

| Option         | Description       | Default / Required |
| -------------- | ----------------- | ------------------ |
| `--preseq`     | Preseq executable | **Required**       |
| `-s, --sample` | Sample name       | sample             |

### Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_preseq_curves.txt` | Raw curve data with added preseq predictions |
| `<prefix>_preseq_summary.txt` | Summary data with added preseq predictions |
| `<prefix>_preseq_pilot_curves.txt` | Pilot curves with added preseq predictions. Only generated if pilot depth analyses were enabled with `scdepth limit` |
| `<prefix>_preseq_pilot_summary.txt` | Pilot summaries with added preseq predictions. Only generated if pilot depth analyses were enabled with `scdepth limit` |

---

## genes

Gene overlap and recovery analysis. For scRNA-seq this requires `scdepth emptydrops` output or equivalent files. For Visium/Visium HD the command requires the `scdepth barcode` output.

### Usage

```
scdepth genes [options] <prefix>
```

### Common Options

| Option          | Description      | Default / Required |
| --------------- | ---------------- | ------------------ |
| `-s, --sample`  | Sample name      | sample             |
| `-t, --threads` | Threads          | 1                  |
| `--max-hist`    | Histogram max    | 50                 |
| `--custom-libs` | Custom libraries | Optional           |
| `-S, --seed`    | Random seed      | 42                 |

### Parameters

| Option          | Description          | Default / Required |
| --------------- | -------------------- | ------------------ |
| `-r, --repeats` | Number of curves     | 10                 |
| `--recoveries`  | Recovery targets     | 30,40,50,60,70     |
| `--write`       | Write overlap data   | False              |
| `--min_obs`     | Minimum observations | 10                 |

### Filtering

| Option            | Description   | Default / Required  |
| ----------------- | ------------- | ------------------- |
| `--min-molecules` | Min molecules | 500 for Visium, 10 (1x1), 20 (2x2 and 3x3) and 50 (4x4) and 100 for everything else<br>See the `--vis-bins argument`). <br>*Not used for scRNA-seq since the emptydrops file is already filtered* |
| `--molecule-mads` | MAD threshold | 3 for Visium 4 for Visium HD. <br>*Not used for scRNA-seq since the emptydrops file is already filtered* |

### Visium HD Options  

| Option            | Description   | Default / Required  |
| ----------------- | ------------- | ------------------- |
| `--vis-rows` | Rows in the slide | 3350 |
| `--vis-cols` | Columns in the slide | 3350 |
| `--vis-bins` | Bin size (BxB bins) | 8 |

### Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_genes_curve_stats.txt` | Raw curve stats for the molecular recovery curves |
| `<prefix>_genes_summary.txt` | Summary data for the gene stability analyses |
| `<prefix>_genes_barcode_summary.svg` | Violin plots of the per barcode gene detection stability for the global and molecule classes |


---

## stability

Barcode clustering stability analysis. For scRNA-seq this requires `scdepth emptydrops` output or equivalent files. For Visium/Visium HD the command requires the `scdepth barcode` output.

### Usage

```
scdepth stability [options] <prefix>
```

### Common Options

| Option          | Description      | Default / Required |
| --------------- | ---------------- | ------------------ |
| `-s, --sample`  | Sample name      | sample             |
| `-t, --threads` | Threads          | 1                  |
| `--max-hist`    | Histogram max    | 50                 |
| `--custom-libs` | Custom libraries | Optional           |
| `-S, --seed`    | Random seed      | 42                 |

### Parameters

| Option          | Description             | Default / Required |
| --------------- | ----------------------- | ------------------ |
| `-r, --repeats` | Number of curves        | 10                 |
| `--recoveries`  | Recovery targets        | 30,40,50,60,70     |
| `--resolution`  | Clustering resolution   | 0.6,0.8,1.0        |
| `--min-cells`   | Min cells for filtering | 10                 |
| `--hvgs`        | Highly variable genes   | 2000               |
| `--scale`       | Normalization scale     | 10000              |
| `--pcs`         | Principal components    | 30                 |
| `--neighbours`  | Neighbors               | 30                 |

### Filtering

| Option            | Description   | Default / Required  |
| ----------------- | ------------- | ------------------- |
| `--min-molecules` | Min molecules | 500 for Visium, 10 (1x1), 20 (2x2 and 3x3) and 50 (4x4) and 100 for everything else<br>See the `--vis-bins argument`). <br>*Not used for scRNA-seq since the emptydrops file is already filtered* |
| `--molecule-mads` | MAD threshold | 3 for Visium 4 for Visium HD. <br>*Not used for scRNA-seq since the emptydrops file is already filtered* |

### Visium HD Options  

| Option            | Description   | Default / Required  |
| ----------------- | ------------- | ------------------- |
| `--vis-rows` | Rows in the slide | 3350 |
| `--vis-cols` | Columns in the slide | 3350 |
| `--vis-bins` | Bin size (BxB bins) | 8 |

### Outputs

| File        | Description |
| ----------- | ----------- |
| `<prefix>_stability_curve_stats.txt` | raw curve stats for the recovery curves |
| `<prefix>_stability_summary.txt` | Summary statistics of the HVGS/KNN/clustering stability analyses |
| `<prefix>_stability_summary.svg` | Clustering stability summary plot |

---


