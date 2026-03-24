# scdepth Binary Tag File Format

This document describes the **binary tag file format** used by `scdepth` for storing per-barcode UMI and gene count data.

---

## Overview

The file is a **BGZF-compressed binary stream** consisting of a sequence of:

```
[BarcodeHeader][UMITag][UMITag]...[UMITag]
```

Each block corresponds to a single barcode and contains:

1. A fixed-size header (`BarcodeHeader`)
2. A variable-length list of `UMITag` records

---

## Compression

* Files are compressed using **BGZF (Blocked GZIP Format)**
* Compatible with standard tools such as:

  * `htslib`
  * `bgzip`
* Supports **random access via virtual offsets**

---

## High-Level Structure

| Component       | Description                     |
| --------------- | ------------------------------- |
| `BarcodeHeader` | Metadata for a barcode block    |
| `UMITag[]`      | Array of UMI/gene count records |

---

## BarcodeHeader (12 bytes)

| Field          | Type     | Size (bytes) | Description                            |
| -------------- | -------- | ------------ | -------------------------------------- |
| `barcode`      | `uint32` | 4            | Encoded barcode identifier             |
| `n_count_tags` | `uint32` | 4            | Number of `UMITag` records             |
| `block_size`   | `uint32` | 4            | Total bytes of following `UMITag` data |

### Notes

* Total size: **12 bytes**
* Standard layout and trivially copyable
* `block_size` must equal:

  ```
  n_count_tags * sizeof(UMITag)
  ```

---

## UMITag (16 bytes)

| Field       | Type     | Size (bytes) | Description     |
| ----------- | -------- | ------------ | --------------- |
| `gene`      | `uint32` | 4            | Gene identifier |
| `umi`       | `uint32` | 4            | UMI identifier  |
| `spliced`   | `uint16` | 2            | Spliced count   |
| `unspliced` | `uint16` | 2            | Unspliced count |
| `ambiguous` | `uint16` | 2            | Ambiguous count |
| `flags`     | `uint16` | 2            | Status flags    |

### Total size: **16 bytes**

---

## Flag Definitions

| Flag              | Bit | Description                       |
| ----------------- | --- | --------------------------------- |
| `F_SPLICED_OVF`   | 0   | Spliced count overflowed `uint16` |
| `F_UNSPLICED_OVF` | 1   | Unspliced count overflow          |
| `F_AMBIG_OVF`     | 2   | Ambiguous count overflow          |
| `F_INVALID`       | 3   | Entry marked invalid              |

---

## Data Semantics

### Counts

* Counts are stored as **uint16 with saturation**
* If a count exceeds `65535`, it is:

  * Clamped to `65535`
  * Corresponding overflow flag is set

### Total Count

Total fragments for a tag:

```
spliced + unspliced + ambiguous
```

---

## File Layout Example

```
[BarcodeHeader]
  barcode = 12345
  n_count_tags = 3
  block_size = 48

[UMITag 1]
[UMITag 2]
[UMITag 3]
```

---

## Reading the File

### Sequential Access

1. Read `BarcodeHeader` (12 bytes)
2. Allocate array of size `n_count_tags`
3. Read `n_count_tags * sizeof(UMITag)` bytes
4. Validate:

   ```
   bytes_read == block_size
   ```

### Example (C++ pseudocode)

```
read BarcodeHeader
for i in n_count_tags:
    read UMITag
```

---

## Random Access via Index

`scdepth` provides a **barcode index file** that enables random access.

### Barcode Index

* Maps `barcode → BGZF virtual offset`
* Allows direct seeking to a barcode block

### Workflow

1. Look up barcode in index
2. Seek using BGZF virtual offset
3. Read `BarcodeHeader`
4. Read associated `UMITag` block

---

## Requirements

To correctly read the format:

* Use a **BGZF-aware reader** (e.g. `htslib`)
* Respect structure sizes:

  * `BarcodeHeader = 12 bytes`
  * `UMITag = 16 bytes`
* Ensure proper handling of:

  * End-of-file conditions
  * Partial reads
  * Block size validation

---

## Error Handling

Readers should validate:

* Header read size == 12 bytes
* Each `UMITag` read size == 16 bytes
* Total bytes read == `block_size`

Invalid conditions should be treated as:

* Corrupt file
* Truncated file
* I/O error

---

## Notes

* The format is designed for:

  * Efficient streaming
  * Multi-threaded processing
  * Compact storage
  * Fast per-barcode access

For example the full C++ header and reading/writing see [umis.hpp](src/umis.hpp) and [umis.cpp](src/umis.cpp)

---
