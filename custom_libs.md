# Custom Libraries (Beta)

> **Experimental Feature**
>
> Custom library support is currently in **beta/testing**.
> The schema and validation rules may change in future releases.
> Use with caution and validate results carefully.

---

## Overview

`scdepth` includes a set of built-in library definitions (e.g. `10X_3p_v3`, `10X_visium`, etc.), but you can extend these by providing your own **custom library specifications**.

Custom libraries are defined in a JSON file and loaded at runtime using:

```
--custom-libs <path>
```

Some commands (e.g. `cache`, `emptydrops`, `fit`, etc.) require this flag when using a non-default library.

---

## File Format

Custom libraries are defined as a JSON object mapping **library name → specification**:

```json
{
  "my_library": {
    "five_prime_like": true,
    "CB_tag": "CB",
    "CB_length": 16,
    "CB_re": "^[ACGT]{16}$",
    "UR_tag": "UR",
    "UR_length": 12,
    "sample_group": "scrna",
    "probe_based": false,
    "has_bins": false,
    "random_hex_re": "",
    "random_hex_value": "",
  }
}
```

---

## Field Definitions

| Field              | Type           | Description                                                                                                                    |
| ------------------ | -------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `five_prime_like`  | bool           | Whether the library behaves like a 5′ protocol                                                                                 |
| `CB_tag`           | string         | BAM tag used for cell barcode                                                                                                  |
| `CB_length`        | int ≥ 0        | Expected barcode length (0 = variable)                                                                                         |
| `CB_re`            | string (regex) | Regex pattern used to validate barcodes                                                                                        |
| `UR_tag`           | string         | BAM tag used for UMI                                                                                                           |
| `UR_length`        | int ≥ 0        | Expected UMI length                                                                                                            |
| `sample_group`     | string         | One of: `scrna`, `visium`, `visium_hd`, etc.                                                                                   |
| `probe_based`      | bool           | Whether the protocol is probe-based                                                                                            |
| `has_bins`         | bool           | Whether spatial binning is used                                                                                                |
| `random_hex_re`    | string (regex) | Regex pattern used to extract if a read is random hexamer based from the qname, there must be a single capture group or blank, this RE should also only accept valid values not blankly take anything as its the primary filter  |
| `random_hex_value` | string         | The expected value from the RE pattern or blank                                                                                |

---

## Validation Rules

When loading custom libraries, `scdepth` performs strict validation:

* All fields listed above are **required**
* `CB_re` must be a **valid regular expression**
* `CB_length` and `UR_length` must be **≥ 0 integers**
* `CB_tag`, `UR_tag`, and `sample_group` must be **non-empty strings**
* `probe_based` and `has_bins` must be **boolean**

Invalid specifications will raise an error and stop execution.

---

## Example

### Example: custom scRNA library

```json
{
  "custom_scrna_v1": {
    "five_prime_like": false,
    "CB_tag": "CB",
    "CB_length": 16,
    "CB_re": "^[ACGT]{16}$",
    "UR_tag": "UR",
    "UR_length": 10,
    "sample_group": "scrna",
    "probe_based": false,
    "has_bins": false
  }
}
```

---

## Using Custom Libraries

### 1. Provide the JSON file

```
scdepth cache \
  -g genes.gtf \
  -l custom_scrna_v1 \
  --custom-libs custom_libs.json \
  input.bam output_prefix
```

### 2. Use in downstream commands

Any command that depends on library metadata must also include:

```
--custom-libs custom_libs.json
```

---

## Overriding Existing Libraries

By default, custom libraries **cannot overwrite built-in libraries**.

If a name conflict occurs, an error will be raised.

Override behavior is only allowed if explicitly enabled (internal / advanced usage):

```
override=True
```

or via CLI (if exposed):

```
--lib-override
```

---

## Summary

Custom libraries provide flexibility for non-standard protocols, but:

> They are **powerful but experimental** — use carefully.
