#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

from types import SimpleNamespace
import json
import re
from pathlib import Path
from typing import Any, Mapping

columns = ['five_prime_like', 'CB_tag', 'CB_length', 'CB_re', 
            'UR_tag', 'UR_length', 'sample_group', 'probe_based',
           'has_bins']
visium_re = r"^s_\d{3}um_\d{5}_\d{5}-\d+$"
seq_re = r"^[ACGT]+-\d+$"
DEFAULT_LIBS = {
    #scRNA 3-prime
    '10X_3p_v2':   (False, 'CB', 0, seq_re, 'UR', 10, 'scrna', False, False),
    '10X_3p_v3':   (False, 'CB', 0, seq_re, 'UR', 12, 'scrna', False, False),
    '10X_3p_v31':  (False, 'CB', 0, seq_re, 'UR', 12, 'scrna', False, False),
    '10X_3p_v4':   (False, 'CB', 0, seq_re, 'UR', 12, 'scrna', False, False),

    #scRNA 5-prime
    '10X_5p_v2':   (True,  'CB', 0, seq_re, 'UR', 10, 'scrna', False, False),
    '10X_5p_v3':   (True,  'CB', 0, seq_re, 'UR', 12, 'scrna', False, False),

    #probe based
    '10X_flex_v1': (True,  'CB', 0, seq_re, 'UR', 12, 'scrna', True, False),
    '10X_visium': (True,  'CB', 0, seq_re, 'UR', 12, 'visium', True, True),
    '10X_visium_hd': (True,  'CB', 0, visium_re, 'UR', 9, 'visium_hd', True, True),

    #non-probe based spatial
    '10X_visium_3p': (False,  'CB', 0, seq_re, 'UR', 12, 'visium', False, True),
    '10X_visium_hd_3p': (False,  'CB', 0, visium_re, 'UR', 9, 'visium_hd', False, True),
}

libs: dict[str, tuple] = dict(DEFAULT_LIBS)

class LibrarySpecError(ValueError):
    pass


def _normalize_and_validate_spec(name: str, spec: Mapping[str, Any]) -> tuple:
    missing = [c for c in columns if c not in spec]
    if missing:
        raise LibrarySpecError(f'{name}: missing fields: {missing}')

    # Basic type/constraints checks (tighten as you like)
    if not isinstance(spec['five_prime_like'], bool):
        raise LibrarySpecError(f'{name}: five_prime_like must be bool')

    for k in ('CB_tag', 'UR_tag', 'sample_group'):
        if not isinstance(spec[k], str) or not spec[k]:
            raise LibrarySpecError(f'{name}: {k} must be a non-empty string')

    for k in ('CB_length', 'UR_length'):
        if not isinstance(spec[k], int) or spec[k] < 0:
            raise LibrarySpecError(f'{name}: {k} must be int >= 0')

    if not isinstance(spec['probe_based'], bool) or not isinstance(spec['has_bins'], bool):
        raise LibrarySpecError(f'{name}: probe_based/has_bins must be bool')

    cb_re = spec['CB_re']
    if not isinstance(cb_re, str) or not cb_re:
        raise LibrarySpecError(f'{name}: CB_re must be a non-empty regex string')
    try:
        re.compile(cb_re)
    except re.error as e:
        raise LibrarySpecError(f'{name}: CB_re is not a valid regex: {e}') from e

    return tuple(spec[c] for c in columns)


def add_custom_libraries(path: str | Path, *, override: bool = False) -> None:
    """
    Load user library definitions from a JSON file and merge into the runtime `libs`.

    File format (recommended):
      {
        'my_lib': {
          'five_prime_like': true,
          'CB_tag': 'CB',
          'CB_length': 16,
          'CB_re': '^[ACGT]{16}$',
          'UR_tag': 'UR',
          'UR_length': 12,
          'sample_group': 'scrna',
          'probe_based': false,
          'has_bins': false
        }
      }
    """
    p = Path(path)
    data = json.loads(p.read_text())

    if not isinstance(data, dict):
        raise LibrarySpecError('Top-level JSON must be an object mapping lib_name -> spec')

    for name, spec in data.items():
        if not isinstance(name, str) or not name:
            raise LibrarySpecError(f'Invalid library name: {name!r}')
        if not isinstance(spec, dict):
            raise LibrarySpecError(f'{name}: spec must be an object')

        if (name in libs) and not override:
            raise LibrarySpecError(
                f'{name}: conflicts with an existing library. '
                f'Use override=True / --lib-override to replace it.'
            )

        libs[name] = _normalize_and_validate_spec(name, spec)

def print_libs():
    import pandas as pd
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', len(libs))
    df = pd.DataFrame.from_dict(libs, orient='index',
                            columns=columns)
    print(df)

def library2args(lib: str) -> tuple:
    return libs[lib]

def library2ns(lib: str) -> SimpleNamespace:
    return SimpleNamespace(library=lib, **{k: v for k, v in zip(columns, libs[lib])})

def is_library(lib: str) -> bool:
    return lib in libs

def is_probe_based(lib: str) -> bool:
    if lib not in libs:
        return False
    return library2ns(lib).probe_based

def has_bins(lib: str) -> bool:
    if lib not in libs:
        return False
    return library2ns(lib).has_bins
