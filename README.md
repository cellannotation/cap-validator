# cap-validator

[![PyPI version](https://img.shields.io/pypi/v/cap-upload-validator)](https://pypi.org/project/cap-upload-validator/)  [![License](https://img.shields.io/github/license/cellannotation/cap-validator)](https://github.com/cellannotation/cap-validator/blob/main/LICENSE)  [![Build Status](https://github.com/cellannotation/cap-validator/actions/workflows/unit_testing.yml/badge.svg)](https://github.com/cellannotation/cap-validator/actions)


## Overview

Python tool for validating H5AD AnnData files before uploading to the Cell Annotation Platform. The same validation code is used in [Cell Annotation Platform](https://celltype.info/) following requirements from the CAP-AnnData schema published [here](https://github.com/cellannotation/cell-annotation-schema/blob/main/docs/cap_anndata_schema.md).

## Features
- ✨ Validates all upload requirements and returns results at once
- 🚀 RAM efficient
- 🧬 Provides a full list of supported ENSEMBL gene IDs for *Homo sapiens* and *Mus musculus*


## Installation
```bash
pip install -U cap-upload-validator
```

## Usage

### Basic usage

```python
from cap_upload_validator import UploadValidator

h5ad_path = "path_to.h5ad"

uv = UploadValidator(h5ad_path)
uv.validate()
```

Full documentation could be found in [GitHub Wiki](https://github.com/cellannotation/cap-validator/wiki)

### CLI interface
```console
$ capval tmp/tmp.h5ad
CapMultiException: 
AnnDataMissingEmbeddings: 
        The embedding is missing or is incorrectly named: embeddings must be a [n_cells x 2] 
        numpy array saved with the prefix X_, for example: X_tsne, X_pca or X_umap.
        
AnnDataMisingObsColumns: 
            Required obs column(s) missing: file must contain 
            'assay', 'disease', 'organism' and 'tissue' fields with valid values.
        
For details visit: 
        https://github.com/cellannotation/cap-validator/wiki/Validation-Errors

$ capval --help
usage: capval [-h] adata_path

CLI tool to validate an AnnData H5AD file before uploading to the Cell Annotation Platform.
The validator will raise CAP-specific errors if the file does not follow the CAP AnnData Schema,
as defined in:
https://github.com/cellannotation/cell-annotation-schema/blob/main/docs/cap_anndata_schema.md

Full documentation, including a list of possible validation errors, is available at:
https://github.com/cellannotation/cap-validator/wiki

Usage Example:
    `capval path/to/anndata.h5ad`

positional arguments:
  adata_path  Path to the AnnData h5ad file.

optional arguments:
  -h, --help  show this help message and exit
```

## License
[BSD 3-Clause License](LICENSE)
