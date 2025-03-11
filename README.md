# cap-validator

[![PyPI version](https://img.shields.io/pypi/v/your-package-name)](https://pypi.org/project/your-package-name/)  
[![License](https://img.shields.io/github/license/cellannotation/cap-validator)](https://github.com/cellannotation/cap-validator/blob/main/LICENSE)  
[![Build Status](https://github.com/cellannotation/cap-validator/actions/workflows/unit_testing.yml/badge.svg)](https://github.com/cellannotation/cap-validator/actions)


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


## License
[BSD 3-Clause License](LICENSE)
