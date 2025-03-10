# cap-validator

[![PyPI version](https://img.shields.io/pypi/v/your-package-name)](https://pypi.org/project/your-package-name/)  
[![License](https://img.shields.io/github/license/cellannotation/cap-validator)](https://github.com/cellannotation/cap-validator/blob/main/LICENSE)  
[![Build Status](https://github.com/cellannotation/cap-validator/actions/workflows/unit_testing.yml/badge.svg)](https://github.com/cellannotation/cap-validator/actions)


## Overview

Python tool for validating H5AD AnnData files before uploading to the Cell Annotation Platform. The same validation code is used in [Cell Annotation Platform](https://celltype.info/) following requirements from the CAP-AnnData schema published [here](https://github.com/cellannotation/cell-annotation-schema/blob/main/docs/cap_anndata_schema.md).

## Features
- ✨ Validating all upload requirements and return results at once
- 🚀 RAM efficient
- 🧬 Provides a full list of supported ENSEMBL gene ids for Homo Sapiens and Mus Musculus


## Installation
```bash
pip install cap_validator
```

## Usage
```python
import your_package_name

your_package_name.do_something()
```

## CLI Usage (if applicable)
```bash
your-package-command --option
```

## Examples
```python
from your_package_name import FeatureClass

obj = FeatureClass(param="value")
obj.run()
```

## Development

## Contributing
1. Fork the repo
2. Create a new branch (`git checkout -b feature-branch`)
3. Make your changes and commit (`git commit -m "Description of changes"`)
4. Push to the branch (`git push origin feature-branch`)
5. Open a Pull Request

## License
[BSD 3-Clause License](LICENSE)

## Acknowledgments
- List any contributors, inspirations, or resources here.
