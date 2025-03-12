from .upload_validator import UploadValidator
from .errors import CapException

from sys import stderr
import argparse

descrption = \
"""
CLI tool to validate an AnnData H5AD file before uploading to the Cell Annotation Platform.
The validator will raise CAP-specific errors if the file does not follow the CAP AnnData Schema,
as defined in:
https://github.com/cellannotation/cell-annotation-schema/blob/main/docs/cap_anndata_schema.md

Full documentation, including a list of possible validation errors, is available at:
https://github.com/cellannotation/cap-validator/wiki

Usage Example:
    `capval path/to/anndata.h5ad`
"""


def validate():
    parser = argparse.ArgumentParser(
        description=descrption,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("adata_path", type=str, help="Path to the AnnData h5ad file.")
    args = parser.parse_args()
    adata_path = args.adata_path
    
    try:
        uv = UploadValidator(adata_path=adata_path)
        uv.validate()
    except CapException as e:
        print(e, file=stderr)
    except Exception as e:
        raise e
