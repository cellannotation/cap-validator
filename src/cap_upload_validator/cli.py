from .upload_validator import UploadValidator
from .errors import CapException
from sys import stderr
from os import PathLike
from typing import Union
import click

@click.command()
@click.argument("adata_path", type=click.Path(exists=True, readable=True, path_type=str))
def validate(
    adata_path: Union[str, PathLike]
) -> None:
    try:
        uv = UploadValidator(adata_path=adata_path)
        uv.validate()
    except CapException as e:
        print(e, file=stderr)
    except Exception as e:
        raise e
