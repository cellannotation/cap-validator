from typing import List, Optional


class CapException(BaseException):
    name = "Unknown"
    message = "Generic CAP exception"

    def __str__(self) -> str:
        return f"{self.name}: {self.message}"

    def __eq__(self, other):
        if isinstance(other, CapException):
            return (self.name, self.message) == (other.name, other.message)
        raise TypeError(f"The __eq__ operation is not defined for CapException and {type(other)}!")

    def __hash__(self):
        return hash((self.name, self.message))


class BadAnnDataFile(CapException):
    name = 'Unknown'
    message = 'The file format is not supported!'


class CapMultiException(CapException):
    """
    Class to raise multiple errors at once
    Usage example
    e = CapMultiException()
    if ...:
        e.append(CapException())  # append error instead of raise it
    if ...:
        e.append(BadAnndataFile())
    
    if e.have_errors():
        raise e
    """
    name = "CapMultiException"
    message = ""
    ex_list: list = None
    raise_on_append: bool = False  # for debug and tests

    def __init__(self, message: str=""):
        """This init is important to be in class, 
        as of it fix the behaviour where python re-use CapMultiException 
        with exist 'ex_list' on each validation call."""
        self.message = message
        self.ex_list = []

    def append(self, other: CapException) -> None:
        if isinstance(other, CapException):
            if self.raise_on_append:
                raise other
            else:
                self.ex_list.append(other)

    def __str__(self) -> str:
        own_str = super().__str__()
        res_list = [own_str] + self.ex_list
        res_message = "\n".join(map(str, res_list))
        res_message += "\nFor details visit: \n\thttps://github.com/cellannotation/cap-validator/wiki/Validation-Errors"
        return res_message

    def have_errors(self) -> bool:
        return len(self.ex_list) > 0

    def __getitem__(self, item: int) -> CapException:
        return self.ex_list[item]

    def to_list(self) -> List[CapException]:
        return self.ex_list


class AnnDataMissingCountMatrix(CapException):
    name = "AnnDataMissingCountMatrix"
    message = "Matrix is missing in both `.X` and `.raw.X`."


class AnnDataInvalidCountMatrix(CapException):
    name = "AnnDataInvalidCountMatrix"
    message = "Raw counts matrix values must be non-negative integers."


class AnnDataMissingEmbeddings(CapException):
    name = "AnnDataMissingEmbeddings"
    message = (
        "Embeddings are missing or incorrectly formatted. "
        "They must be stored in `.obsm` as [n_cells × 2] datasets "
        "with names starting with 'X_' (e.g., X_umap, X_tsne)."
    )


class AnnDataMissingObs(CapException):
    name = "AnnDataMissingObs"
    message = "The `.obs` is missing."


class AnnDataMissingObsColumns(CapException):
    name = "AnnDataMissingObsColumns"

    def __init__(self, missing_columns: list[str] = None):
        msg = (
            "Required `.obs` metadata is missing. The file must contain "
            "'assay', 'disease', 'organism', 'tissue' "
            "or corresponding '<x>_ontology_term_id' fields."
        )
        if missing_columns:
            cols_str = ", ".join(missing_columns)
            self.message = f"{msg}\nMissing columns: {cols_str}"
        else:
            self.message = msg


class AnnDataEmptyOrNoneInGeneralMetadata(CapException):
    name = "AnnDataEmptyOrNoneInGeneralMetadata"

    def __init__(
        self,
        none_columns: Optional[List[str]] = None,
        empty_columns: Optional[List[str]] = None,
    ):
        msg = (
            "Required `.obs` metadata metadata contains missing or invalid values.\n"
            "All required fields must be filled with valid values."
        )

        if none_columns:
            msg += "\nColumns with None / NaN values: " + ", ".join(sorted(none_columns))

        if empty_columns:
            msg += "\nColumns with empty values: " + ", ".join(sorted(empty_columns))

        self.message = msg


class AnnDataNonStandardVarError(CapException):
    name = "AnnDataNonStandardVarError"


class AnnDataMissingVarIndex(AnnDataNonStandardVarError):
    name = "AnnDataMissingVarIndex"
    message = "The `.var.index` is missing or empty."


class AnnDataNumericVarIndex(AnnDataNonStandardVarError):
    name = "AnnDataNumericVarIndex"
    message = "The `.var.index` contains numeric values instead of gene identifiers."


class AnnDataVarNotSubsetOfRawVar(AnnDataNonStandardVarError):
    name = "AnnDataVarNotSubsetOfRawVar"
    message = "`var.index` must be a subset of `raw.var.index`."


class AnnDataGeneIndexIsNotUnique(AnnDataNonStandardVarError):
    name = "AnnDataGeneIndexIsNotUnique"
    message = "`var.index` must not contain duplicates."


class AnnDataUnsupportedGenes(AnnDataNonStandardVarError):
    name = "AnnDataUnsupportedGenes"

    def __init__(self, missing_genes_count: int = None):
        msg = (
            "File does not contain valid ENSEMBL terms in var.\n"
            "We currently support only Homo sapiens and Mus musculus.\n"
            "In the case of multiple species in the dataset, orthologous Homo sapiens genes are required.\n"
            "If there are other species you wish to upload to CAP, please contact "
            "support@celltype.info and we will work to accommodate your request."
        )
        if missing_genes_count is not None:
            msg += f"\nNumber of unsupported genes found: {missing_genes_count}"
        self.message = msg


class CSCMatrixInX(CapException):
    name = "CSCMatrixInX"

    def __init__(self, locations: list[str]):
        """
        locations: list of matrix locations, e.g. ['X'], ['raw.X'], or ['X', 'raw.X']
        """
        super().__init__()
        self.locations = locations
        loc_str = " and ".join(locations)
        self.message = (
            f"A CSC matrix is found in {loc_str}. "
            "The gene expression matrix must be stored in CSR or dense format!"
        )
