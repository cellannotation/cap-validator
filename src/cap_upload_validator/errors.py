from typing import List


class CapException(BaseException):
    name = "Unknown"
    message = "Useless CAP exception"

    def __str__(self) -> str:
        return f"{self.name}: {self.message}"

    def __eq__(self, other):
        if isinstance(other, CapException):
            return (self.name, self.message) == (other.name, other.message)
        raise TypeError(f"The __eq__ operation doesn't defined for CapException and {type(other)}!")

    def __hash__(self):
        return hash((self.name, self.message))


class BadAnnDataFile(CapException):
    name = 'Unknown'
    message = 'File format is not supported!'


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


class AnnDataFileMissingCountMatrix(CapException):
    name = "AnnDataFileMissingCountMatrix"
    message = "DataFile Incorrect format: raw data matrix is missing in .raw.X or .X."


class AnnDataMissingEmbeddings(CapException):
    name = "AnnDataMissingEmbeddings"
    message = \
        """
        The embedding is missing or is incorrectly named: embeddings must be a [n_cells x 2] 
        numpy array saved with the prefix X_, for example: X_tsne, X_pca or X_umap.
        """


class AnnDataMisingObsColumns(CapException):
    name = "AnnDataMisingObsColumns"
    message = \
        """
            Required obs column(s) missing: file must contain 
            'assay', 'disease', 'organism' and 'tissue' fields with valid values.
        """


class AnnDataNonStandardVarError(CapException):
    name = "AnnDataNonStandardVarError"
    message = \
        """
            File does not contain valid ENSEMBL terms in var.
            We currently support Homo sapiens and Mus musculus.
            In the case of multiple species in the dataset, orthologous Homo sapiens genes are required.
            If there are other species you wish to upload to CAP, please contact
            support@celltype.info and we will work to accommodate your request.
        """ 
