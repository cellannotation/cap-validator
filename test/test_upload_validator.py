import pytest
import numpy as np
import pandas as pd

import anndata as ad
from packaging import version
if version.parse(ad.__version__) >= version.parse("0.11.0"):
    ad.settings.allow_write_nullable_strings = True

import scipy.sparse as sp
from pathlib import Path
import tempfile
from cap_anndata import CapAnnDataDF, read_h5ad
from contextlib import nullcontext

from cap_upload_validator.upload_validator import (
    UploadValidator,
    GENERAL_METADATA,
    ORGANISM_COLUMN,
    ORGANISM_ONT_ID_COLUMN,
)
from cap_upload_validator.gene_mapping import (
    GeneMap,
    HomoSapiens,
)
from cap_upload_validator.errors import (
    AnnDataMissingEmbeddings,
    AnnDataMissingObs,
    AnnDataMissingObsColumns,
    AnnDataNonStandardVarError,
    CapMultiException,
    AnnDataEmptyOrNoneInGeneralMetadata,
    CSCMatrixInX,
    AnnDataMultipleOntologyIDs,
    AnnDataInvalidDiseaseOntologyForHuman,
    AnnDataGeneIndexIsNotUnique,
    AnnDataUnsupportedGenes,
)

TMP_DIR = Path(tempfile.mkdtemp())

@pytest.mark.parametrize("expected_with_data", [
    (True, np.array([
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 0]
    ])),
    (False, np.array([
        [0.0, 1.5, 2.5],
        [3.5, 4.5, 5.5],
        [6.5, 7.5, 0.0]
    ])),
    (False, np.array([
        [0.0, -1.0, -2.0],
        [-3.0, 4.0, -5.0],
        [6.0, 7.0, 0.0]
    ])),
])
@pytest.mark.parametrize("sparse", [True, False])
def test_is_positive_integers(sparse, expected_with_data):
    expected, X = expected_with_data
    X = X.astype(np.float32)
    adata = ad.AnnData(X=sp.csr_matrix(X, dtype=np.float32) if sparse else X)
    v = UploadValidator(adata_path=None)
    assert v._check_is_positive_integers(adata) == expected, "Incorrect X matrix validation!"


def test_has_embeddings():
    file_path = TMP_DIR / "test_has_embeddings.h5ad"
    emb_name = "X_test"
    
    adata = ad.AnnData(X=np.eye(10))
    adata.obsm[emb_name] = np.ones(shape=(adata.shape[0], 2))
    adata.write_h5ad(file_path)

    with read_h5ad(file_path, edit=False) as cap_adata:
        v = UploadValidator(adata_path=file_path)
        v._multi_exception.raise_on_append = True
        try:
            v._check_obsm(cap_adata)
        except:
            assert False, "Must be embeddings in file!"

        del cap_adata.obsm[emb_name]

        try:
            v._check_obsm(cap_adata)
            assert False, "Must not be embeddings in file!"
        except AnnDataMissingEmbeddings:
            pass


def test_obs():
    file_path = TMP_DIR / "test_obs.h5ad"

    adata = ad.AnnData(X=np.eye(10))
    for col in GENERAL_METADATA:
        adata.obs[col] = "test_value"
    
    adata.write_h5ad(file_path)
    
    with read_h5ad(file_path, edit=False) as cap_adata:
        v = UploadValidator(adata_path=file_path)
        v._multi_exception.raise_on_append = True
        cap_adata.read_obs()
        df = cap_adata.obs.copy()

        def check_obs(ca, correct_expected: bool):
            try:
                v._check_obs(ca)
                if not correct_expected:
                    assert False, "Must not be correct obs!"
            except (AnnDataMissingObsColumns, AnnDataMissingObs):
                assert not correct_expected, "Unexpected result"

        check_obs(cap_adata, True)

        for col in GENERAL_METADATA:
            cap_adata.obs = CapAnnDataDF.from_df(df.drop(col, axis=1, inplace=False))
            check_obs(cap_adata, False)

        cap_adata.obs = CapAnnDataDF.from_df(df[[]])
        check_obs(cap_adata, False)


def test_var_index():    
    v = UploadValidator(None)
    gene_map = GeneMap.data_frame()
    
    adata = ad.AnnData(X=np.eye(10))
    n_genes = adata.shape[1]
    
    file_path = TMP_DIR / "test_adata.h5ad"
    def check_var_index():
        with read_h5ad(file_path, edit=False) as cap_adata:
            cap_adata.read_var()
            cap_adata.read_obs(columns=[ORGANISM_COLUMN])
            v._check_var_index(cap_adata)

    adata.write_h5ad(filename=file_path)
    
    # proper organism and genes
    adata.obs[ORGANISM_COLUMN] = HomoSapiens.name
    adata.var.index = gene_map.ENSEMBL_gene[:n_genes]
    adata.write_h5ad(filename=file_path)
    check_var_index()

    # proper organism and genes with version suffixes
    adata.obs[ORGANISM_COLUMN] = HomoSapiens.name
    adata.var.index = [f"{g}.{i}" for i, g in enumerate(gene_map.ENSEMBL_gene[:n_genes])]
    adata.write_h5ad(filename=file_path)
    check_var_index()

    # proper organism and bad genes
    adata.var.index = map(str, range(n_genes))
    adata.write_h5ad(filename=file_path)
    try:
        check_var_index()
    except AnnDataNonStandardVarError:
        pass
    except Exception as e:
        assert False, f"Unpredicted error: {e}"

    # unsuported organism and bad genes
    adata.obs.organism = 'unsuported'
    adata.write_h5ad(filename=file_path)
    try:
        check_var_index()
    except:
        assert False, "Wrong validation failure for unsuported organism!"

    # multiple organisms and non ENSG genes
    adata.obs['organism'] = adata.obs['organism'].cat.add_categories(['new organism'])
    adata.obs.loc['0', 'organism'] = 'new organism'

    adata.write_h5ad(filename=file_path)
    try:
        check_var_index()
    except AnnDataNonStandardVarError:
        pass
    except Exception as e:
        assert False, f"Unpredicted error: {e}"

    # multiple organisms and ENSG genes
    adata.var.index = gene_map.ENSEMBL_gene[:n_genes]

    adata.write_h5ad(filename=file_path)
    try:
        check_var_index()
    except Exception as e:
        assert False, f"Unpredicted error: {e}"


@pytest.mark.parametrize(
    ("organisms", "var_names", "expected_error"),
    [
        (
            "unsupported organism",
            ["ENSG000001.1", "ENSG000001.2"],
            AnnDataGeneIndexIsNotUnique,
        ),  # duplicate after removing version suffix
        (
            "unsupported organism",
            ["TP53", "TP53"],
            AnnDataGeneIndexIsNotUnique,
        ),  # duplicate gene symbol
        (
            "unsupported organism",
            ["unknown_gene_1", "unknown_gene_2"],
            None,
        ),  # unique genes for a single unsupported organism
        (
            [HomoSapiens.name, "unsupported organism"],
            ["ENSG00000290825", "ENSG00000223972"],
            None,
        ),  # known human genes for multiple organisms
        (
            [HomoSapiens.name, "unsupported organism"],
            ["unknown_gene_1", "unknown_gene_2"],
            AnnDataUnsupportedGenes,
        ),  # genes outside the human gene map for multiple organisms
    ],
)
def test_var_validation_with_unsupported_organisms(organisms, var_names, expected_error):
    adata = ad.AnnData(X=np.eye(len(var_names)))
    adata.var_names = var_names
    adata.obs[ORGANISM_COLUMN] = organisms

    validator = UploadValidator(None)
    validator._multi_exception.raise_on_append = True

    if expected_error:
        with pytest.raises(expected_error):
            validator._check_var_index(adata)
    else:
        assert validator._check_var_index(adata) is None


@pytest.mark.parametrize("set_organism", [False, True, "ont"])
def test_validator(set_organism):
    x = np.eye(10) + 0.1  # not a counts
    adata = ad.AnnData(X=x) 
    if set_organism is True:
        adata.obs[ORGANISM_COLUMN] = HomoSapiens.name
    elif set_organism == "ont":
        adata.obs[ORGANISM_ONT_ID_COLUMN] = HomoSapiens.ontology_id
    
    file_path = TMP_DIR / "bad_adata.h5ad"

    adata.write_h5ad(filename=file_path)
    del adata

    # Wrong X, wrong var, wrong obs, wrong obsm
    validator = UploadValidator(file_path)
    
    try:
        validator.validate()
    except CapMultiException as e:
        if set_organism:
            expected_errors = 4  # gene ids will be validated
        else:
            expected_errors = 3  # gene ids won't be validated
        assert len(e.ex_list) == expected_errors, "Wrong multi exception content!"
    except Exception as e:
        assert False, "Unpredicted exception while the validation!"


def test_df_in_obsm():
    file_path = TMP_DIR / "test_df_in_obsm"
    adata = ad.AnnData(X=np.eye(10))
    adata.obsm["X_df"] = pd.DataFrame(
        data = {
            "x": 1,
            "y": 2,
        },
        index=adata.obs.index,
    )
    adata.write_h5ad(file_path)
    assert adata.obsm["X_df"].shape == (10,2)
    v = UploadValidator(file_path)
    v._multi_exception.raise_on_append = True
    with pytest.raises(AnnDataMissingEmbeddings):
        with read_h5ad(file_path) as adata:
            v._check_obsm(adata)


@pytest.mark.parametrize("with_none", [True, False])
@pytest.mark.parametrize("names_provided", [True, False])
def test_ontology_id_instead_general_metadata(names_provided, with_none):
    file_path = TMP_DIR / "test_ontology_id_instead_general_metadata.h5ad"
    adata = ad.AnnData(X=np.eye(10))

    for col in GENERAL_METADATA:
        ont_id_col = col + "_ontology_term_id"
        adata.obs[ont_id_col] = ont_id_col
        if names_provided:
            adata.obs[col] = col
    
    for col in adata.obs.columns:
        adata.obs[col] = pd.Categorical(adata.obs[col])

    if with_none:
        adata.obs.iloc[5:, :] = None

    adata.write_h5ad(file_path)
    v = UploadValidator(file_path)
    v._multi_exception.raise_on_append = True

    if with_none:
        context = pytest.raises(AnnDataEmptyOrNoneInGeneralMetadata)
    else:
        context = nullcontext()
    
    with read_h5ad(file_path) as adata:
        with context:
            adata.read_obs(GENERAL_METADATA)
            v._check_obs(adata)


def write_adata_with_matrix(path, X, raw_X=None):
    adata = ad.AnnData(X=X)

    if raw_X is not None:
        adata.raw = ad.AnnData(X=raw_X)

    adata.write_h5ad(path)


def test_csc_in_x_raises(tmp_path):
    p = tmp_path / "test.h5ad"

    X = sp.csc_matrix(np.eye(5)) # CSC
    write_adata_with_matrix(p, X=X)

    v = UploadValidator(p)

    with pytest.raises(CSCMatrixInX) as e:
        with read_h5ad(p, edit=False) as cap_adata:
            v._validate_x_and_raw_x_formats(cap_adata)

    assert "X" in e.value.message


def test_csc_in_raw_x_raises(tmp_path):
    p = tmp_path / "test.h5ad"

    X = sp.csr_matrix(np.eye(5)) # valid CSR
    raw_X = sp.csc_matrix(np.eye(5)) # invalid CSC

    write_adata_with_matrix(p, X=X, raw_X=raw_X)

    v = UploadValidator(p)

    with pytest.raises(CSCMatrixInX) as e:
        with read_h5ad(p, edit=False) as cap_adata:
            v._validate_x_and_raw_x_formats(cap_adata)

    assert "raw.X" in e.value.message


def test_csc_in_both_raises(tmp_path):
    p = tmp_path / "test.h5ad"

    X = sp.csc_matrix(np.eye(5))
    raw_X = sp.csc_matrix(np.eye(5))

    write_adata_with_matrix(p, X=X, raw_X=raw_X)

    v = UploadValidator(p)

    with pytest.raises(CSCMatrixInX) as e:
        with read_h5ad(p, edit=False) as cap_adata:
            v._validate_x_and_raw_x_formats(cap_adata)

    assert "`X` and `raw.X`" in e.value.message


def test_dense_and_csr_pass(tmp_path):
    p = tmp_path / "test.h5ad"

    X = np.random.rand(5, 5) # valid dense
    raw_X = sp.csr_matrix(np.eye(5)) # valid CSR

    write_adata_with_matrix(p, X=X, raw_X=raw_X)

    v = UploadValidator(p)

    with read_h5ad(p, edit=False) as cap_adata:
        v._validate_x_and_raw_x_formats(cap_adata) # should not raise


def add_required_obs_columns(adata):
    adata.obs[ORGANISM_COLUMN] = pd.Categorical(
        [HomoSapiens.name] * adata.n_obs
    )
    adata.obs["assay_ontology_term_id"] = pd.Categorical(
        ["EFO:0000001"] * adata.n_obs
    )
    adata.obs["organism_ontology_term_id"] = pd.Categorical(
        [HomoSapiens.ontology_id] * adata.n_obs
    )
    adata.obs["disease_ontology_term_id"] = pd.Categorical(
        ["MONDO:0000001"] * adata.n_obs
    )
    adata.obs["tissue_ontology_term_id"] = pd.Categorical(
        ["UBERON:0000001"] * adata.n_obs
    )


def test_multiple_ontology_ids_raises(tmp_path, monkeypatch):
    file_path = tmp_path / "test_multiple_ids.h5ad"

    adata = ad.AnnData(X=np.eye(3))
    add_required_obs_columns(adata)

    # Add invalid multi-ID category
    adata.obs["tissue_ontology_term_id"] = (
        adata.obs["tissue_ontology_term_id"]
        .cat.add_categories(["UBERON:0001,UBERON:0002"])
    )

    adata.obs.iloc[0, adata.obs.columns.get_loc("tissue_ontology_term_id")] = \
        "UBERON:0001,UBERON:0002"

    adata.write_h5ad(file_path)

    validator = UploadValidator(file_path)
    validator._multi_exception.raise_on_append = True

    # Mock var validation
    def mock_check_var_index(self, cap_adata):
        self._organism = HomoSapiens

    monkeypatch.setattr(
        UploadValidator,
        "_check_var_index",
        mock_check_var_index,
    )

    with read_h5ad(file_path) as cap_adata:
        cap_adata.read_obs()
        validator._check_var_index(cap_adata)

        with pytest.raises(AnnDataMultipleOntologyIDs):
            validator._check_obs(cap_adata)


def test_invalid_disease_prefix_for_human_raises(tmp_path, monkeypatch):
    file_path = tmp_path / "test_invalid_disease_prefix.h5ad"

    adata = ad.AnnData(X=np.eye(3))
    add_required_obs_columns(adata)

    # Replace disease column with invalid prefix
    adata.obs["disease_ontology_term_id"] = pd.Categorical(
        ["DOID:1234"] * adata.n_obs
    )

    adata.write_h5ad(file_path)

    validator = UploadValidator(file_path)
    validator._multi_exception.raise_on_append = True

    # Mock var validation
    def mock_check_var_index(self, cap_adata):
        self._organism = HomoSapiens

    monkeypatch.setattr(
        UploadValidator,
        "_check_var_index",
        mock_check_var_index,
    )

    with read_h5ad(file_path) as cap_adata:
        cap_adata.read_obs()
        validator._check_var_index(cap_adata)

        with pytest.raises(AnnDataInvalidDiseaseOntologyForHuman):
            validator._check_obs(cap_adata)
