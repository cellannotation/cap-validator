import pytest
import numpy as np
import anndata as ad
import scipy.sparse as sp
from pathlib import Path
import tempfile
from cap_anndata import CapAnnDataDF, read_h5ad

from cap_validator.upload_validator import UploadValidator, OBS_COLUMNS_REQUIRED, ORGANISM_COLUMN
from cap_validator.gene_mapping import GeneMap, EnsemblOrganism
from cap_validator.errors import (
    AnnDataMissingEmbeddings,
    AnnDataMisingObsColumns,
    AnnDataNonStandardVarError,
    CapMultiException
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
        v.multi_exception.raise_on_append = True
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
    for col in OBS_COLUMNS_REQUIRED:
        adata.obs[col] = "test_value"
    
    adata.write_h5ad(file_path)
    
    with read_h5ad(file_path, edit=False) as cap_adata:
        v = UploadValidator(adata_path=file_path)
        v.multi_exception.raise_on_append = True
        cap_adata.read_obs()
        df = cap_adata.obs.copy()

        def check_obs(ca, correct_expected: bool):
            try:
                v._check_obs(ca)
                if not correct_expected:
                    assert False, "Must not be correct obs!"
            except AnnDataMisingObsColumns:
                assert not correct_expected, "Unexpected result"

        check_obs(cap_adata, True)

        for col in OBS_COLUMNS_REQUIRED:
            cap_adata.obs = CapAnnDataDF.from_df(df.drop(col, axis=1, inplace=False))
            check_obs(cap_adata, False)

        cap_adata.obs = CapAnnDataDF.from_df(df[[]])
        check_obs(cap_adata, False)


def test_var_index():
    v = UploadValidator(None)
    gene_map = GeneMap.data_frame()
    n_rows = 100
    n_genes = 100
    X = np.ones(shape=(n_rows, n_genes), dtype=np.float32)
    adata = ad.AnnData(X=X)
    
    file_path = TMP_DIR / "test_var_index.h5ad"

    def check_var_index():
        with read_h5ad(file_path, edit=False) as cap_adata:
            cap_adata.read_var()
            cap_adata.read_obs(columns=[ORGANISM_COLUMN])
            v._check_var_index(cap_adata)

    adata.write_h5ad(filename=file_path)
    
    adata.obs[ORGANISM_COLUMN] = EnsemblOrganism.HUMAN.value
    adata.var.index = gene_map.ENSEMBL_gene[:n_genes]
    adata.write_h5ad(filename=file_path)
    check_var_index()

    adata.var.index = [f"{g}.{i}" for i, g in enumerate(gene_map.ENSEMBL_gene[:n_genes])]
    adata.write_h5ad(filename=file_path)
    check_var_index()

    adata.var.index = map(str, range(n_genes))
    adata.write_h5ad(filename=file_path)
    try:
        check_var_index()
    except AnnDataNonStandardVarError:
        pass
    
    adata.obs.organism = 'unsuported'
    adata.write_h5ad(filename=file_path)
    try:
        check_var_index()
    except:
        assert False, "Wrong validation failure for unsupported organism!"

    adata.obs['organism'] = adata.obs['organism'].cat.add_categories(['new organism'])
    adata.obs.loc['0', 'organism'] = 'new organism'

    adata.write_h5ad(filename=file_path)
    try:
        check_var_index()
    except AnnDataNonStandardVarError:
        pass

    adata.var.index = gene_map.ENSEMBL_gene[:n_genes]
    adata.write_h5ad(filename=file_path)
    try:
        check_var_index()
    except Exception as e:
        assert False, f"Unexpected error: {e}"

@pytest.mark.parametrize("set_organism", [False, True])
def test_validator(set_organism):
    file_path = TMP_DIR / "test_validator.h5ad"

    adata = ad.AnnData(X=np.eye(10).astype(np.float32))
    if set_organism:
        adata.obs[ORGANISM_COLUMN] = EnsemblOrganism.HUMAN.value
    
    adata.write_h5ad(filename=file_path)
    del adata

    validator = UploadValidator(file_path)
    
    try:
        validator.validate()
    except CapMultiException as e:
        expected_errors = 4 if set_organism else 3
        assert len(e.ex_list) == expected_errors, "Wrong multi exception content!"
    except Exception as e:
        assert False, "Unexpected exception while validation!"
