import pandas as pd
import numpy as np
from scipy.sparse import issparse
from cap_anndata import CapAnnData, read_h5ad
import logging
from h5py import Dataset, Group

from .gene_mapping import (
    GeneMap,
    HomoSapiens,
    MusMusculus,
    MultiSpecies,
    Organism,
    str_to_organism,
    ontology_id_to_organism,
)
from .errors import (
    CapMultiException,
    AnnDataMissingCountMatrix,
    AnnDataInvalidCountMatrix,
    AnnDataMissingEmbeddings,
    AnnDataMissingObs,
    AnnDataMissingObsColumns,
    AnnDataMultipleDiseaseOntologyIDs,
    AnnDataInvalidDiseaseOntologyForHuman,
    AnnDataMissingVarIndex,
    AnnDataNumericVarIndex,
    AnnDataVarNotSubsetOfRawVar,
    AnnDataGeneIndexIsNotUnique,
    AnnDataUnsupportedGenes,
    BadAnnDataFile,
    AnnDataEmptyOrNoneInGeneralMetadata,
    CSCMatrixInX,
)
from typing import Optional

logger = logging.getLogger(__name__)

MAX_OBS_ROWS_TO_CHECK = 100
EMBEDDING_PREFIX = "X_"
ORGANISM_COLUMN = "organism"
ORGANISM_ONT_ID_COLUMN = f"{ORGANISM_COLUMN}_ontology_term_id"
GENERAL_METADATA = ["assay", "disease", ORGANISM_COLUMN, "tissue"]
DISEASE_ONTOLOGY_HUMAN_PREFIXES = ("MONDO", "PATO")

class UploadValidator:

    """The class perform the validation """
    def __init__(self, adata_path: str) -> None:
        self._adata_path = adata_path
        # Create container to raise a multiple errors for client
        self._multi_exception = CapMultiException()    
        self._organism: Organism = None  
        self._ensembl_ids = None

    @property
    def adata_path(self) -> str:
        return self._adata_path
    
    @property
    def organism(self) -> str:
        return self._organism

    @property
    def ensembl_ids(self) -> pd.Index:
        return self._ensembl_ids

    def validate(self, report_success: bool = True) -> None:
        """The method validates the input AnnData on the upload stage and raises exceptions if something wrong."""
        logger.debug("Begin anndata file validation...")

        if not str(self._adata_path).endswith(".h5ad"):
            raise BadAnnDataFile()
        
        with read_h5ad(self._adata_path, edit=False) as cap_adata:
            cap_adata.read_obs(columns=GENERAL_METADATA)  # TODO: read all columns?
            cap_adata.read_var(columns=[])
            if cap_adata.raw is not None:
                cap_adata.raw.read_var(columns=[])

            self._validate_x_and_raw_x_formats(cap_adata)
            self._check_X(cap_adata)
            self._check_obsm(cap_adata)
            self._check_var_index(cap_adata)
            self._check_obs(cap_adata) # Must be called after organism detection in _check_var_index 

        # Check any errors were during validation stage and raise them
        if self._multi_exception.have_errors():
            raise self._multi_exception
        if report_success:
            print("Validation passed!")
        logger.debug("Finish anndata file validation")

    def find_missing_genes(self) -> Optional[pd.DataFrame]:
        """The method finds missing genes from the gene map for the validated AnnData."""
        logger.debug("Begin finding missing genes...")
        
        if not str(self._adata_path).endswith(".h5ad"):
            raise BadAnnDataFile()
        
        missing_genes = None
        with read_h5ad(self._adata_path, edit=False) as cap_adata:
            cap_adata.read_obs(columns=GENERAL_METADATA)  # TODO: read all columns?
            cap_adata.read_var()
            if cap_adata.raw is not None:
                cap_adata.raw.read_var()
            
            missing_genes_mask = self._check_var_index(cap_adata=cap_adata)
            if missing_genes_mask is not None:
                missing_genes = cap_adata.var.loc[missing_genes_mask]
        logger.debug("Finished finding missing genes!")
        return missing_genes
        

    def _check_X(self, cap_adata: CapAnnData) -> None:
        logger.debug("Begin checking X")
        X = cap_adata.raw.X if cap_adata.raw is not None else cap_adata.X

        if X is None:
            self._multi_exception.append(AnnDataMissingCountMatrix())
            return

        if not self._check_is_positive_integers(cap_adata):
            self._multi_exception.append(AnnDataInvalidCountMatrix())

        logger.debug("Finish checking X")

    @staticmethod
    def has_only_integers(arr: np.ndarray) -> bool:
        return np.all((arr - arr.astype(int)) == 0)

    @staticmethod
    def has_negative_values(arr: np.ndarray) -> bool:
        return np.any(arr < 0)

    def _check_is_positive_integers(self, cap_adata: CapAnnData) -> bool:
        n_cells = cap_adata.shape[0]
        max_rows = min(n_cells, MAX_OBS_ROWS_TO_CHECK)

        X = cap_adata.raw.X if cap_adata.raw is not None else cap_adata.X

        is_sparse = issparse(X[0])
        logger.debug(f"Checking elements in X to be positive integers, with n_cells = {n_cells}, max_rows = {max_rows}, is_sparse = {is_sparse}!")

        arr = X[0:max_rows].data if is_sparse else X[0:max_rows]
        if not self.has_only_integers(arr) or self.has_negative_values(arr):
            logger.debug(f"There are not positive integers found in X!")
            return False

        return True

    def _check_obsm(self, cap_adata: CapAnnData) -> None:
        logger.debug("Begin checking obsm")
        if self._has_embeddings(cap_adata) is False:
            self._multi_exception.append(AnnDataMissingEmbeddings())
        logger.debug("Finished checking obsm!")

    def _has_embeddings(self, cap_adata: CapAnnData) -> bool:
        if cap_adata.obsm is None:
            logger.debug("Obsm is not found in anndata!")
            return False
        
        n_cells = cap_adata.shape[0]

        for field in cap_adata.obsm_keys():
            if field.startswith(EMBEDDING_PREFIX):
                entity = cap_adata.obsm[field]
                if isinstance(entity, Dataset) and entity.shape == (n_cells, 2):
                    # looking for dense matrix of N x 2 shape 
                    return True

        logger.debug(f"Embeddings not found in obsm_keys = {cap_adata.obsm_keys()}!")
        return False

    def _check_obs(self, cap_adata: CapAnnData) -> None:
        logger.debug("Start checking obs")

        obs_keys = list(cap_adata.obs_keys())
        logger.debug(f"Checking obs_columns = {obs_keys} for required {GENERAL_METADATA}!")

        if cap_adata.obs is None or not obs_keys:
            logger.debug(".obs is missing completely.")
            self._multi_exception.append(AnnDataMissingObs())
            return

        missing_columns: list[str] = []
        none_columns: set[str] = set()
        empty_columns: set[str] = set()

        def _check_column(series: pd.Series, name: str):
            has_none, has_empty = self._classify_missing(series)
            if has_none:
                none_columns.add(name)
            if has_empty:
                empty_columns.add(name)

        for col in GENERAL_METADATA:
            ont_id_col = f"{col}_ontology_term_id"

            col_in_obs = col in obs_keys
            ont_id_col_in_obs = ont_id_col in obs_keys

            # Column missing entirely
            if not (col_in_obs or ont_id_col_in_obs):
                missing_columns.append(col)
                continue

            # Validate regular column values
            if col_in_obs:
                _check_column(cap_adata.obs[col], col)

            # Validate ontology column values
            if ont_id_col_in_obs:
                if ont_id_col not in cap_adata.obs.columns:
                    cap_adata.read_obs(columns=[ont_id_col])

                _check_column(cap_adata.obs[ont_id_col], ont_id_col)

                if col == "disease":
                    self._validate_disease_ontology(cap_adata.obs[ont_id_col])

        # Report missing columns
        if missing_columns:
            logger.debug("Missing required obs columns: " + ", ".join(missing_columns))
            self._multi_exception.append(AnnDataMissingObsColumns(missing_columns=missing_columns))

        # Report empty/None columns
        if none_columns or empty_columns:
            logger.debug(
                "Required obs columns contain empty: %s or None: %s values.",
                ", ".join(sorted(empty_columns)) if empty_columns else "—",
                ", ".join(sorted(none_columns)) if none_columns else "—",
            )
            self._multi_exception.append(
                AnnDataEmptyOrNoneInGeneralMetadata(
                    none_columns=list(none_columns),
                    empty_columns=list(empty_columns),
                )
            )

        logger.debug("Finished checking obs!")

    def _validate_disease_ontology(self, series: pd.Series) -> None:
        if series is None:
            return

        has_multiple_ids = False
        has_invalid_prefix_for_human = False

        for value in series.dropna():
            value_str = str(value).strip()

            if not value_str:
                continue

            # Multiple IDs restriction
            if "," in value_str:
                has_multiple_ids = True
                continue

            # Human-specific restriction
            if self._organism is HomoSapiens:
                delimiter = ":"
                if delimiter not in value_str:
                    continue  # format validation not in scope

                prefix = value_str.split(delimiter, 1)[0]
                if prefix not in DISEASE_ONTOLOGY_HUMAN_PREFIXES:
                    has_invalid_prefix_for_human = True

        # Append errors only once
        if has_multiple_ids:
            self._multi_exception.append(AnnDataMultipleDiseaseOntologyIDs())

        if has_invalid_prefix_for_human:
            self._multi_exception.append(AnnDataInvalidDiseaseOntologyForHuman())

    @staticmethod
    def _classify_missing(series: pd.Series) -> tuple[bool, bool]:
        """
        Returns:
            has_none  -> True if None / NaN values are present
            has_empty -> True if empty or whitespace-only strings are present
        """
        has_none = series.isna().any()

        # Only check empties on non-null values
        non_null = series.dropna()
        has_empty = non_null.astype(str).str.match(r"^\s*$").any()

        return has_none, has_empty

    def _check_var_index(self, cap_adata: CapAnnData) -> Optional[pd.Series]:
        logger.debug("Start checking var index...")

        index = cap_adata.var.index

        if index is None or index.empty:
            self._multi_exception.append(AnnDataMissingVarIndex())
            return

        if pd.api.types.is_any_real_numeric_dtype(index):
            self._multi_exception.append(AnnDataNumericVarIndex())
            return

        clean_index = self._remove_gene_version(index)
        self._ensembl_ids = clean_index

        if not clean_index.is_unique:
            logger.debug("There are non unique gene ids in .var.index!")
            self._multi_exception.append(AnnDataGeneIndexIsNotUnique())
            return

        # Check if the var.index is a subset of raw.var.index
        if cap_adata.raw is not None and cap_adata.raw.var is not None:
            logger.debug("As of raw exists, checking that var.index is a subset of raw.var.index!")
            if not index.isin(cap_adata.raw.var.index).all():
                self._multi_exception.append(AnnDataVarNotSubsetOfRawVar())
                return

        # Check the number of organisms in the dataset
        known_organisms = {HomoSapiens.name, MusMusculus.name} # Only Human and Mouse supported this moment
        obs_keys = cap_adata.obs_keys()

        if ORGANISM_COLUMN in obs_keys:
            dataset_organisms = cap_adata.obs[ORGANISM_COLUMN].unique().tolist()
            dataset_organisms = [o for o in dataset_organisms if o]
            dataset_organisms = list(map(str_to_organism, dataset_organisms))

        elif ORGANISM_ONT_ID_COLUMN in obs_keys:
            if ORGANISM_ONT_ID_COLUMN not in cap_adata.obs.columns:
                cap_adata.read_obs(columns=[ORGANISM_ONT_ID_COLUMN])

            org_ids = cap_adata.obs[ORGANISM_ONT_ID_COLUMN].unique().tolist()
            org_ids = [o for o in org_ids if o]
            dataset_organisms = list(map(ontology_id_to_organism, org_ids))

        else:
            dataset_organisms = []

        logger.debug(f"Organism(s) in dataset = {dataset_organisms}, known organisms = {known_organisms}")

        missing_genes_mask = None
        # Check ENSEMBL ids for supported organism
        if len(dataset_organisms) == 1:
            organism = dataset_organisms[0]
            self._organism = organism

            if organism.name in known_organisms:
                logger.debug("Single known organism found, validating gene IDs.")
                missing_genes_mask = self._validate_gene_ids(clean_index, organism)
            else:
                logger.debug("Unknown organism found, skipping gene validation.")
        elif len(dataset_organisms) > 1:
            logger.debug("There are multiple organisms in dataset")
            self._organism = MultiSpecies
            missing_genes_mask = self._validate_gene_ids(clean_index, MultiSpecies)

        logger.debug("Finished checking var index!")
        return missing_genes_mask

    def _validate_gene_ids(
        self,
        ens_ids: pd.Index,
        organism: Organism,
    ) -> Optional[pd.Series]:
        """
        The method finds missing genes from gene map for given organism. 
        Return None if all genes are valid. 
        Else return pd.Series of boolean mask of missing genes.
        """
        if ens_ids.empty:
            self._multi_exception.append(AnnDataMissingVarIndex())
            logger.debug("Gene names are missed!")
            return

        if pd.api.types.is_any_real_numeric_dtype(ens_ids):
            logger.debug("Gene names are numeric!")
            self._multi_exception.append(AnnDataNumericVarIndex())
            return

        df = GeneMap.data_frame(organisms=organism)
        missing_mask = ~ens_ids.isin(df["ENSEMBL_gene"])

        if missing_mask.any():
            missing_genes_count = missing_mask.sum()
            logger.debug(f"{missing_genes_count} gene(s) are not standard!")
            self._multi_exception.append(AnnDataUnsupportedGenes(missing_genes_count=missing_genes_count))
            return missing_mask

        return None

    @staticmethod
    def _remove_gene_version(ensemble_ids: pd.Index) -> pd.Index:
        """
        The static methods removes the version suffixes from the ensemble ids.
        Examples:
            ENSG0001.8 -> ENSG0001
            ENSG0001   -> ENSG0001
        """
        clean_index = ensemble_ids.to_series().apply(lambda x: x.split(".")[0])
        clean_index = pd.Index(clean_index)
        return clean_index

    def _is_csc(self, group_or_dataset) -> bool:
        """
        Returns True if HDF5 object represents a CSC sparse matrix.
        """
        if not isinstance(group_or_dataset, Group):
            return False

        encoding = group_or_dataset.attrs.get("encoding-type", None)
        return encoding == "csc_matrix"

    def _is_csr(self, group_or_dataset) -> bool:
        if not isinstance(group_or_dataset, Group):
            return False
        return group_or_dataset.attrs.get("encoding-type", None) == "csr_matrix"

    def _is_dense(self, group_or_dataset) -> bool:
        return isinstance(group_or_dataset, Dataset)

    def _validate_x_and_raw_x_formats(self, cap_adata: CapAnnData) -> None:
        """
        Validate that X and raw.X (if exists) are dense or CSR.
        Raise CSCMatrixInX otherwise.
        """
        locations = []
    
        f = cap_adata.file
    
        # X
        x = f["X"]
        if self._is_csc(x):
            locations.append("X")
        elif not (self._is_dense(x) or self._is_csr(x)):
            locations.append("X")

        # raw.X
        if "raw" in f and "X" in f["raw"]:
            raw_x = f["raw/X"]
            if self._is_csc(raw_x):
                locations.append("raw.X")
            elif not (self._is_dense(raw_x) or self._is_csr(raw_x)):
                locations.append("raw.X")

        if locations:
            raise CSCMatrixInX(locations=locations)
