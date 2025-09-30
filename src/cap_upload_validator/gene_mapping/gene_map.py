from dataclasses import dataclass
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
HUMAN_GENE_MAP_PATH = HERE / "data/homo_sapiens.csv"
MOUSE_GENE_MAP_PATH = HERE / "data/mus_musculus.csv"


@dataclass(frozen=True)
class _BasicOrganism:
    name: str
    ontology_id: str | None
    gene_prefix: str | None
    gene_map_path: str | None

@dataclass(frozen=True)
class HomoSapiens(_BasicOrganism):
    name = "Homo sapiens"
    ontology_id = "NCBITaxon:9606"
    gene_prefix = "ENSG"
    gene_map_path = HUMAN_GENE_MAP_PATH

@dataclass(frozen=True)
class MusMusculus(_BasicOrganism):
    name = "Mus musculus"
    ontology_id = "NCBITaxon:10090"
    gene_prefix = "ENSMUSG"
    gene_map_path = MOUSE_GENE_MAP_PATH

@dataclass(frozen=True)
class MultiSpecies(HomoSapiens):
    name = "Multi species"
    ontology_id = None


def str_to_organism(organism_str: str) -> _BasicOrganism | None:
    clean_str = organism_str.strip().lower()
    if clean_str == "homo sapiens":
        return HomoSapiens
    if clean_str == "mus musculus":
        return MusMusculus
    if clean_str == "multi species":
        return MultiSpecies
    return None

class GeneMap:

    @staticmethod
    def data_frame(organisms: str | _BasicOrganism = None, index_col=None):
        if organisms is None:
            organisms = [HomoSapiens, MusMusculus]
        
        if isinstance(organisms, str):
            # the single string value is given
            organisms = str_to_organism(organisms)

        dfs = []
        for organism in organisms:
            if isinstance(organism, _BasicOrganism):
                fp = organism.gene_map_path
                df = pd.read_csv(fp, sep=',', header=0, index_col=index_col)  # index=0 to make Ensemble ids index
                dfs.append(df)
        if len(dfs) > 0:
            return pd.concat(dfs, axis=0)
        else:
            return None
