from dataclasses import dataclass
import pandas as pd
from pathlib import Path
from typing import Union

HERE = Path(__file__).parent
HUMAN_GENE_MAP_PATH = HERE / "data/homo_sapiens.csv"
MOUSE_GENE_MAP_PATH = HERE / "data/mus_musculus.csv"


@dataclass(frozen=True)
class Organism:
    name: str
    ontology_id: Union[str, None]
    gene_prefix: Union[str, None]
    gene_map_path: Union[str, None]


@dataclass(frozen=True)
class HomoSapiens(Organism):
    name = "Homo sapiens"
    ontology_id = "NCBITaxon:9606"
    gene_prefix = "ENSG"
    gene_map_path = HUMAN_GENE_MAP_PATH

@dataclass(frozen=True)
class MusMusculus(Organism):
    name = "Mus musculus"
    ontology_id = "NCBITaxon:10090"
    gene_prefix = "ENSMUSG"
    gene_map_path = MOUSE_GENE_MAP_PATH

@dataclass(frozen=True)
class MultiSpecies(HomoSapiens):
    name = "Multi species"
    ontology_id = None

@dataclass(frozen=True)
class UnsupportedOrganism(Organism):
    name = "Unsupported"
    ontology_id = None
    gene_prefix = None
    gene_map_path = None


def str_to_organism(organism_str: str) -> Organism:
    clean_str = organism_str.strip().lower()
    if clean_str == "homo sapiens":
        return HomoSapiens
    if clean_str == "mus musculus":
        return MusMusculus
    if clean_str == "multi species":
        return MultiSpecies
    return UnsupportedOrganism


def ontology_id_to_organism(ont_id_str: str) -> Organism:
    for organism in [HomoSapiens, MusMusculus]:
        if ont_id_str == organism.ontology_id:
            return organism
    return UnsupportedOrganism


class GeneMap:

    @staticmethod
    def data_frame(
        organisms: Union[str, Organism] = None,
        index_col=None,
    ):
        if organisms is None:
            organisms = [HomoSapiens, MusMusculus]
        
        if isinstance(organisms, str):
            # the single string value is given
            organisms = str_to_organism(organisms)
        if not isinstance(organisms, list):
            organisms = [organisms]
        
        dfs = []
        for organism in organisms:
            if issubclass(organism, Organism):
                fp = organism.gene_map_path
                if fp is not None:
                    df = pd.read_csv(fp, sep=',', header=0, index_col=index_col)  # index=0 to make Ensemble ids index
                    dfs.append(df)
        if len(dfs) > 0:
            return pd.concat(dfs, axis=0)
        else:
            return None
