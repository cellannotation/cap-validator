from enum import Enum
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
HUMAN_GENE_MAP_PATH = HERE / "data/homo_sapiens.csv"
MOUSE_GENE_MAP_PATH = HERE / "data/mus_musculus.csv"


class EnsemblOrganism(Enum):
    HUMAN = "Homo sapiens"
    MOUSE = "Mus musculus"
    MULTI_SPECIES = "Multi species"

    @staticmethod
    def from_str(value):
        if value == EnsemblOrganism.HUMAN.value:
            return EnsemblOrganism.HUMAN
        if value == EnsemblOrganism.MOUSE.value:
            return EnsemblOrganism.MOUSE
        if value == EnsemblOrganism.MULTI_SPECIES.value:
            return EnsemblOrganism.MULTI_SPECIES
    
    def prefix(self):
        if self.value in [
            EnsemblOrganism.HUMAN.value,
            EnsemblOrganism.MULTI_SPECIES.value,
        ]:
            return "ENSG"
        if self.value == EnsemblOrganism.MOUSE.value:
            return "ENSMUSG"

    def map_file_name(self):
        if self.value in [
            EnsemblOrganism.HUMAN.value,
            EnsemblOrganism.MULTI_SPECIES.value,
        ]:
            return HUMAN_GENE_MAP_PATH
        if self.value == EnsemblOrganism.MOUSE.value:
            return MOUSE_GENE_MAP_PATH

    @staticmethod
    def supported(organism: str) -> bool:
        return organism in (o.value for o in EnsemblOrganism)


class GeneMap:

    @staticmethod
    def data_frame(organisms=None, index_col=None):
        if organisms is None:
            organisms = [EnsemblOrganism.HUMAN.value, EnsemblOrganism.MOUSE.value]
        
        if isinstance(organisms, str):
            # the single string value is given
            organisms = [organisms]

        dfs = []
        for organism in organisms:
            eo = EnsemblOrganism.from_str(organism)
            if eo is not None:
                fp = eo.map_file_name()
                df = pd.read_csv(fp, sep=',', header=0, index_col=index_col)  # index=0 to make Ensemble ids index
                dfs.append(df)
        if len(dfs) > 0:
            return pd.concat(dfs, axis=0)
        else:
            return None
