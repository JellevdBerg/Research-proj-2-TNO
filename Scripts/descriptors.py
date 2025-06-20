import os
import random
import tempfile
from typing import Any, Iterable

import numpy as np
import pandas as pd
from plip.structure.preparation import PDBComplex
from qsprpred.data.descriptors.sets import DescriptorSet
from rdkit import Chem
from rdkit.Chem import Mol

from spock.storage.tabular import SpockProtein
#from spock.utils.formats.pdb import read_pdb_block

def read_pdb_block(pdb_block):
    prot = Chem.MolFromPDBBlock(
        pdb_block, sanitize=True, removeHs=False, proximityBonding=True
    )
    assert prot is not None, "Could not read protein from PDB block: \n" + pdb_block
    #prot = protein_sanitizer(prot)
    return prot

def calc_plip(mols: Iterable[tuple[str, Chem.Mol]], pdb_block: str):
    import warnings

    warnings.filterwarnings("ignore")
    prot = read_pdb_block(pdb_block)
    dfs = []
    for pose_id, mol in mols:
        complex = Chem.CombineMols(prot, mol)
        assert complex is not None, "Failed to combine protein and ligand."
        hash = random.getrandbits(128)
        temp_file = "%032x_complex.pdb" % hash
        temp_file = os.path.join(tempfile.gettempdir(), temp_file)
        Chem.MolToPDBFile(complex, temp_file, flavor=4)
        mol = PDBComplex()
        mol.load_pdb(temp_file)
        mol.analyze()
        longnames = [x.longname for x in mol.ligands]
        bsids = [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands]
        indices = [j for j, x in enumerate(longnames) if x == "UNL"]
        bsid = bsids[indices[0]]
        interactions = mol.interaction_sets[bsid].all_itypes
        df_mol = pd.DataFrame()
        for interaction in interactions:
            name = interaction.__class__.__name__.replace("_interaction", "")
            if hasattr(interaction, "protisdon"):
                donor = "d" if interaction.protisdon else "a"
            else:
                donor = ""
            res_name = f"{name}{donor}_{interaction.restype}_{interaction.resnr}_{interaction.reschain}"
            df_mol[res_name] = True
        df_mol.loc[0, df_mol.columns] = True
        df_mol["Pose_ID"] = pose_id
        dfs.append(df_mol)
        os.remove(temp_file)
    df = pd.concat(dfs)
    df.fillna(False, inplace=True)
    warnings.filterwarnings("default")
    return df


class PLIPIFP(DescriptorSet):

    def __init__(
        self,
        protein: SpockProtein,
        n_poses: int = 1,
        id_prop: str | None = None,
    ):
        super().__init__(id_prop=id_prop)
        self.protein = protein
        self.nPoses = n_poses
        self._descriptors = []

    @property
    def descriptors(self) -> list[str]:
        return self._descriptors

    @descriptors.setter
    def descriptors(self, value: list[str]):
        self._descriptors = value

    def __str__(self):
        return "PLIPIFP"

    def parsePropsAndMols(
        self, mols: list[str | Mol], props: dict[str, list[Any]] | None
    ) -> tuple[list[Mol], dict[str, list[Any]]]:
        return mols, props

    def getDescriptors(
        self, mols: list[Mol], props: dict[str, list[Any]], *args, **kwargs
    ) -> np.ndarray:
        to_calculate = []
        parents = []
        mol_ids = []
        for mol in mols:
            for rep in mol.representations:
                assert hasattr(rep, "poses"), (
                    f"No poses attribute found in representation: {rep}. "
                    f"Make sure you are working with a 'DockableStore'."
                )
                poses = rep.poses[self.protein.id]
                poses = [(pose.id, pose.as_rd_mol()) for pose in poses[: self.nPoses]]
                to_calculate.extend(poses)
                parents.extend([rep.id] * len(poses))
                mol_ids.extend([mol.id] * len(poses))
        df = calc_plip(
            to_calculate,
            self.protein.as_pdb(),
        )
        df["Parent_ID"] = parents
        df[self.idProp] = mol_ids
        # df = df.groupby("Pose_ID").sum()
        self.descriptors = self.descriptors + df.columns.tolist()
        return df.values