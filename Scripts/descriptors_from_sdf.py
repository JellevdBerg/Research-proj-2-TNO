def calc_plip_from_dir(sdf_dir: str, pdb_path: str):
    import warnings, tempfile, os, random
    from rdkit import Chem
    import pandas as pd
    from plip.structure.preparation import PDBComplex

    warnings.filterwarnings("ignore")
    with open(pdb_path) as f:
        pdb_block = f.read()
    prot = Chem.MolFromPDBBlock(pdb_block, removeHs=False)

    dfs = []
    for fname in os.listdir(sdf_dir):
        if not fname.endswith(".sdf"):
            continue
        path = os.path.join(sdf_dir, fname)
        pose_id = os.path.splitext(fname)[0]
        mol = Chem.SDMolSupplier(path, removeHs=False)[0]
        if mol is None:
            continue

        complex = Chem.CombineMols(prot, mol)
        assert complex is not None, "Failed to combine protein and ligand."

        hash = random.getrandbits(128)
        temp_file = os.path.join(tempfile.gettempdir(), f"{hash:032x}_complex.pdb")
        Chem.MolToPDBFile(complex, temp_file, flavor=4)

        mol_plip = PDBComplex()
        mol_plip.load_pdb(temp_file)
        mol_plip.analyze()

        longnames = [x.longname for x in mol_plip.ligands]
        bsids = [":".join([x.hetid, x.chain, str(x.position)]) for x in mol_plip.ligands]
        indices = [j for j, x in enumerate(longnames) if x == "UNL"]
        if not indices:
            continue
        bsid = bsids[indices[0]]

        interactions = mol_plip.interaction_sets[bsid].all_itypes
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
