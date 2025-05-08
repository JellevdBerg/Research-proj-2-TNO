import pandas as pd

# Load the files
df1 = pd.read_csv("Data/csv data/2.raw_genes_long.csv")  # Contains Symbol, Uniprot, PDB-ID
df2 = pd.read_csv("Data/PDB ligands/PDB_ligands_combined.csv")  # Contains @structureId (PDB-ID) and @chemicalID

# Merge file2 into file1 based on PDB-ID
merged = df2.merge(df1, left_on="@structureId", right_on="PDB-ID", how="outer")

# Rename @chemicalID to PDB_ligand
merged = merged.rename(columns={"@chemicalID": "pdb_ligand_ID"})

# Keep only the desired columns
merged = merged[["Symbol", "Uniprot", "PDB-ID", "pdb_ligand_ID"]]

# Save the result
merged.to_csv("Data/csv data/3.pdbligands&genes_combined.csv", index=False)