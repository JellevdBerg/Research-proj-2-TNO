import pandas as pd

# Load the CSV
df = pd.read_csv("Data/csv data/1.raw_genes.csv")

# Drop rows where 'PDB-ID' is empty or missing
df = df[df["PDB-ID"].notna() & (df["PDB-ID"] != "")]

# Split PDB-IDs and explode the column
df["PDB-ID"] = df["PDB-ID"].str.split(",")
df = df.explode("PDB-ID")

# Save the new CSV
df.to_csv("Data/csv data/2.raw_genes_long.csv", index=False)