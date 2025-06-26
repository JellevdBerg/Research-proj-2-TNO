import pandas as pd

# === Load Datasets ===
df = pd.read_csv("./Data/docs/4.papyrus_present_genes&compounds.csv")
pdb_df = pd.read_csv("./Data/docs/3.pdbligands&genes_combined.csv")
df_raw = pd.read_csv("./Data/docs/1.raw_genes.csv")

# === Classify Compounds ===
df['Active'] = df['pchembl_value_Mean'] > 6.5
df['Inactive'] = df['pchembl_value_Mean'] <= 6.5

# === Summarize by Symbol ===
summary = df.groupby('Symbol').agg(
    Datapoints=('SMILES', 'count'),
    Active=('Active', 'sum'),
    Inactive=('Inactive', 'sum'),
    pchembl_min=('pchembl_value_Mean', 'min'),
    pchembl_max=('pchembl_value_Mean', 'max')
).reset_index()

summary['pchembl_min'] = summary['pchembl_min'].round(3)
summary['pchembl_max'] = summary['pchembl_max'].round(3)
summary['Ratio'] = (summary['Active'] / summary['Datapoints']).round(3)
summary['Activity range'] = summary['pchembl_min'].astype(str) + ' to ' + summary['pchembl_max'].astype(str)

# === Merge PDB Ligand Counts ===
pdb_counts = pdb_df.groupby('Symbol')['pdb_ligand_ID'].nunique().reset_index()
pdb_counts.rename(columns={'pdb_ligand_ID': 'PDB_ligands'}, inplace=True)
summary = summary.merge(pdb_counts, on='Symbol', how='left')
summary['PDB_ligands'] = summary['PDB_ligands'].fillna(0).astype(int)

# === Filter Symbols by Actives/Inactives ===
symbols_less = summary[(summary['Active'] < 30) & (summary['Inactive'] < 30)]['Symbol'].tolist()
symbols_more = summary[(summary['Active'] > 30) & (summary['Inactive'] > 30)]['Symbol'].tolist()

print(f"Number of targets with less than 30 actives or 30 inactives: {len(symbols_less)}")
print(f"Number of Symbols with more than 30 active and inactive: {len(symbols_more)}")

# === Process Sufficient Data ===
filtered_more = df_raw[df_raw['Symbol'].isin(symbols_more)].copy()
filtered_more['PDB-ID Count'] = filtered_more['PDB-ID'].apply(lambda x: len(x.split(',')) if pd.notnull(x) else 0)
filtered_more = filtered_more[filtered_more['PDB-ID Count'] > 0]

final_df = pd.merge(filtered_more, summary, on='Symbol', how='left')
final_df = final_df[final_df['PDB_ligands'] > 0].sort_values(by='Datapoints', ascending=False)

# === Process Insufficient Data ===
filtered_less = df_raw[df_raw['Symbol'].isin(symbols_less)].copy()
filtered_less['PDB-ID Count'] = filtered_less['PDB-ID'].apply(lambda x: len(x.split(',')) if pd.notnull(x) else 0)
filtered_less = filtered_less[filtered_less['PDB-ID Count'] > 0]

final_df_alt = pd.merge(filtered_less, summary, on='Symbol', how='left')
final_df_alt = final_df_alt[final_df_alt['PDB_ligands'] > 0].sort_values(by='Datapoints', ascending=False)

# === Save Output ===
final_df.to_excel("./Data/docs/TOI_sufficient_data.xlsx", index=False)
final_df_alt.to_excel("./Data/docs/TOI_insufficient_data.xlsx", index=False)

print("Saved")
