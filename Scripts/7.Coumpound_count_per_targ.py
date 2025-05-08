import pandas as pd

# Load the datasets
papyrus_data = pd.read_csv('Data/Papyrus_output/Papyrus_dataset.tsv', sep='\t')
gene_list_raw = pd.read_csv('Data/csv data/1.raw_genes.csv')

# Merge the datasets on a common column (replace 'common_column' with the actual column to merge on)
merged_data = pd.merge(papyrus_data, gene_list_raw[['Symbol', 'Uniprot']], left_on='accession', right_on='Uniprot', how='inner')

# Keep only the specified columns
merged_data = merged_data[['Symbol', 'source', 'CID', 'Uniprot','InChIKey', 'SMILES', 'pchembl_value_Mean']]

# Save to CSV
merged_data.to_csv('Data/csv data/4.papyrus_present_genes&compounds.csv', index=False)