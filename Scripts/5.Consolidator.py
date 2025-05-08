import pandas as pd
import os

input_directory = "Data/PDB ligands"
output_path = "Data/PDB ligands/PDB_ligands_combined.csv"

csv_files = [f for f in os.listdir(input_directory) if f.endswith('.csv')]

df_list = [pd.read_csv(os.path.join(input_directory, file)) for file in csv_files]
combined_df = pd.concat(df_list, ignore_index=True)

combined_df.to_csv(output_path, index=False)
print(f"Combined CSV saved to {output_path}")