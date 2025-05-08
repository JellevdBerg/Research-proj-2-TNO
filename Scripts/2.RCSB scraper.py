import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor
import numpy as np

# Load CSV file
file_path = "Data/csv data/1.raw_genes.csv"  # Change this to your actual file path
df = pd.read_csv(file_path)

def get_pdb_ids(uniprot_id):
    """Fetch all PDB IDs from UniProt API using a given UniProt ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            pdb_ids = [
                xref["id"] for xref in data.get("uniProtKBCrossReferences", []) 
                if xref["database"] == "PDB"
            ]
            return ",".join(pdb_ids) if pdb_ids else np.nan  # Return NaN if no PDB IDs found
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
    
    return np.nan  # Return NaN on failure

def update_pdb_id(row):
    """Update the PDB-ID column for a given row."""
    index, uniprot_id = row
    # Always fetch PDB IDs, even if they are already filled in
    pdb_ids = get_pdb_ids(uniprot_id)
    df.at[index, "PDB-ID"] = pdb_ids
    print(f"Updated {uniprot_id} → PDB: {pdb_ids}")

# Use multi-threading for faster API calls
with ThreadPoolExecutor(max_workers=11) as executor:
    executor.map(update_pdb_id, df[["Uniprot"]].dropna().itertuples(index=True))

# Save updated CSV
df.to_csv("Data/csv data/1.raw_genes.csv", index=False)
print("PDB IDs updated!")