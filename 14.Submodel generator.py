import pandas as pd
import json
import os

csv_path = "./Results/pharmacophores/no_cluster/consensus/no_clusters.csv"
include = {
    "Aromatic": [1,3],
    "HydrogenAcceptor": [10,21],
    "HydrogenDonor": [7],
    "Hydrophobic": []
}

# Extract cluster number from filename
output_path = "./Results/pharmacophores/no_cluster/consensus/val_model.json"

# Load & filter
df = pd.read_csv(csv_path)
df = df[df.apply(lambda row: row["cluster"] in include.get(row["name"], []), axis=1)]

# Convert & save
points = df.to_dict(orient="records")
os.makedirs("./Results/pharmacophores/no_cluster", exist_ok=True)
with open(output_path, "w") as f:
    json.dump({"points": points}, f, indent=4)

print(f"Saved to {output_path}")
