import pandas as pd
import json
import os

cluster_number = 7  # Example cluster number, can be changed as needed
model_type = "Scaffold"  # Example model type, can be changed as needed

csv_path = f"./Results/pharmacophores/important features/consensus_weights_cluster_{cluster_number}.csv"
include = {
    "Aromatic": [1],
    "HydrogenAcceptor": [5,6,7,8],
    "HydrogenDonor": [1],
    "Hydrophobic": []
}

# Extract cluster number from filename
output_path = f"./Results/{model_type}_models/{model_type}_cluster_{cluster_number}.json"

# Load & filter
df = pd.read_csv(csv_path)
df = df[df.apply(lambda row: row["cluster"] in include.get(row["name"], []), axis=1)]

# Convert & save
points = df.to_dict(orient="records")
os.makedirs(f"./Results/{model_type}_models", exist_ok=True)
with open(output_path, "w") as f:
    json.dump({"points": points}, f, indent=4)

print(f"Saved to {output_path}")
