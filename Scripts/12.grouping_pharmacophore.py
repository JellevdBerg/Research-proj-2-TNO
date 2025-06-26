# store pharmacophores based on clusters
import os
import shutil
import pandas as pd

df_clusters = pd.read_csv("./Data/docs/11.PLIPIFP_clusters.csv")  # replace with actual path

# Directory containing the JSON files
json_dir = "./Results/pharmacophores"  # replace with actual path

# Ensure cluster folders exist and move files
for _, row in df_clusters.iterrows():
    cluster = row['cluster']
    name_prefix = row['name'][:4]  # first 3-4 characters
    cluster_folder = os.path.join(json_dir, f"cluster_{cluster}")
    os.makedirs(cluster_folder, exist_ok=True)

    # Find and move matching json files
    for file in os.listdir(json_dir):
        if file.endswith(".json") and file.startswith(name_prefix):
            src = os.path.join(json_dir, file)
            dst = os.path.join(cluster_folder, file)
            shutil.move(src, dst)