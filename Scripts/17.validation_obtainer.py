import os
import pandas as pd

data = []
version = "PLIP" # PLIP or Scaffold
folder = f'./Results/pharmacophores/{version}_based/{version}_results'

for file in os.listdir(folder):
    if not file.startswith('val_cluster_'):
        continue
    cluster = f"cluster_{file.split('_')[-1]}"
    with open(os.path.join(folder, file)) as f:
        content = f.read()

    blocks = content.split('$$$$')
    for block in blocks:
        block = block.strip()
        if not block:
            continue
        lines = block.splitlines()
        name = lines[0].strip()
        if '_-' in name:
            name = name.split('_-')[0] + '_VINA'

        affinity = None
        rmsd = None
        for i, line in enumerate(lines):
            if line.startswith('>  <minimizedAffinity>'):
                affinity = lines[i+1].strip()
            if line.startswith('>  <minimizedRMSD>'):
                rmsd = lines[i+1].strip()
        data.append({
            'cluster': cluster,
            'filename': name,
            'affinity': affinity,
            'rmsd': rmsd
        })

df = pd.DataFrame(data)
df.to_csv(f'./Results/pharmacophores/{version}_based/{version}_results/merged_output.csv', index=False)
