import os
import platform
import pandas as pd
import matplotlib.pyplot as plt
from pymol import cmd
from conphar.Pharmacophores import parse_json_pharmacophore, save_pharmacophore_to_pymol, save_pharmacophore_to_json, compute_concensus_pharmacophore, get_ligand_receptor_pharmacophore

# --- Paths ---
receptor_path = "./Data/target/AF-Q13705-ACVR2B_-_prepared.pdb"
ligand_dir = "./Data/other/SDF poses/"
results_dir = "./Results/pharmacophores"
os.makedirs(results_dir, exist_ok=True)

# --- Generate pharmacophores (Linux only) ---
if platform.system() == "Linux":
    for file in os.listdir(ligand_dir):
        if file.endswith(".sdf"):
            ligand_path = os.path.join(ligand_dir, file)
            output_name = os.path.splitext(file)[0]
            out_path = os.path.join(results_dir, output_name)
            get_ligand_receptor_pharmacophore(receptor_path, ligand_path, out_path)
else:
    print("Skipping pharmacophore generation — not on Linux!")

# --- Iterate over clusters ---
for cluster in sorted(os.listdir(results_dir)):
    cluster_path = os.path.join(results_dir, cluster)
    if not os.path.isdir(cluster_path) or "cluster" not in cluster.lower():
        continue

    print(f"Processing {cluster}...")

    # --- Parse pharmacophore data ---
    p4_table = pd.DataFrame()
    for file in os.listdir(cluster_path):
        if file.endswith(".json"):
            try:
                p4, lig, rec = parse_json_pharmacophore(os.path.join(cluster_path, file))
                p4['ligand'] = file.replace('.json', '')
                p4_table = pd.concat([p4_table, p4], ignore_index=True)
            except Exception:
                pass

    if p4_table.empty:
        print(f"No pharmacophore data found in {cluster}")
        continue

    # --- Clean & save all features ---
    p4_table['color'] = p4_table['color'].replace({'navy': 'blue', 'white': 'yellow'})
    p4_table = p4_table[~p4_table['name'].isin(['NegativeIon', 'PositiveIon'])]

    cons_path = os.path.join(cluster_path, "consensus")
    os.makedirs(cons_path, exist_ok=True)

    save_pharmacophore_to_pymol(p4_table, f"{cons_path}/all_features_{cluster}.pse")
    save_pharmacophore_to_json(p4_table, f"{cons_path}/all_features_{cluster}.json")

    # --- Compute consensus ---
    consensus, links = compute_concensus_pharmacophore(
        p4_table,
        save_data_per_descriptor=True,
        out_folder=cons_path,
        cmap_plots="viridis",
        method='complete',
        h_dist=2.8
    )

    os.makedirs(f'{results_dir}/important features', exist_ok=True)
    save_pharmacophore_to_pymol(consensus, f'{cons_path}/{cluster}_concensus.pse', select='concensus')
    save_pharmacophore_to_json(consensus, f'{cons_path}/{cluster}_concensus.json')

    # --- Normalize weights ---
    max_amount_cluster = len(p4_table['ligand'].unique())
    consensus['frequency'] = consensus['weight'] / max_amount_cluster
    consensus.to_csv(f'{results_dir}/important features/consensus_weights_{cluster}.csv', index=False)

    # --- PyMOL visual for each descriptor ---
    subset_names = ['Aromatic', 'Hydrophobic', 'HydrogenAcceptor', 'HydrogenDonor']

    for name in subset_names:
        subset = consensus[consensus['name'] == name]
        cmd.reinitialize()

        for idx, row in subset.iterrows():
            atom_name = f"{row['cluster']}_{idx}"
            cmd.pseudoatom(object=atom_name, pos=[row['x'], row['y'], row['z']], vdw=row['radius'], b=row['frequency'], color=row['color'])
            cmd.label(atom_name, f'"{row["cluster"]}"')

        cmd.spectrum("b", palette=f"white {row['color']}", selection="*")
        cmd.group(name, "*")
        cmd.center("all")
        cmd.show("spheres")
        cmd.save(f"{cons_path}/{name}_clusters_by_weight.pse")

    # --- Submodel creation: top features ---
    selection_rules = {
        'Aromatic': 1,
        'HydrogenDonor': 2,
        'HydrogenAcceptor': 2
    }

    valid_pairs = set()
    for name, top_n in selection_rules.items():
        subset = consensus[consensus['name'] == name]
        top_rows = subset.nlargest(top_n, 'frequency')
        for _, row in top_rows.iterrows():
            valid_pairs.add((row['name'], row['cluster']))

    submodel_df = consensus[[ (row['name'], row['cluster']) in valid_pairs for _, row in consensus.iterrows() ]]
    save_pharmacophore_to_json(submodel_df, out_file=f'{cons_path}/main_model_{cluster}.json')
    save_pharmacophore_to_pymol(submodel_df, out_file=f'{cons_path}/main_model_{cluster}.pse')

    print(f"Done with {cluster}\n")
    plt.close('all')
