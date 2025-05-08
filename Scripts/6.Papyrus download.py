from qsprpred.data.sources.papyrus import Papyrus
import pandas as pd
import numpy as np

# path to CSV data
csv_path = "./Data/csv data/1.raw_genes.csv"

uniprot_col = "Uniprot" # read the CSV Uniprot column
acc_keys = pd.read_csv(csv_path)[uniprot_col].dropna().unique().tolist() # for now limit to top 10 rows
dataset_name = "Papyrus_dataset"  # name of the file to be generated
quality = "high"  # choose minimum quality from {"high", "medium", "low"}
papyrus_version = "05.7"  # Papyrus database version
data_dir = "./Data/Papyrus"  # directory to store the Papyrus data
output_dir = "./Data/Papyrus_output"  # directory to store the generated dataset

# Create a Papyrus object, which specifies the version and directory to store the payrus data
papyrus = Papyrus(
    data_dir=data_dir,
    version=papyrus_version,
    stereo=False,
    plus_only=False,
)

# Create subset of payrus data for the given accession keys, returns a MoleculeTable
mt = papyrus.getData(
    dataset_name,
    acc_keys,
    quality,
    output_dir=output_dir,
    use_existing=True,
    activity_types=["Ki", "IC50", "Kd"]
)
mt.getDF()