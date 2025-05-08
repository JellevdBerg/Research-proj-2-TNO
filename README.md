# RP2-TNO
## Overview of project aim
This project is focused on identifying novel therapeutic targets within Sarcopenia by applying a diverse set of advanced data science techniques. After pinpointing these potential targets, we will implement a classical in-silico workflow, which involves multiple stages such as training machine learning models to predict biological activity, verifying the activity of identified compounds, and performing computational docking analyses to explore potential binding interactions. Additionally, molecular dynamics simulations will be carried out to further investigate the stability and behavior of compound-target interactions over time. The ultimate objective of this project is to discover a novel compound targeting an uncharacterized or underexplored target, validate its binding poses through docking studies, and experimentally verify its affinity in in-vitro analyses, ultimately contributing to the development of new therapeutic strategies for Sarcopenia.

### Dependencies
All data is stored within the “data” folder. The scripts are designed to run seamlessly as long as the repository is cloned correctly, with no need for modifications to file paths or directories. The execution order of the scripts and notebooks follows a clear, numbered sequence, ensuring the workflow operates smoothly.

> [!IMPORTANT]
>This repository relies on a very specific set of python packages and versions in order to function properly, the following step by step plan for recreation is advised:
1. **Create a virtual environment** with Python 3.12.2
   ```
   create conda -n spock python=3.12.2
   ```
2. **Install Vina** using:  
   ```
   conda install -c conda-forge vina
   ```
3. **Clone the [Spock](https://github.com/martin-sicho/spock/tree/dev) repository** (developed by M. Sicho):  
   ```
   git clone https://github.com/martin-sicho/spock/tree/dev
   ```
4. **Clone the [QSPRpred](https://github.com/CDDLeiden/QSPRpred/tree/feature/representations_storage) repository** (developed by the LACDR-CDD group):  
   ```
   git clone https://github.com/CDDLeiden/QSPRpred/tree/feature/representations_storage
   ```
5. **Navigate to each repository directory** using `cd "INSERT_PATH_TO_REPO"` and install each package with:  
   ```
   pip install .
   ```  
   Make sure **Spock** is installed before **QSPRpred**.
6. **Install the correct version of RDKit**:  
   ```
   pip install rdkit==2023.9.1
   ```
7. **If the workflow reports an issue with Dimorphite**, install the required version:  
   ```
   conda install -c conda-forge dimorphite-dl=1.3.2
   ```

Following these steps will create a workable conda environment capable of running all scripts in this repository.

### Script overview
The scripts and notebooks are numbered sequentially in the repository, which dictates the order in which they should be run. Each script corresponds to a specific step in the workflow. The first script requires a **"Starting List.xlsx"** file containing a column named **"Symbol"**, which lists the gene symbols of the targets of interest.

1. **Convert Gene Symbols to Uniprot IDs**: This script converts the gene symbols to Uniprot IDs, which are more widely applicable for data science tasks, and saves them to a file.
2. **Scrape PDB-IDs from RCSB PDB.com**: Retrieves PDB-IDs based on the Uniprot IDs and stores them in a comma-separated dataframe.
3. **Pivot PDB Data**: Converts the previously saved PDB data into a long format.
4. **Scrape unique ligands from RCSB PDB.com**: Fetches the unique ligands associated with the obtained PDB-IDs (currently a work in progress).
5. **Download Papyrus Database Subset**: Uses QSPRpred and the Uniprot IDs of the targets to download a subset of the Papyrus database.
6. **Store Papyrus Data**: Organizes the obtained Papyrus data into a more accessible format.
   - **6b**: Filters and visualizes the obtained activity data, combining it with the PDB data for easy target selection.

> [!NOTE]
> This project is still a work in progress, many of the scripts are still to be added, modified, or relocated.