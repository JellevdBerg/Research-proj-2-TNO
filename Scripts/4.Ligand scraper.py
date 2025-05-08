import requests
import pandas as pd
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor
import time  # Import the time module
from tenacity import retry, stop_after_attempt, wait_exponential  # Import tenacity for retry logic

# Function to get ligands from a PDB entry
def get_ligands(pdb_id):
    try:
        pdb_info = _fetch_pdb_nonpolymer_info(pdb_id)
        
        # Check if nonpolymer_entities is None or empty
        if not pdb_info or not pdb_info.get('data', {}).get('entry', {}).get('nonpolymer_entities'):
            print(f"No ligands found for PDB ID {pdb_id}")
            return None  # or you can return an empty dict, depending on your preference
        
        nonpolymer_entities = pdb_info['data']['entry'].get('nonpolymer_entities', [])
        
        ligand_expo_ids = [
            nonpolymer_entities_item["pdbx_entity_nonpoly"]["comp_id"]
            for nonpolymer_entities_item in nonpolymer_entities
        ]
        
        if not ligand_expo_ids:
            return None  # Return None if no ligands are found

        ligands = {}
        with ThreadPoolExecutor(max_workers=11) as executor:
            # Fetch ligand data asynchronously
            future_to_ligand = {executor.submit(_fetch_ligand_expo_info, ligand_expo_id): ligand_expo_id for ligand_expo_id in ligand_expo_ids}
            for future in future_to_ligand:
                ligand_expo_id = future_to_ligand[future]
                ligand_expo_info = future.result()
                ligands[ligand_expo_id] = ligand_expo_info
                time.sleep(0)  # Add a delay between requests

        return ligands
    except Exception as e:
        print(f"Error processing PDB ID {pdb_id}: {e}")
        return None  # Return None if any exception occurs

# Fetch non-polymer data (ligand info) from RCSB using GraphQL
@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=4, max=10))
def _fetch_pdb_nonpolymer_info(pdb_id):
    query = (
        """{
        entry(entry_id: "%s") {
            nonpolymer_entities {
                pdbx_entity_nonpoly {
                    comp_id
                    name
                    rcsb_prd_id
                }
            }
        }
    }"""
    % pdb_id
)

    query_url = f"https://data.rcsb.org/graphql?query={query}"
    response = requests.get(query_url)
    response.raise_for_status()
    info = response.json()
    return info

# Fetch detailed ligand data from ligand-expo website using HTML parsing
@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=4, max=10))
def _fetch_ligand_expo_info(ligand_expo_id):
    r = requests.get(f"http://ligand-expo.rcsb.org/reports/{ligand_expo_id[0]}/{ligand_expo_id}/")
    r.raise_for_status()
    html = BeautifulSoup(r.text, 'html.parser')
    info = {}
    for table in html.find_all("table"):
        for row in table.find_all("tr"):
            cells = row.find_all("td")
            if len(cells) != 2:
                continue
            key, value = cells
            if key.string and key.string.strip():
                info[key.string.strip()] = "".join(value.find_all(string=True))

    # Postprocess some known values (convert to appropriate data types)
    info["Molecular weight"] = float(info["Molecular weight"].split()[0])
    info["Formal charge"] = int(info["Formal charge"])
    info["Atom count"] = int(info["Atom count"])
    info["Chiral atom count"] = int(info["Chiral atom count"])
    info["Bond count"] = int(info["Bond count"])
    info["Aromatic bond count"] = int(info["Aromatic bond count"])
    return info

df = pd.read_csv('Data/csv data/2.raw_genes_long.csv')

# Grab the top 10 values from a specific column
selected_pdb_ids = df['PDB-ID'].tolist()

# Split the PDB IDs into chunks of 4000
chunk_size = 4000
chunks = [selected_pdb_ids[i:i + chunk_size] for i in range(0, len(selected_pdb_ids), chunk_size)]

# Columns for the output DataFrame
columns = [
    "@structureId", 
    "@chemicalID", 
    "@type", 
    "@molecularWeight", 
    "InChIKey", 
    "smiles"
]

# Process each chunk of pdb_ids
for chunk_index, chunk in enumerate(chunks):
    rows = []
    with ThreadPoolExecutor(max_workers=11) as executor:
        futures = {executor.submit(get_ligands, pdb_id): pdb_id for pdb_id in chunk}
        for future in futures:
            pdb_id = futures[future]
            ligands = future.result()
            if ligands is None:  # Skip PDBs that do not have ligands
                continue
            # Filter ligands with molecular weight > 100
            ligands = {ligand_id: properties for ligand_id, properties in ligands.items() if properties["Molecular weight"] > 100}
            if not ligands:  # Skip if no ligands are left after filtering
                continue
            # Select the largest ligand by molecular weight
            ligand_id, properties = max(ligands.items(), key=lambda kv: kv[1]["Molecular weight"])
            
            # Check if any required property is missing
            if not all(key in properties for key in ["Component type", "Molecular weight", "InChIKey descriptor", "Stereo SMILES (OpenEye)"]):
                print(f"Skipping PDB ID {pdb_id} due to missing properties")
                continue  # Skip PDB ID if any property is missing

            # Add the ligand data to the rows if all required properties are present
            rows.append([
                pdb_id,
                ligand_id,
                properties["Component type"],
                properties["Molecular weight"],
                properties["InChIKey descriptor"],
                properties["Stereo SMILES (OpenEye)"],
            ])
            
    
    # Convert the ligand data to a DataFrame for easy inspection
    ligands_df = pd.DataFrame(rows, columns=columns)

    # Save each chunk to a separate CSV file
    chunk_filename = f"Data/PDB Ligands/ligands_from_pdb_chunk_{chunk_index + 1}.csv"
    ligands_df.to_csv(chunk_filename, index=False)
    print(f"Ligand data for chunk {chunk_index + 1} saved to '{chunk_filename}'")
    
# running with head
# ook ff naar selenium kijken
