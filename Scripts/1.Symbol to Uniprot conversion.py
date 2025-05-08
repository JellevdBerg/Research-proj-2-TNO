import pandas as pd
import mygene

# Load Excel file (replace 'input.xlsx' with your actual file)
file_path = 'Data/1.RAW.xlsx'  # Update with the actual file path
df = pd.read_excel(file_path)

gene_symbols = df['Symbol'].dropna().tolist()  # Ensure no NaN values

def convert_symbols_to_uniprot(gene_symbols):
    mg = mygene.MyGeneInfo()
    result = mg.querymany(gene_symbols, scopes='symbol', fields='uniprot', species='human')
    
    symbol_to_uniprot = {}
    for entry in result:
        gene_symbol = entry['query']
        uniprot_id = ''
        if 'uniprot' in entry:
            uniprot_info = entry['uniprot']
            if isinstance(uniprot_info, list):
                uniprot_id = uniprot_info[0].get('Swiss-Prot', '')
            else:
                uniprot_id = uniprot_info.get('Swiss-Prot', '')
        symbol_to_uniprot[gene_symbol] = uniprot_id if uniprot_id else None
    
    return symbol_to_uniprot

# Convert symbols to UniProt IDs
symbol_to_uniprot = convert_symbols_to_uniprot(gene_symbols)

df['Uniprot'] = df['Symbol'].map(symbol_to_uniprot)

# Save output to CSV
df.to_csv('Data/csv data/1.raw_genes.csv', index=False)