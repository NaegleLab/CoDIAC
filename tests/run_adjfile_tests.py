import os
import sys
import pandas as pd
from CoDIAC import PDB, AdjacencyFiles
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def main():
    PDB_list = ['4JGH','3TL0','5EHP','6R5G']
    adj_path = './Data/Adjacency_files'
    filepath = './Data/'
    pdb_ref_csv='temp_pdb_structure.csv'

    current_dir = os.getcwd()
    if current_dir not in sys.path:
        sys.path.insert(0, current_dir)

    PDB.download_cifFile(PDB_list, filepath)

    for PDB_ID in PDB_list:
        cif_file = os.path.join(filepath, f'{PDB_ID}.cif')
        json_file = os.path.join(filepath, f'{PDB_ID}.json')
        newpath = os.path.join(adj_path, PDB_ID)

        os.makedirs(newpath, exist_ok=True)

        AdjacencyFiles.ProcessArpeggio(json_file, newpath, cif_file, small_molecule=False)
        AdjacencyFiles.BinaryFeatures(PDB_ID, newpath, translate_resid=True)

        subdir = os.path.join(adj_path, PDB_ID)
        binary_file = os.path.join(subdir, PDB_ID + '_BF.txt')
        mmcif_file = os.path.join(filepath, f'{PDB_ID}.cif')
        
        # Assert that the directory exists
        assert os.path.isdir(subdir), f"[ASSERTION FAILED] Directory not found: {subdir}"
    
        # Assert that the binary BF file exists
        assert os.path.isfile(binary_file), f"[ASSERTION FAILED] {pdb_id}_BF.txt not found in: {subdir}"

        entity_ids_cif = extract_entity_ids_from_cif(mmcif_file)
        # print(f"Entity IDs from CIF: {entity_ids_cif}")
        entity_ids_bf = extract_entities_from_binaryfile(binary_file)
        # print(f"Entity IDs from binary_bf.txt: {entity_ids_bf}")
        missing_in_bf = set(entity_ids_cif) - set(entity_ids_bf)
        assert not missing_in_bf, f"[ASSERTION FAILED] For {PDB_ID} Entity IDs in CIF but missing in binary file: {missing_in_bf}"

        pdb_df = pd.read_csv(pdb_ref_csv)
        if PDB_ID in pdb_df['PDB_ID'].tolist():
            entities_from_pdb = extract_entities_from_PDBref(pdb_ref_csv, PDB_ID)
    
            pdb_and_bf = set(entities_from_pdb) - set(entity_ids_bf)
            assert not pdb_and_bf, f"[ASSERTION FAILED] For {PDB_ID} Entity IDs in PDB_REF but missing in binary file: {pdb_and_bf}"
    
            bf_and_pdb = set(entity_ids_bf) - set(entities_from_pdb)
            assert not bf_and_pdb, f"[ASSERTION FAILED] For {PDB_ID} Entity IDs in binary file but not in PDB_REF: {bf_and_pdb}"

        

def extract_entities_from_binaryfile(binary_bf_path):
    df = pd.read_csv(binary_bf_path, sep='\t')
    unique_values = pd.unique(df[['Entity1', 'Entity2']].values.ravel())
    return list(unique_values)
    
def extract_entity_ids_from_cif(cif_path):
    polymer_entity_ids = []
    pdb_info = MMCIF2Dict(cif_path)
    for i in range(len(pdb_info["_entity.id"])):  # Iterate through the entity IDs
        if pdb_info["_entity.type"][i] == "polymer": # Check entity type
            polymer_entity_ids.append(pdb_info["_entity.id"][i]) # Add the ID if it's a non-polymer
    return [int(x) for x in polymer_entity_ids]

def extract_entities_from_PDBref(pdb_ref_csv, pdb_id):
    pdb_df = pd.read_csv(pdb_ref_csv)
    return pdb_df.loc[pdb_df['PDB_ID'] == pdb_id, 'ENTITY_ID'].tolist() 
    
if __name__ == "__main__":
    main()

