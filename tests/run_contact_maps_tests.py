import os
import pandas as pd
from Bio import SeqIO
import CoDIAC
from CoDIAC import contactMap as cm
from CoDIAC import PDBHelper
from CoDIAC import pTyrLigand_helpers as pTyr_helpers

def main():
       
    
    # === CONFIGURATION ===
    PATH = 'Data/Adjacency_files/'
    ann_file = 'temp_pdb_filtered.csv'
    PDB_ID = '5EHP'
    entity = 1
    main = pd.read_csv(ann_file)
    
    # === DOMAINS OF INTEREST ===
    ROI_1 = [110, 216]  # SH2
    ROI_2 = [246, 523]  # PTP_cat
    fastaHeader = 'SH2|PTP_cat'
    
    # === INITIAL OFFSETS ===
    N_offset_ROI_1 = C_offset_ROI_1 = 0
    N_offset_ROI_2 = C_offset_ROI_2 = 0
    
    # === STEP 1: Load Annotations ===
    assert os.path.isfile(ann_file), f"PDB Reference file not found: {ann_file}"
    ann = pd.read_csv(ann_file)
    
    # === STEP 2: Initialize Entity & Domain ===
    entities = PDBHelper.PDBEntitiesClass(ann, PDB_ID)
    assert entity in entities.pdb_dict, f"Entity {entity} not found"
    pdbClass = entities.pdb_dict[entity]
    
    # === STEP 3: Construct Chain Map ===
    chain_A = cm.chainMap(PDB_ID, entity)
    chain_A.construct(PATH)
    chain_A_aligned = cm.translate_chainMap_to_RefSeq(chain_A, pdbClass)
    
    # === STEP 4: Alignment Assertions ===
    assert all(i in chain_A_aligned.unmodeled_list for i in [141, 155, 315, 321]), \
    f"Unexpected unmodeled_list: {chain_A_aligned.unmodeled_list}"
    assert chain_A_aligned.ERROR_CODE == 0, \
        f"Unexpected ERROR_CODE: {chain_A_aligned.ERROR_CODE}"
    
    if hasattr(chain_A_aligned, 'refseq'):
        value = True
        assert isinstance(chain_A_aligned.refseq, str), "refseq should be a string"
    else:
        value = False
        assert isinstance(chain_A_aligned.structSeq, str), "structSeq should be a string"
    
    # === STEP 5: Domain Assertions ===
    # If there are multiple domains, we expect to be analyzing both
    num_of_doms = len(pdbClass.domains)
    if num_of_doms > 1:
        dom_keys = list(pdbClass.domains[0].keys())
        assert 'SH2' in dom_keys, f"'SH2' not found in domain keys: {dom_keys}"
        
    assert ROI_1 == pdbClass.domains[1]['SH2'][:2], (
        f"ROI_1 mismatch: expected {ROI_1}, found {pdbClass.domains[0]['SH2'][:2]}"
    )
    
    # === STEP 6: Write FASTA + Feature Files (No Offset) ===
    filename = 'test_interdomain'
    chain_A_aligned.print_fasta_feature_files(
        ROI_1[0], N_offset_ROI_1,
        ROI_1[1], C_offset_ROI_1,
        ROI_2[0], N_offset_ROI_2,
        ROI_2[1], C_offset_ROI_2,
        'Q06124|PTPN11|2|110|216',
        'PTP_cat',
        filename,
        append=False,
        use_ref_seq_aligned=value
    )
    
    # === STEP 7: Check File Creation ===
    expected_files = [filename + ".fasta", filename + ".fea"]
    for f in expected_files:
        assert os.path.isfile(f), f"Expected file missing: {f}"
    
    # === STEP 8: Write FASTA with C-terminal Offset ===
    filename_offset = 'temp_interdomain_offset'
    C_offset_ROI_1 = 3  # Change offset
    
    chain_A_aligned.print_fasta_feature_files(
        ROI_1[0], N_offset_ROI_1,
        ROI_1[1], C_offset_ROI_1,
        ROI_2[0], N_offset_ROI_2,
        ROI_2[1], C_offset_ROI_2,
        'Q06124|PTPN11|2|110|216',
        'PTP_cat',
        filename_offset,
        append=False,
        use_ref_seq_aligned=value
    )
    
    seq_default = load_fasta_sequence(filename + '.fasta')
    seq_offset = load_fasta_sequence(filename_offset + '.fasta')

    seq_id = 'Q06124|PTPN11|2|110|216'
    len_diff = len(seq_offset[seq_id]) - len(seq_default[seq_id])
    assert len_diff == C_offset_ROI_1, \
        f"Expected sequence length difference of {C_offset_ROI_1}, got {len_diff}"
    
    # === STEP 10: Compare Feature File Line Counts ===
    with open(filename_offset + '.fea') as f1, open(filename + '.fea') as f2:
        lines_offset = [line for line in f1 if line.strip()]
        lines_default = [line for line in f2 if line.strip()]
    
    assert len(lines_offset) > len(lines_default), \
        f"Expected more lines in offset feature file, got {len(lines_offset)} vs {len(lines_default)}"

    print("***Successfully tested Domain-Domain Contact Mapping (within entity extraction)***")
    
    PDB_ID = '6R5G'
    FASTA_HEADER = 'Q06124|PTPN11|2|110|216|lig_5'
    PTM = 'PTR'
    PTM_FILE = 'temp_ptm'
    SH2_FILE = 'temp_dom'
    LIG_ENTITY = 2
    SH2_ENTITY = 1
    N_offset_ROI_1 = C_offset_ROI_1 = 0
    N_offset_ROI_2 = C_offset_ROI_2 = 0
    
    # === Build Contact Maps for Ligand and SH2 Domain ===
    dict_of_lig = CoDIAC.contactMap.return_single_chain_dict(main, PDB_ID, PATH, LIG_ENTITY)
    dict_of_SH2 = CoDIAC.contactMap.return_single_chain_dict(main, PDB_ID, PATH, SH2_ENTITY)
    
    cm_aligned = dict_of_lig['cm_aligned']
    aligned_str = cm_aligned.structSeq
    aligned_ptm_dict = cm_aligned.transDict
    pdb_ptm_dict = dict_of_lig['pdb_class'].transDict
    
    # === Assertions on PTMs and Sequence Matching ===
    assert PTM in aligned_ptm_dict.values() and PTM in pdb_ptm_dict.values(), \
        f"Expected PTM '{PTM}' not found in PTM dictionaries"
    
    assert aligned_str in dict_of_lig['pdb_class'].ref_seq_mutated, \
        "Sequence mismatch between PDB ref file and aligned sequence"
    
    # === PTM Site Analysis ===
    for res in cm_aligned.transDict:
        if res in cm_aligned.resNums and PTM in cm_aligned.transDict[res]:
    
            # First region of interest (7 N-term, 5 C-term)
            n_term1, c_term1 = 7, 5
            res_start, res_end, aligned_str1, tick_labels = pTyr_helpers.return_pos_of_interest(
                cm_aligned.resNums, cm_aligned.structSeq, res, n_term_num=n_term1, c_term_num=c_term1, PTR_value='y')
            expected_len1 = n_term1 + 1 + c_term1
    
            # Second region of interest (2 N-term, 4 C-term)
            n_term2, c_term2 = 2, 4
            res_start, res_end, aligned_str2, tick_labels = pTyr_helpers.return_pos_of_interest(
                cm_aligned.resNums, cm_aligned.structSeq, res, n_term_num=n_term2, c_term_num=c_term2, PTR_value='y')
            expected_len2 = n_term2 + 1 + c_term2
    
            # Assertions on sequence generation
            assert len(aligned_str1) == expected_len1, "Length mismatch in first aligned sequence"
            assert len(aligned_str2) == expected_len2, "Length mismatch in second aligned sequence"
            assert aligned_str1 != aligned_str2, "Offsets did not produce different sequences"
    
    # === Generate Inter-Chain Adjacency Lists ===
    adjList, arr = CoDIAC.contactMap.return_interChain_adj(PATH, dict_of_lig, dict_of_SH2)
    adjList_alt, arr_alt = CoDIAC.contactMap.return_interChain_adj(PATH, dict_of_SH2, dict_of_lig)
    
    # === Validate Residues of Interest Are in Sequence ===
    assert all(i in cm_aligned.resNums for i in adjList.keys()), \
        "Region of interest is outside aligned structure sequence"
    
    # === Extract SH2 Domain Region of Interest ===
    domains = dict_of_SH2['pdb_class'].domains

    SH2_start = int(domains[0]['SH2'][0])
    SH2_stop = int(domains[0]['SH2'][1])

    arr_sub, list_aa_from_sub, list_to_aa_sub = CoDIAC.contactMap.return_arr_subset_by_ROI(arr,
                                                             res_start, res_end, dict_of_lig['cm_aligned'].return_min_residue(),
                                                             SH2_start, SH2_stop, dict_of_SH2['cm_aligned'].return_min_residue())
    
    # === Print FASTA Feature Files ===

    CoDIAC.contactMap.print_fasta_feature_files(arr_alt, dict_of_SH2['cm_aligned'].refseq,
                                                        SH2_start,N_offset_ROI_1, SH2_stop, C_offset_ROI_1,
                                                        dict_of_SH2['cm_aligned'].return_min_residue(),
                                                        res_start,C_offset_ROI_1, res_end, C_offset_ROI_2,
                                                        dict_of_lig['cm_aligned'].return_min_residue(),
                                                        FASTA_HEADER,'pTyr', PTM_FILE, threshold=1, append=True)
    CoDIAC.contactMap.print_fasta_feature_files(arr, dict_of_lig['cm_aligned'].structSeq,
                                                             res_start, N_offset_ROI_1,res_end, C_offset_ROI_1,
                                                            dict_of_lig['cm_aligned'].return_min_residue(),
                                                             SH2_start,N_offset_ROI_2, SH2_stop, C_offset_ROI_2,
                                                            dict_of_SH2['cm_aligned'].return_min_residue(),
                                                            FASTA_HEADER,'SH2', SH2_FILE, threshold=1, append=True )

    
    # === Validate FASTA Output ===
    records = SeqIO.parse(SH2_FILE + '.fasta', 'fasta')
    for fasta in records:
        name, sequence = fasta.id, str(fasta.seq)
        assert sequence == aligned_str2.upper(), \
            f"Mismatch between FASTA output and expected sequence: {sequence} vs. {aligned_str2.upper()}"

    print("***Successfully tested Domain-Ligand Contact Mapping (across entities extraction)***")

# === STEP 9: Compare FASTA Sequence Lengths ===
def load_fasta_sequence(file_path):
    records = SeqIO.parse(file_path, 'fasta')
    return {record.id: str(record.seq) for record in records}

if __name__ == "__main__":
    main()