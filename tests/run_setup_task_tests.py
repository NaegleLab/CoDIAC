 
from CoDIAC import UniProt, InterPro
import sys
import os
 
# Add the directory containing this script to sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
import test_UniProt
import test_PDB_fetch_annotate_filter

def main():

    domain_test_file = 'temp_uniprot_reference.csv' 
    # run the UniProt reference file creation in test 
    print("Creating UniProt reference file and tests...")
    test_UniProt.main(domain_test_file)

    print("Testing PDB fetch, annotation, filtering...")
    PDB_test_file='temp_pdb_structure.csv'
    annotated_test_file='temp_pdb_annotated.csv'
    Interpro_ID = 'IPR000980'  # SH2 domain
    PDB_file_filtered='temp_pdb_filtered.csv'
    test_PDB_fetch_annotate_filter.main(domain_test_file, PDB_test_file, annotated_test_file, Interpro_ID, PDB_file_filtered)


    print("Testing UniProt domain fasta creation...")
    # let's now make an SH2 fasta file and test that it works correctly, including expanding the domain boundaries
    fasta_test_file = 'temp_sh2.fasta'
    Interpro_ID = 'IPR000980'  # SH2 domain

    N_OFFSET = 0
    C_OFFSET = 0
    sh2_fasta_dict = UniProt.print_domain_fasta_file(domain_test_file, Interpro_ID, fasta_test_file, N_OFFSET, C_OFFSET, APPEND=False)
    assert len(sh2_fasta_dict) == 3, "Expected 3 SH2 domains in the fasta file, got {}".format(len(sh2_fasta_dict))

    print("Testing fasta file creation with different offsets...")
    N_OFFSET = 5
    C_OFFSET = 5
    sh2_fasta_dict_n_term_offset = UniProt.print_domain_fasta_file(domain_test_file, Interpro_ID, fasta_test_file, N_OFFSET, C_OFFSET, APPEND=False)
    assert len(sh2_fasta_dict_n_term_offset) == 3, "Expected 3 SH2 domains in the fasta file, got {}".format(len(sh2_fasta_dict_n_term_offset))
    # check that the sequences are longer than the original ones
    for uniprot_id in sh2_fasta_dict_n_term_offset:
        assert len(sh2_fasta_dict_n_term_offset[uniprot_id]) == len(sh2_fasta_dict[uniprot_id])+10, f"Expected sequence for {uniprot_id} to be longer with offsets, but it was not."
    

if __name__ == "__main__":
    main()