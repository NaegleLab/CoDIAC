from CoDIAC import PDB, IntegrateStructure_Reference
import pandas as pd



def main(domain_test_file='temp_uniprot_reference.csv', PDB_test_file='temp_pdb_structure.csv', annotated_test_file='temp_pdb_annotated.csv', Interpro_ID='IPR000980', PDB_file_filtered='temp_pdb_filtered.csv'):

    #domain_test_file = 'temp_uniprot_reference.csv'
    PDB.generateStructureRefFile_fromUniprotFile(domain_test_file, PDB_test_file)

    # the total number of PDB IDs is found in the PDB Ids colulmn of the domain test file
    domain_df = pd.read_csv(domain_test_file)
    pdb_list = []

    # It's hard to code exactly the number of lines expected, since the PDB may grow and because some structures will
    # have multiple lines associated with them (multiple chais). Instead, we'll read this into a dataframe and 
    # check for a few things that are certain to be known. 
    pdb_df = pd.read_csv(PDB_test_file)

    PDB_IDs_to_check = ['2SHP', '4JE4', '8U7X']
    for pdb_id in PDB_IDs_to_check:
        assert pdb_id in pdb_df['PDB_ID'].values, f"Expected {pdb_id} to be in the PDB IDs, but it was not found."

    # check the columns include expected columns
    expected_columns = ['PDB_ID', 'rcsb_uniprot_protein_sequence', 'CHAIN_ID', 'title']
    for col in expected_columns:
        assert col in pdb_df.columns, f"Expected column {col} to be in the PDB DataFrame, but it was not found."

    struct_df_out = IntegrateStructure_Reference.add_reference_info_to_struct_file(PDB_test_file, domain_test_file, annotated_test_file, INTERPRO=True, verbose=False)
    struct_df_filter =IntegrateStructure_Reference.filter_structure_file(annotated_test_file, Interpro_ID, PDB_file_filtered)
    # 

    # open and test that there are lines, that it's equal or less than the number of lines in unfiltered file
    with open(PDB_file_filtered, 'r') as f:
        lines = f.readlines()
        assert len(lines) > 0, "Expected some lines in the filtered PDB file, but got none."
        assert len(lines) <= len(pd.read_csv(annotated_test_file)), "Filtered PDB file should not have more lines than the annotated file."
        # assert that the first line is a header
        assert lines[0].startswith('PDB_ID'), "First line of the filtered PDB file should be a header."



    # # let's check the file
    # with open(domain_test_file, 'r') as f:
    #     lines = f.readlines()
    #     assert len(lines) == 4, "Expected 4 lines in the file, got {}".format(len(lines))
    #     assert lines[0].startswith('UniProt ID'), "First line should be header"
    #     assert len(lines[0].split(',')) == 9, "Expected 9 columns in the header, got {}".format(len(lines[0].split(',')))
    #     assert 'PTPN11' in lines[1], "Expected PTPN11 domain information in the second line"
    #     assert 'SH2|SH2|PTP_cat' in lines[1], "Expected PTPN11 domain types in the second line"
    #     assert 'CLNK' in lines[2], "Expected CLNK domain information in the third line"

    


if __name__ == "__main__":
    main()