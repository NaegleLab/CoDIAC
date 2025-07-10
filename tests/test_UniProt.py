from CoDIAC import UniProt



def main(domain_test_file='temp_uniprot_reference.csv'):


    uniprot_ids = ['Q06124', 'Q7Z7G1', 'Q9HB58-6'] # Q06124-PTPN11 expecting two domains of one type Q7Z7G1 is CLNK, just one domain; Q9HB58-6 an isoform wtih bromodomain
    UniProt.makeRefFile(uniprot_ids, domain_test_file)

    # let's check the file
    with open(domain_test_file, 'r') as f:
        lines = f.readlines()
        assert len(lines) == 4, "Expected 4 lines in the file, got {}".format(len(lines))
        assert lines[0].startswith('UniProt ID'), "First line should be header"
        assert len(lines[0].split(',')) == 9, "Expected 9 columns in the header, got {}".format(len(lines[0].split(',')))
        assert 'PTPN11' in lines[1], "Expected PTPN11 domain information in the second line"
        assert 'SH2|SH2|PTP_cat' in lines[1], "Expected PTPN11 domain types in the second line"
        assert 'CLNK' in lines[2], "Expected CLNK domain information in the third line"
        assert 'Bromodomain' in lines[3], "Expected Bromodomain domain information in the fourth line"



if __name__ == "__main__":
    main()