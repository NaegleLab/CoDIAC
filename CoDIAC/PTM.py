
from CoDIAC import IntegrateStructure_Reference, UniProt, proteomeScoutAPI, PhosphoSitePlus_Tools
import pandas as pd
import os 

package_directory = os.path.dirname(os.path.abspath(__file__))
PROTEOMESCOUT_DATA = package_directory + '/data/proteomescout_everything_20190701/data.tsv'
PTM_API = proteomeScoutAPI.ProteomeScoutAPI(PROTEOMESCOUT_DATA)


def write_PTM_features(Interpro_ID, uniprot_ref_file, feature_dir, mapping_file ='', n_term_offset=0, c_term_offset=0, gap_threshold=0.7, num_PTM_threshold = 5, PHOSPHOSITE_PLUS=False):
    """
    Writes all PTM features from ProteomeScout or PhosphoSitePlus on Interpro domains from a uniprot reference file, if there are more
    than num_PTM_threshold that occur across all domains of that type in the reference. Returns the ptm_count_dict for reference
    and the feature dict that is generated to write the files. Files are named Interpro_ID_<PTM_Type>.feature
    and the reference fasta is also generated so that it is clear the features are attached to that 
    particular run of the domains. 

    Parameters
    ----------
    Interpro_ID: string 
        Interpro ID - for example in a reference line such as 
        SH3_domain:IPR001452:82:143; SH2:IPR000980:147:246; Prot_kinase_dom:IPR000719:271:524
        the interpro ID for the SH3_domain is IPR001452; for the SH2 domain is IPR000980
    uniprot_reference_file: string
        File location that contains the reference of interest (like produced from Uniprot.makeRefFile)
    feature_dir: string
        Feature Directory to place files in
    mapping_file: string
        A CSV file location, if wanted, that holds a translation of the long header into a shorter header
        If this is an empty string, then it will not attempt mapping
    n_term_offset: int
        Number of amino acids to extend in the n-term direction (up to start of protein)
    c_term_offset: int
        Number of amino acids to extend in the c-term direction (up to end of protein)
    gap_threshold: float
        fraction gap allowed before dispanding with PTM translation from ProteomeScout
    num_PTM_threshold: int
        Number of PTMs in all domains of a type required to generate a feature file
    PHOSPHOSITE_PLUS: bool
        If True, will generate PTMs from PhosphoSitePlus instead of ProteomeScout. 
        See PhosphoSitePlus_Tools.py convert_pSiteDataFiles can be used to update or create the 
        API-formatted files. These resources are stored in GitHub LFS.
    

    Returns
    -------
    file_list: list
        List of files generated as features
    ptm_count_dict: dict
        keys are the modification type and values are total number encountered
    domain_feature_dict: dict of dicts
        keys are the fasta headers 
        Inner dict keys are the modification types and the values of this is a list of zero-based positions. 

    
    """
    ptm_feature_file_list = []
    if not os.path.isdir(feature_dir):
        print("Directory for target files does not exist, please create it.")
        return ptm_feature_file_list

    #print("DEBUG: PhosphoSitePlus is %s"%(PHOSPHOSITE_PLUS))
    #write the fasta into the directory - this actually doesn't make sense, it's the wrong headers if mappingDict is incorrect
    #fasta_output_file = feature_dir + Interpro_ID + '.fasta'
    domain_dict = UniProt.make_domain_fasta_dict(uniprot_ref_file, Interpro_ID, n_term_offset, c_term_offset)

    # get the mapping dictionary
    mapping_dict = {}
    if mapping_file == '':
        for header in domain_dict:
            mapping_dict[header] = header
    else:
        #test if the mapping file looks correct
        mapping_df = pd.read_csv(mapping_file)
        if 'short' in mapping_df.columns and 'full' in mapping_df.columns:
            for index, row in mapping_df.iterrows():
                long_header = '>'+row['full']
                short = '>'+row['short']
                mapping_dict[long_header] = short
        else:
            print("ERROR: Mapping file is not correctly shapped, skipping mapping")
            for header in domain_dict:
                mapping_dict[header] = header


    color = '3546b5'
    ptm_count_dict, ptm_feature_dict = get_Interpro_PTMs(Interpro_ID, uniprot_ref_file, PHOSPHOSITE_PLUS, n_term_offset, c_term_offset, gap_threshold)

    if num_PTM_threshold <=0:
        num_PTM_threshold = 1 
        print("ERROR: Resetting the PTM number threshold to 1, this value is an int and should be greater than 0")

    ptm_list = []
    for ptm_type in ptm_count_dict:
        if ptm_count_dict[ptm_type] >= num_PTM_threshold:
            ptm_list.append(ptm_type)

    # for each feature type, create the feature file, printing features for each protein that has them. 
    for ptm_type in ptm_list:
        ptm_dict_temp = {}
        ptm_feature_file = feature_dir + Interpro_ID + '_' + ptm_type + '.feature'
        for header in ptm_feature_dict:
            if ptm_type in ptm_feature_dict[header]:
                if header in mapping_dict:
                    short_header = mapping_dict[header] #if a mapping file was passed then we'll use it, unless that header had no map
                else:
                    print("DID not find short header for %s"%(header))
                    short_header = header
                ptm_dict_temp[short_header] = ptm_feature_dict[header][ptm_type]
        write_PTM_feature_file(ptm_type, color, ptm_dict_temp, ptm_feature_file)
        ptm_feature_file_list.append(ptm_feature_file)

    return ptm_feature_file_list, ptm_count_dict, ptm_feature_dict, mapping_dict

def translate_PTMs(uniprot_ID, uniprot_seq, gap_threshold, PHOSPHOSITE_PLUS):
    """
    Given an Uniprot ID and a uniprot sequence, return PTMs on the uniprot sequence positions. Uses the pairwise alignment
    of proteomescout sequence to the uniprot sequence to build a position map and translates PTMs into 
    the uniprot reference positions.

    Parameters
    ----------
    uniprot_ID: str
        Uniprot ID
    uniprot_seq: str
        Uniprot sequence
    gap_threshold: float
        Allowed gap fraction before considering too poor an alignment
    PHOSPHOSITE_PLUS: bool
        Will use PhosphoSitePlus data if True, ProteomeScout if False
    
    Returns
    -------
    errors: string
        A string error description if there was an issue with alignment or 
    translated_PTMs: tuples of PTMs
        tuples described by (pos, aa, mod_type)
        position is in ones-based counting in the uniprot sequence
    failed_PTMs: tuples
        Keeps track of PTMs that did not translate and the reason
    
    """
    #print("DEBUG: PhosphoSitePlus is %s"%(PHOSPHOSITE_PLUS))
    if not PHOSPHOSITE_PLUS:
        proteomescout_seq, PTMs = get_PTMS_proteomeScout(uniprot_ID)
    else:
        proteomescout_seq, PTMs = get_PTMS_phosphoSitePlus(uniprot_ID)
        #print("DEBUG: Working on %s"%(uniprot_ID))
        #print(proteomescout_seq)
        #print(PTMs)
    #proteomescout_seq = PTM_API.get_sequence(uniprot_ID)
    #PTMs = PTM_API.get_PTMs(uniprot_ID)
    errors = ""
    if proteomescout_seq == '-1':
        errors = "Error: PTM record not found by %s"%(uniprot_ID)
        print(errors)
        return errors, [], []
    aln, struct_sequence_ref_spanning, from_start, from_end, range, pos_diff, diffList, gaps_ref_to_struct, gaps_struct_to_ref = IntegrateStructure_Reference.return_mapping_between_sequences(proteomescout_seq, uniprot_seq, 1, 1, len(uniprot_seq))
    
    #mapToRef = aln.get_gapped_seq('reference').gap_maps()[1]
    map_to_ref = aln.get_gapped_seq('structure').gap_maps()[0]
    numGaps = aln.seqs[0].count_gaps()

    #check that the alignment is of sufficient quality
    if (numGaps/len(uniprot_seq) >= gap_threshold):
        errors = "Error: Percent gap too large %s, has %0.2f percent gap"%(uniprot_ID, numGaps/len(uniprot_seq))
        return errors, [], []

    structure_aln = aln.get_degapped_relative_to('structure')
    aa_list  = list(structure_aln)
    translated_PTMs = []
    failed_PTMs = []
    proteomescout_ind = 0
    uniprot_ind = 1
    for a_PTM in PTMs:
        #print(a_PTM)
        pos, aa, ptm_type = a_PTM
        pos = int(pos) #this is ones-based counted
        if aa_list[pos-1][proteomescout_ind]!=aa:
            #print("error: proteomescout position is %s not %s"%(aa_list[int(pos)-1][proteomescout_ind], aa))
            failed_PTMs.append((pos, aa, ptm_type, 'PTM reference issue, non-matching amino acid'))
        else:
            if aa_list[pos-1][uniprot_ind]!=aa:
                #print("Skipping %s, proteomescout and uniprot don't match"%(pos))
                failed_PTMs.append((pos, aa, ptm_type, 'different amino acid'))
            else:
                pos_translated = map_to_ref[pos-1]+1
                translated_PTMs.append((str(pos_translated), aa, ptm_type))

    return(errors, translated_PTMs, failed_PTMs)


def get_PTMS_proteomeScout(uniprot_ID):
    """
    Parameters
    ----------
    uniprot_ID: str
        Uniprot ID
    
    Returns
    -------
    seq: str
        Sequence of the protein. 
    PTMs: list of tuples
        List of tuples with PTM information in [('position', 'amino acid', 'modification')]
    
    """

    proteomescout_seq = PTM_API.get_sequence(uniprot_ID)
    PTMs = PTM_API.get_PTMs(uniprot_ID)

    return proteomescout_seq, PTMs


def get_PTMS_phosphoSitePlus(uniprot_ID):
    """
    Found cases where PhosphoSitePlus had no modification annotations under the canonical number (e.g. Q9HBL0, TNS1), but had them under a 
    isoform number Q9HBL0-1. So added the ability to check isoforms, until no isoforms are found. This now returns a sequence and PTMs for the record found. 
    The isoform differences should be handled by alignment to reference. If the isoform with PTMs (starting with 1 and going forward) is very different, 
    failure will occur in mapping and if slightly different, numbering will be corrected. 

    Parameters
    ----------
    uniprot_ID: str
        Uniprot ID
    
    Returns
    -------
    seq: str
        Sequence of the protein. Please note that we may have moved to an isoform in the search for PTMs in PSitePlus.
        Sequence returned will reflect the sequence PTMs are attached to. 
    PTMs: list of tuples
        List of tuples with PTM information. [('position', 'amino acid', 'modification')]
    
    """

    PTMs, ID_used = PhosphoSitePlus_Tools.get_PTMs(uniprot_ID)

    seq = PhosphoSitePlus_Tools.get_sequence(ID_used) #make sure if the ID changed, we get the right sequence
    return seq, PTMs

def get_Interpro_PTMs(Interpro_ID, uniprot_reference_file, PHOSPHOSITE_PLUS, n_term_offset=0, c_term_offset=0, gap_threshold=0.7):
    """
    Given an uniprot reference file and a particular InterproID of interest, get all the PTMs that exist
    within the domains in the uniprot reference file that have that interpro ID. 
    Creates a dictionary count of all PTMs encoungered across all domains of that type in the file

    Parameters
    ----------
    Interpro_ID: string 
        Interpro ID - for example in a reference line such as 
        SH3_domain:IPR001452:82:143;SH2:IPR000980:147:246;Prot_kinase_dom:IPR000719:271:524
        the interpro ID for the SH3_domain is IPR001452; for the SH2 domain is IPR000980
    uniprot_reference_file: string
        File location that contains the reference of interest (like produced from Uniprot.makeRefFile)
    n_term_offset: int
        Number of amino acids to extend in the n-term direction (up to start of protein)
    c_term_offset: int
        Number of amino acids to extend in the c-term direction (up to end of protein)
    gap_threshold: float
        fraction gap allowed before dispanding with 
    PHOSPHOSITE_PLUS: bool
        If True, will generate PTMs from PhosphoSitePlus instead of ProteomeScout. 
        This requires running your own local setup of PhosphoSite files. See PhosphoSitePlus_Tools.py convert_pSiteDataFiles
    
    Returns
    -------
    ptm_count_dict: dict
        keys are the modification type and values are total number encountered
    domain_feature_dict: dict of dicts
        keys are the fasta headers 
        Inner dict keys are the modification types and the values of this is a list of zero-based positions. 


    """
    uniprot_df = pd.read_csv(uniprot_reference_file)
    domain_fasta_dict = UniProt.make_domain_fasta_dict(uniprot_reference_file, Interpro_ID, n_term_offset, c_term_offset)

    if gap_threshold < 0 or gap_threshold > 1:
        gap_threshold = 0.7
        print("Gap threshold a number between 0 and 1, setting to default of 0.7")

    ptm_domain_dict = {}
    translated_PTMs_dict = {}
    failed_PTMs_dict = {}
    errors_dict = {}
    for index, row in uniprot_df.iterrows():
        seq = row['Ref Sequence']
        uniprot_id = row['UniProt ID']
        domains = row['Interpro Domains']
        errors_dict[uniprot_id], translated_PTMs_dict[uniprot_id], failed_PTMs_dict[uniprot_id] = translate_PTMs(uniprot_id, seq, gap_threshold, PHOSPHOSITE_PLUS)

        #Writing/printing a report of global issues encountered.
        if errors_dict[uniprot_id]:
            print(errors_dict[uniprot_id])
        elif failed_PTMs_dict[uniprot_id]:
            print(failed_PTMs_dict[uniprot_id])

    # Next, walk through each fasta sequence, and find the translated IDs, if they exist in the domains of interest. 
    domain_feature_dict = {}
    for header in domain_fasta_dict:
        domain_feature_dict[header] = {}
        header_temp = header.strip(">")
        u_id_from_header, gene_name, species, domain_name, domainNum, IP_id_from_header, start_from_header, end_from_header = header_temp.split("|")
        translated_PTMs = translated_PTMs_dict[u_id_from_header]
        if translated_PTMs:
            for a_PTM in translated_PTMs:
                pos, aa, ptm_type = a_PTM
                pos = int(pos)
                if pos >= int(start_from_header) and pos < int(end_from_header):
                        if ptm_type not in ptm_domain_dict:
                            ptm_domain_dict[ptm_type] = 0
                        ptm_domain_dict[ptm_type] +=1
                        dom_pos = pos-int(start_from_header) #this is now a zero-based counting
                        if ptm_type not in domain_feature_dict[header]:
                            domain_feature_dict[header][ptm_type] = []
                        # try:
                        #     if domain_fasta_dict[header][dom_pos] != aa:
                        #         print("ERROR: trying to place a feature in domain for %s at %d for domain number %d"%(u_id_from_header, dom_pos, domainNum))
                        #         #uni_seq = uniprot_df[uniprot_df['UniProt ID']==u_id_from_header]['Ref Sequence']
                        #         #print("aa should be %s and in the full protein it is %s"%(aa, uni_seq[pos-1]))
                        # except:
                        #     print("Exception: feature is %d and length of domain is %d for %s"%(dom_pos, len(domain_fasta_dict[header]), header))

                        domain_feature_dict[header][ptm_type].append(dom_pos+1) #one's based position for the features
    

    # elif translated_PTMs:
    #     for domain in domains.split(';'):
    #         domain_name, IP_id, start, stop = domain.split(':')
    #         if IP_id == Interpro_ID:
    #             for a_PTM in translated_PTMs:
    #                 pos, aa, ptm_type = a_PTM
    #                 if pos >= start and pos <= stop:
    #                     if ptm_type not in ptm_domain_dict:
    #                         ptm_domain_dict[ptm_type] = 0
    #                     ptm_domain_dict[ptm_type] +=1
    return ptm_domain_dict, domain_feature_dict
            
def write_PTM_feature_file(PTM_name, PTM_color, feature_dict, output_file):
    """
    Given a PTM that should be named 'PTM_name' and the color it should be assigned
    write a Jalview compatible feature file that sets that feature name and color and then
    prints for every header in feature_dict the positions of that feature type to output_file

    Parameters
    ----------
    PTM_name: string
        Feature name
    PTM_color: string
        Color string (compatible with Jalview conventions) for the color to assigne for a feature
    feature_dict: dict
        Has keys equal to the headers to be used and values is a list of positions to print as features
    output_file: string
        file that will be written.. overwrites an existing file.

    """

    with open(output_file, 'w') as file:
        file.write('{}\t{}\n'.format(PTM_name, PTM_color))
        for header in feature_dict:
            for feature_pos in feature_dict[header]:
                file.write('{}\t{}\t-1\t{}\t{}\t{}\n'.format(PTM_name, header.strip(">"), feature_pos, feature_pos, PTM_name))
    file.close()

