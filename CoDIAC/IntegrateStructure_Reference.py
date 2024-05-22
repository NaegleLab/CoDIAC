from CoDIAC import alignmentTools
import bisect
import pandas as pd


def return_mapping_between_sequences(struct_sequence, ref_sequence, ref_start, pdb_start, ref_length):
    """
    Return information about the pairwise alignment between a structure sequence (mapping from this sequence) 
    to a reference sequence (to_sequence). Specifically, this assumes that the structure sequence is some shortened 
    from of the reference and that you would like to understand where it lands in the reference and if that 
    numbering is similar between structure and reference. 

    Parameters
    ----------
    struct_sequence: str
        structure sequence
    ref_sequence: str
        reference sequence
    ref_start: int
        the one's based counting of where the ref_start is documented as being
    pdb_start: int 
        the one's based counting of where pdb begins matching reference
    ref_length: int 
        the length of the sequence in the experimental structure that matches the reference
    
    Returns
    -------
    aln: cogent3 alignment object
        alignment object with keys 'structure' and 'reference'
    start: int
        start position (ones based counting) of where structure starts in the reference
    end: int
        end position (ones based counting) of where structure ends in the reference
    range: str
        a printable range of the start-end 
    pos_diff: int
        the difference between structure and reference based counting, according to the ref_start given
    diffList: list
        a list of mutations denoting aa1POSaa2, aa1 the original amino acid, pos the reference-based numbering, aa2 the mutation
    gaps_ref_to_struct: list
        a list of gaps from the reference frame of view (i.e. where reference has insertions, relative to structure) 
    gaps_struct_to_ref:
        a list of gaps from the strcture frame of view (i.e. where stucture has insertions, relative to reference)
    """
    #pscout here is reference instead
    fromName = 'structure'
    toName = 'reference'
    struct_sequence_ref_spanning = struct_sequence[pdb_start-1:pdb_start+ref_length+1]
    aln  = alignmentTools.returnAlignment(struct_sequence_ref_spanning, ref_sequence, fromName, toName)
    mapToRef = aln.get_gapped_seq('structure').gap_maps()[1]
    fromSeqValues = list(mapToRef.keys())
    from_start = min(fromSeqValues)
    from_end = max(fromSeqValues)
    #print("Range is %d to %d"%(from_start, from_end))
    range = "%d-%d"%(from_start+1, from_end+1)

    pos_diff = ref_start-(from_start+1) #move from_start to ones-base 
    diffList = alignmentTools.findDifferencesBetweenPairs(aln, from_start, from_end, ref_start, toName, fromName)

    #print("DEBUG: ref_start=%s, pdb_start=%d, ref_length=%d, length of subseq=%d, from_start_mapped=%d, pos_diff=%d"%(ref_start, pdb_start, ref_length, len(struct_sequence_ref_spanning), from_start, pos_diff))
    #print("DEBUG: fromSeqValues")
    #print(fromSeqValues)
    #if diffList:
        #print("Diff between Structure and Reference: %s"%(';'.join(diffList)))

    gaps_ref_to_struct = aln[from_start:from_end].count_gaps_per_seq()[fromName]
    gaps_struct_to_ref =  aln[from_start:from_end].count_gaps_per_seq()[toName]
    
    return aln, struct_sequence_ref_spanning, from_start+1, from_end+1, range, pos_diff, diffList, gaps_ref_to_struct, gaps_struct_to_ref


def return_domains_tuple(domain_str, INTERPRO):
    """
    Given the domain string from a reference file, split this into tuples

    Parameters
    ----------
    domain_str: str
        domain string that is ';' separated with domain name:start:stop
    INTERPRO: boolean
        If True, assumes the string to parse is INTERPRO (contains an extra ID), else it's Uniprot
    
    Returns
    -------
    domain_tuple: list
        Tuples of [domain name, start, stop] if INTERPRO=False
        Tuples of [domain name, Interpro_ID, start, stop] if INTERPRO=True
    """

    domain_list = domain_str.split(';')
    domain_tuple = []
    for domain_vals in domain_list:
        if INTERPRO:
            domain_name, Interpro_ID, start, stop = domain_vals.split(':')
            domain_tuple.append([domain_name, Interpro_ID, int(start), int(stop)])

        else:
            domain_name, start, stop = domain_vals.split(':')
            domain_tuple.append([domain_name, int(start), int(stop)])

    return domain_tuple

def returnDomainStruct(aln, ref_start, ref_stop, domains, diffList, domainExclValue=10):
    """
    Given a cogent3 alignment, from return_mapping_between_sequences and a list of domains in reference, test whether
    the alignment is of good enough quality in the region and return a dictionary of mapped elements. 
    Domains is a list of tuples in the form of [[]'Name', 'start', 'stop'], ['Name2', 'start', 'stop'], ..], where start and stop
    define the region of the domains in the reference sequence space. 
    
    Parameters
    ----------
    aln: cogent3.make_aligned_seqs
        Has two sequences, first entry is name1, aln_seq_1 and other is name2, aln_seq_2
    domains: list of tuples
        List of domains in format [('domainName', 'start', 'stop')] 
        or [('domainName', 'InterPro_ID', 'start', 'stop)]
        this is relative to the fromName sequence (i.e. name1 or name2 used in alignment creation)
    diffList: list
        list of mutations in aaPosaa 
    domainExclValue: int
        The window where you would allow a domain to map within this range of the start and stop
    
    Returns
    -------
    domainStruct: dict
        This dictionary has the domain tuples and maps to a tuple of [toName_sequence_start, toName_sequence_stop, numGaps, domain_num]
        using a unique key according to the start position of the domain in the alignment.
        if INTERPRO is detected then the tuple is [Interpro_ID, 
        toName_sequence_start, toName_sequence_stop, numGaps, domain_name]
        numGaps says how many gaps existed in the alignment. 
        Returns -1 if the region could not be mapped. It returns an empty dictionary if the alignment did not meet a gap threshold of less than 30%
    """
    ref_start = int(ref_start)
    ref_stop = int(ref_stop)
    fromName = 'structure'
    toName = 'reference'
    #check how to map from positions in struct to alignment
    #seq_to_aln_map = aln.get_gapped_seq(toName).gap_maps()[0]
    pscout_to_aln_map = aln.get_gapped_seq(fromName).gap_maps()[0]
    mapToStruct = aln.get_gapped_seq(toName).gap_maps()[1]
    #print(mapToStruct)
    INTERPRO=False
    mutPositions = returnMutPosList(diffList)
    gap_threshold = 0.7
    domainStruct = {}
    for domain in domains:
        #first check to see if the domain is within the range of the mapping. 
        if len(domain)==4:
            INTERPRO = True
            domain_name = domain[0]
            interpro_ID = domain[1]
            start_domain = int(domain[2])-1
            stop_domain = int(domain[3])-1

        else:
            domain_name = domain[0] #Uniprot indexes positions by 1, so have to remove
            start_domain = int(domain[1])-1
            stop_domain = int(domain[2])-1
        #print("DEBUG: checking if domain %d start and %d stop is between positions %d and %d"%(start_domain, stop_domain, ref_start, ref_stop))
        if start_domain < (ref_start-domainExclValue):
            #print("\tit is not")
            continue
        if stop_domain > (ref_stop + domainExclValue):
            #print("\tit is not")
            continue
        
        #print("\tFound domain")
        if start_domain < ref_start:
            start = ref_start
        else:
            start = start_domain
        if stop_domain > ref_stop:
            stop = ref_stop
        else:
            stop = stop_domain
  
        #domain = aln.get_seq(fromName).add_annotation(Feature, 'domain', domain_name, [(start, end)])
        
        #There is a very infrequent, but possible event, that a domain boundary exists at a gap, leading to a key issue
        if start not in mapToStruct:
            #how about is a start +/-1 OK?
            if start-1 in mapToStruct:
                start = start-1
            elif start+1 in mapToStruct:
                start = start+1
            else: 
                print("ERROR: domain start boundary not found due to gaps")
        
        if stop not in mapToStruct:
            #how about is a stop +/-1 OK?
            if stop+1 in mapToStruct:
                stop = stop+1
            elif stop-1 in mapToStruct:
                stop = stop-1
            else: 
                print("ERROR: domain stop boundary not found due to gaps")
        

        
        start_aln = mapToStruct[start]
        stop_aln = mapToStruct[stop]
        numGaps = aln.seqs[1][start_aln:stop_aln].count_gaps()
        numMuts = 0
        if(diffList):
            numMuts = countMutsInRange(mutPositions, start, stop) 
        #count the number of gaps in each sequence:
    # print('Domain: %s\t Length: %d \t Num Gaps: %d'%(domain_name, stop_aln-start_aln, numGaps ))
        if (numGaps/(stop_aln-start_aln) <= gap_threshold):
            start_aln_val = start_aln
            stop_aln_val = stop_aln
            try:
                if start_aln_val in domainStruct:
                    print("ERROR: more than one domain with the same start value was found. First will be overwritten")
                if not INTERPRO:
                    #domainStruct[domain_name] = [mapToStruct[start_aln_val]+1,  mapToStruct[stop_aln_val]+1, numGaps, numMuts]
                    domainStruct[start_aln_val] = [mapToStruct[start_aln_val]+1,  mapToStruct[stop_aln_val]+1, numGaps, numMuts, domain_name]
                else:
                    #domainStruct[domain_name] = [interpro_ID, mapToStruct[start_aln_val]+1,  mapToStruct[stop_aln_val]+1, numGaps, numMuts]
                    domainStruct[start_aln_val] = [interpro_ID, mapToStruct[start_aln_val]+1,  mapToStruct[stop_aln_val]+1, numGaps, numMuts, domain_name]
                #domain_name may already exist, so add a new value to it
            except:
                #print("ERROR: could not map domains")
                return -1
    return domainStruct


def countMutsInRange(posList, start, end):
    """
    Given a region of a protein, defined by start, end - count the number of mutations in that region
    This assumes all sequencing is in the same base (i.e. 0-based or 1-based counting)

    Parameters
    ----------
    posList: list of ints
        sorted integer list of ints from returnMutPosList(diffList)
    start: int
        position of start of range
    end: int
        position of end of range

    Returns
    -------
    numMuts: int
        number of mutations found in range start to end

    """

    numMuts = 0
    i = bisect.bisect_left(posList, start)
    g = bisect.bisect_right(posList, end)
    #if i != len(posList) and g != len(posList):
    return len(posList[i:g])
    #else return 0

    

def returnMutPosList(diffList):
    """
    Given a diffList (such as found in return_mapping_between_sequences), return just the positions that have mutations

    Parameters
    ----------
    diffList: list
        list of strings, each string is <aa><pos><aa>


    Returns
    -------
    mut_positions: list
        list of ints, just the positions mutations exist in, sorted

    """
    mut_positions = []
    for mutation in diffList:
        pos = ''.join(c for c in mutation if c.isdigit())
        mut_positions.append(int(pos))
    mut_positions.sort()
    return mut_positions


def return_reference_information(reference_df, uniprot_id, struct_seq, ref_seq_pos_start, pdb_pos_start, ref_length, INTERPRO):
    """
    Given inforation from a structure reference line, for a uniprot_id, the structure sequence
    and the mapped reference position, return the string-based information for appending to the
    structure file, including domains, mutations, sequence range, etc.

    Parameters
    ----------
    reference_df: pandas DataFrame
        loaded from a reference file
    uniprot_id: str
        uniprot ID of the structure sequence of interest
    struct_seq: str
        the structure sequence, includes mutations and modifications translated to single aa codes
    ref_seq_pos_start: int  
        position in reference of structure where experimental structure begins coverage
    pdb_pos_start: int
        position of structure where the reference position begins, may be different than 1 when tags exist experimentally
    INTERPRO: boolean   
        If True, uses the Intepro domain information, otherwise use Uniprot
    Returns
    -------
    gene_name: str
        The gene name found in the Uniprot Reference (reports N/A (Not Found in Reference) if this is not a domain-containing sequence found in the reference)
    struct_seq_ref_spanning: str
        ?
    rangeStr: str
        A string that denotes the start-end in the reference that was captured in the experiment. 
        For example 30-150 indicates amino acids 30 to 150 of the reference appeared in the structure.
    pos_diff: int
        The integer offset between the number of the reference and the numbering indicated in the structure file.
    diffStr: str
    
    gaps_ref_to_struct: int
        Number of insertions in the full range covering the reference that exist in strucutre that are not in reference
    gaps_struct_to_ref: int
        Number of insertions in the full range covering the reference that exist in reference that are not in structure
    domainStr: str
        A domain information string that follows the following rules
        Domains are separated by ';'
        Each domain set is separated by ':' and includes 'name:InterproID:info_string'
        info_string is ',' separated and has 4 entries 
            start - start position of domain IN STRUCTURE SEQUENCE NUMBERING
            end - end position fo the domain IN STRUCTURE SEQUENCE NUMBERING
            number_variations - number of variant positions (positions that differ between structure and reference) within the domain spanning region
            number_gaps - number of gaps/insertions within the domain range.
    structure_arch: str
        This is an easily readable string of domain names, in the order they appear in the protein from N- to C-terminal that was
        covered by the experiment. Domains are separated by '|'. For example SH3_domain|SH2|Prot_kinase_dom means that the 
        structure fully covered the SH3_domain, SH2, and Protein_kinase_dom of a SRC family kinase. 
    full_protein_domain: str
        For easy reference, this includes the domains and domain boundaries of the full protein. If domainStr does not contain 
        the protein domain, it's due to lack of full inclusion. 
    full_domain_arch: str
        For easy reference, this includes the architecture of the full protein so it can be seen what part of the whole protein
        was studied in the experiment.


    """


    diffStr = 'N/A'
    domainStr = 'N/A'
    gene_name = 'N/A (Not Found In Reference)'
    pos_diff = 0
    domainStr = 'N/A'
    structure_arch = 'N/A'
    full_protein_domain = 'N/A'
    full_domain_arch = 'N/A'


    if ref_seq_pos_start == 'not found': 
        # if not found, then we will set the ref_struct range and position

        rangeStr = '1-'+str(len(struct_seq)) #set the range starting from start to end
        
        struct_seq_ref_spanning = struct_seq #what it is based on existing.
        return gene_name, struct_seq_ref_spanning, rangeStr, pos_diff, diffStr, 0, 0, domainStr, structure_arch, full_protein_domain, full_domain_arch
    else:
        ref_seq_pos_start = int(ref_seq_pos_start)
        ref_length = int(ref_length)
        pdb_pos_start = int(pdb_pos_start)
    rangeStr = str(ref_seq_pos_start)+'-'+str(ref_seq_pos_start+ref_length) #start the default, to be updated later 
    struct_seq_ref_spanning = struct_seq[pdb_pos_start-1:pdb_pos_start+ref_length+1]

    #First find the protein information in the reference file based on uniprot_id
    protein_rec = reference_df[reference_df['UniProt ID']==uniprot_id]
    if len(protein_rec.index) < 1 or uniprot_id == 'not found':
        #print("NOTE: Encountered Uniprot %s in PDB, not in reference"%(uniprot_id))
        #return default information here
    
        return gene_name, struct_seq_ref_spanning, rangeStr, pos_diff, diffStr, 0, 0, domainStr, structure_arch, full_protein_domain, full_domain_arch
    elif len(protein_rec.index) > 1:
        print("ERROR: Found more than one record for %s in reference"%(uniprot_id))
    else:
        if INTERPRO:
            domains = list(protein_rec['Interpro Domains'])[0]
            full_protein_domain = domains
            full_domain_arch = list(protein_rec['Interpro Domain Architecture'])[0]

        else:
            domains = list(protein_rec['Uniprot Domains'])[0]
            full_protein_domain = domains
            full_domain_arch = list(protein_rec['Uniprot Domain Architecture'])[0]

            
        domain_tuple = return_domains_tuple(domains, INTERPRO) #domain_tuple is a different element if INTERPRO (includes ID)
        reference_seq = list(protein_rec['Ref Sequence'])[0]
        gene_name = list(protein_rec['Gene'])[0]
    
    aln, struct_seq_ref_spanning, from_start, from_end, rangeStr, pos_diff, diffList, gaps_ref_to_struct, gaps_struct_to_ref = return_mapping_between_sequences(struct_seq, reference_seq, ref_seq_pos_start, pdb_pos_start, ref_length)
    domainStruct = returnDomainStruct(aln, from_start, from_end, domain_tuple, diffList)
    #make the domainStr
    domainList = []
    domainDict_forArch = {}
    if isinstance(domainStruct, dict):
        for start_val in domainStruct:
            if INTERPRO:
                interproID, start, end, numGaps, numMuts, domain_name = domainStruct[start_val]
            else:
                start, end, numGaps, numMuts, domain_name = domainStruct[start_val]
            valStr = str(start)+','+str(end)+','+str(numGaps)+','+str(numMuts)
            domainDict_forArch[start] = domain_name
            if INTERPRO:
                domainList.append(domain_name+':'+ interproID + ':'+valStr)
            else:
                domainList.append(domain_name+':'+valStr)
    domainStr = ';'.join(domainList)

    arch_list = []
    starts = sorted(list(domainDict_forArch.keys()))
    for start in starts:
        arch_list.append(domainDict_forArch[start])
    structure_arch = '|'.join(arch_list)


    diffStr = ";".join(diffList)
    
    #make gaps 
    
    return gene_name, struct_seq_ref_spanning, rangeStr, pos_diff, diffStr, gaps_ref_to_struct, gaps_struct_to_ref, domainStr, structure_arch, full_protein_domain, full_domain_arch

def add_reference_info_to_struct_file(struct_file, ref_file, out_file, INTERPRO=True, verbose=False):
    """
    Given a PDB meta structure file and a Uniprot reference, integrate the two pieces to add information from reference
    
    Parameters
    ----------
    struct_file: str
        Name of structure reference file
    ref_file: str 
        Name of reference file
    out_file: str   
        name of output file to write
    INTERPRO: boolean
        If True, uses Interpro, otherwise appends Uniprot from reference file. 
        Recommended behavior is to use Interpro - it is more inclusive of domain boundaries and has better naming
        conventions, along with perserving ability to use the Interpro ID for filtering strucutres containing domains of interest. 
    verbose: boolean
        Print information about processing. Default is False.
    
    Returns
    -------
    out_struct: pandas dataframe
        the appended dataframe of the structure (also written to out_file)
    """
    struct_df = pd.read_csv(struct_file)
    reference_df = pd.read_csv(ref_file)

    for index, row in struct_df.iterrows():
        uniprot_id = row['database_accession']
        struct_seq = row['pdbx_seq_one_letter_code_can']
        ref_seq_pos_start = row['ref_beg_seq_id']
        pdb_seq_pos_start  = row['entity_beg_seq_id']
        ref_length = row['aligned_regions_length']
        if verbose:
            print("Working on %s for protein %s"%(row['PDB_ID'], uniprot_id))
        information_list = return_reference_information(reference_df, uniprot_id, struct_seq, ref_seq_pos_start, pdb_seq_pos_start, ref_length, INTERPRO)
        gene_name, struct_seq_ref_spanning, rangeStr, pos_diff, diffStr, gaps_ref_to_struct, gaps_struct_to_ref, domainStr, structure_arch, full_protein_domain, full_domain_arch = information_list
        struct_df.loc[index,'ref:gene name'] = gene_name
        struct_df.loc[index, 'ref:struct/ref sequence'] = struct_seq_ref_spanning
        struct_df.loc[index, 'ref:reference range'] = rangeStr
        struct_df.loc[index, 'ref:start position offset'] = pos_diff
        struct_df.loc[index, 'Gaps ref:struct'] = gaps_ref_to_struct
        struct_df.loc[index, 'Gaps struct:ref'] = gaps_struct_to_ref
        struct_df.loc[index, 'ref:variants'] = diffStr
        struct_df.loc[index, 'ref:domains'] = domainStr
        struct_df.loc[index, 'ref:struct domain architecture'] = structure_arch
        struct_df.loc[index, 'ref:full protein domain'] = full_protein_domain
        struct_df.loc[index, 'ref:protein domain architecture'] = full_domain_arch
    struct_df.to_csv(out_file, index=False)
    return struct_df


def filter_structure_file(appended_structure_file, Interpro_ID, filtered_structure_file):
    """
    Given an annotated structure file, keep only structures that have at least one chain that contain the 
    Interpro_ID of interest. 

    Prints the filtered structure file to filtered_structure_file

    Parameters
    ----------
    appended_structure_file: str
        location of the UniProt reference file that also has been appended using InterPro.appendRefFile to add
        interpro structures. 
    Interpro_ID: str
        Interpro ID (controlled Interpro ID database identifier)
    filtered_structure_file: str
        location to write the output file to - all the same column fields, but reducing the rows to only those
        that contain an Interpro_ID of interest. 

    """

    PDB_df = pd.read_csv(appended_structure_file)
    filtered_PDB_list = []
    total_PDBs = 0
    for name, group in PDB_df.groupby('PDB_ID'):
        total_PDBs+=1
        for index, row in group.iterrows():
            if isinstance(row['ref:domains'], str) and Interpro_ID in row['ref:domains']:
                filtered_PDB_list.append(name)
    PDB_df_sub = PDB_df[PDB_df['PDB_ID'].isin(filtered_PDB_list)]
    PDB_df_sub.to_csv(filtered_structure_file, index=False)
    print("Made %s file: %d structures retained of %d total"%(filtered_structure_file, len(filtered_PDB_list), total_PDBs))
        
