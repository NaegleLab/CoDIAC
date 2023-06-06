from CoDAC import alignmentTools
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


def return_domains_tuple(domain_str):
    """
    Given the domain string from a reference file, split this into tuples

    Parameters
    ----------
    domain_str: str
        domain string that is ';' separated with domain name:start:stop
    
    Returns
    -------
    domain_tuple: list
        Tuples of [domain name, start, stop]
    """

    domain_list = domain_str.split(';')
    domain_tuple = []
    for domain_vals in domain_list:
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
        List of domains in format [('domainName', 'start', 'stop')] this is relative to the fromName sequence (i.e. name1 or name2 used in alignment creation)
    diffList: list
        list of mutations in aaPosaa 
    domainExclValue: int
        The window where you would allow a domain to map within this range of the start and stop
    
    Returns
    -------
    domainStruct: dict
        This dictionary has the name of domainName tuples and maps to a tuple of [toName_sequence_start, toName_sequence_stop, numGaps]
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

    mutPositions = returnMutPosList(diffList)
    gap_threshold = 0.7
    domainStruct = {}
    for domain in domains:
        #first check to see if the domain is within the range of the mapping. 
        domain_name = domain[0] #Uniprot indexes positions by 1, so have to remove
        start_domain = int(domain[1])-1
        stop_domain = int(domain[2])-1
        #print("DEBUG: checking if domain %d start is between positions %d and %d"%(start_domain, ref_start, ref_stop))
        if start_domain < (ref_start-domainExclValue):
            print("\tit is not")
            continue
        if stop_domain > (ref_stop + domainExclValue):
            print("\tit is not")
            continue
        
        if start_domain < ref_start:
            start = ref_start
        else:
            start = start_domain
        if stop_domain > ref_stop:
            stop = ref_stop
        else:
            stop = stop_domain
  
        #domain = aln.get_seq(fromName).add_annotation(Feature, 'domain', domain_name, [(start, end)])
        start_aln = mapToStruct[start_domain]
        stop_aln = mapToStruct[stop_domain]
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
                if domain_name in domainStruct:
                    domain_name = "%s_%d"%(domain_name, 2)
                domainStruct[domain_name] = [mapToStruct[start_aln_val]+1,  mapToStruct[stop_aln_val]+1, numGaps, numMuts]
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


def return_reference_information(reference_df, uniprot_id, struct_seq, ref_seq_pos_start, pdb_pos_start, ref_length):
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

    Returns
    -------

    """
    rangeStr = str(ref_seq_pos_start)+'-'+str(ref_seq_pos_start+ref_length) #start the default, to be updated later 
    diffStr = '-1'
    domainStr = ''
    gene_name = '-1'
    pos_diff = 0
    domainStr = '-1'
    structure_arch = '-1'
    full_domain_arch = '-1'
    struct_seq_ref_spanning = struct_seq[pdb_pos_start-1:pdb_pos_start+ref_length+1]

    #First find the protein information in the reference file based on uniprot_id
    protein_rec = reference_df[reference_df['UniProt ID']==uniprot_id]
    if len(protein_rec.index) < 1:
        print("ERROR: Did not find %s in reference"%(uniprot_id))
        #return default information here
        return gene_name, struct_seq_ref_spanning, rangeStr, pos_diff, diffStr, 0, 0, domainStr, structure_arch, full_domain_arch
    elif len(protein_rec.index) > 1:
        print("ERROR: Found more than one record for %s in reference"%(uniprot_id))
    else:
        domains = list(protein_rec['Domains'])[0]
        domain_tuple = return_domains_tuple(domains)
        full_domain_arch = list(protein_rec['Domain Architecture'])[0]
        reference_seq = list(protein_rec['Ref Sequence'])[0]
        gene_name = list(protein_rec['Gene'])[0]
    
    aln, struct_seq_ref_spanning, from_start, from_end, rangeStr, pos_diff, diffList, gaps_ref_to_struct, gaps_struct_to_ref = return_mapping_between_sequences(struct_seq, reference_seq, ref_seq_pos_start, pdb_pos_start, ref_length)
    domainStruct = returnDomainStruct(aln, from_start, from_end, domain_tuple, diffList)
    #make the domainStr
    domainList = []
    domainDict_forArch = {}
    if isinstance(domainStruct, dict):
        for domain_name in domainStruct:
            start, end, numGaps, numMuts = domainStruct[domain_name]
            valStr = str(start)+','+str(end)+','+str(numGaps)+','+str(numMuts)
            domainDict_forArch[start] = domain_name
            domainList.append(domain_name+':'+valStr)
    domainStr = ';'.join(domainList)

    arch_list = []
    starts = sorted(list(domainDict_forArch.keys()))
    for start in starts:
        arch_list.append(domainDict_forArch[start])
    structure_arch = '|'.join(arch_list)


    diffStr = ";".join(diffList)
    
    #make gaps 
    
    return gene_name, struct_seq_ref_spanning, rangeStr, pos_diff, diffStr, gaps_ref_to_struct, gaps_struct_to_ref, domainStr, structure_arch, full_domain_arch

def add_reference_info_to_struct_file(struct_file, ref_file, out_file, verbose='False'):
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
    verbose: boolean
        Print information about processing
    
    Returns
    -------
    out_struct: pandas dataframe
        the appended dataframe of the structure (also written to out_file)
    """
    struct_df = pd.read_csv(struct_file)
    reference_df = pd.read_csv(ref_file)

    for index, row in struct_df.iterrows():
        uniprot_id = row['ACCESS']
        struct_seq = row['CANNONICAL_REF_SEQ']
        ref_seq_pos_start = row['CANNONICAL_SEQ_BEG_POSITION']
        pdb_seq_pos_start  = row['PDB_SEQ_BEG_POSITION']
        ref_length = row['REF_SEQ_LENGTH']
        if verbose:
            print("Working on %s"%(uniprot_id))
        information_list = return_reference_information(reference_df, uniprot_id, struct_seq, ref_seq_pos_start, pdb_seq_pos_start, ref_length)
        gene_name, struct_seq_ref_spanning, rangeStr, pos_diff, diffStr, gaps_ref_to_struct, gaps_struct_to_ref, domainStr, structure_arch, full_domain_arch = information_list
        struct_df.loc[index,'gene name'] = gene_name
        struct_df.loc[index, 'struct/ref sequence'] = struct_seq_ref_spanning
        struct_df.loc[index, 'reference range'] = rangeStr
        struct_df.loc[index, 'pos diff'] = pos_diff
        struct_df.loc[index, 'Gaps Ref:Struct'] = gaps_ref_to_struct
        struct_df.loc[index, 'Gaps Struct:Ref'] = gaps_struct_to_ref
        struct_df.loc[index, 'mutations'] = diffStr
        struct_df.loc[index, 'domains'] = domainStr
        struct_df.loc[index, 'struct domain architecture'] = structure_arch
        struct_df.loc[index, 'protein domain architecture'] = full_domain_arch
    struct_df.to_csv(out_file, index=False)
    return struct_df



