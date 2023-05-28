from CoDAC import alignmentTools


def return_mapping_between_sequences(struct_sequence, ref_sequence, ref_start):
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
    aln  = alignmentTools.returnAlignment(struct_sequence, ref_sequence, fromName, toName)
    mapToRef = aln.get_gapped_seq('structure').gap_maps()[1]
    fromSeqValues = list(mapToRef.keys())
    from_start = min(fromSeqValues)
    from_end = max(fromSeqValues)
    #print("Range is %d to %d"%(from_start, from_end))
    range = "%d-%d"%(from_start+1, from_end+1)

    pos_diff = mapToRef[from_start] - (ref_start-1) #checking the frame of reference   
    diffList = alignmentTools.findDifferencesBetweenPairs(aln, from_start, from_end, ref_start, toName, fromName)

    #if diffList:
        #print("Diff between Structure and Reference: %s"%(';'.join(diffList)))

    gaps_ref_to_struct = aln[from_start:from_end].count_gaps_per_seq()[fromName]
    gaps_struct_to_ref =  aln[from_start:from_end].count_gaps_per_seq()[toName]
    
    return aln, from_start+1, from_end+1, range, pos_diff, diffList, gaps_ref_to_struct, gaps_struct_to_ref


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

def returnDomainStruct(aln, domains):
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

    Returns
    -------
    domainStruct: dict
        This dictionary has the name of domainName tuples and maps to a tuple of [toName_sequence_start, toName_sequence_stop, numGaps]
        numGaps says how many gaps existed in the alignment. 
        Returns -1 if the region could not be mapped. It returns an empty dictionary if the alignment did not meet a gap threshold of less than 30%
    """

    fromName = 'structure'
    toName = 'reference'
    #check how to map from positions in struct to alignment
    #seq_to_aln_map = aln.get_gapped_seq(toName).gap_maps()[0]
    pscout_to_aln_map = aln.get_gapped_seq(fromName).gap_maps()[0]
    mapToStruct = aln.get_gapped_seq(toName).gap_maps()[1]

    gap_threshold = 0.7
    domainStruct = {}
    for domain in domains:
        domain_name = domain[0] #ProteomeScout indexes positions by 1, so have to remove
        start = int(domain[1])-1
        stop = int(domain[2])-1
        #domain = aln.get_seq(fromName).add_annotation(Feature, 'domain', domain_name, [(start, end)])

        start_aln = pscout_to_aln_map[start]
        stop_aln = pscout_to_aln_map[stop]
        numGaps = aln.seqs[1][start_aln:stop_aln].count_gaps()
        #count the number of gaps in each sequence:
       # print('Domain: %s\t Length: %d \t Num Gaps: %d'%(domain_name, stop_aln-start_aln, numGaps ))
        if (numGaps/(stop_aln-start_aln) <= gap_threshold):
            start_aln_val = start_aln
            stop_aln_val = stop_aln
            #if start_aln > min(mapToStruct.keys()):
            #    start_aln_val = min(mapToStruct.keys())
            #if stop_aln > max(mapToStruct.keys()):
            #    stop_aln_val = max(mapToStruct.keys())
            #print("Domain: %s is MATCHED: start=%d, stop=%d"%(domain_name, mapToStruct[start_aln_val], mapToStruct[stop_aln_val]))
            try:
                if domain_name in domainStruct:
                    domain_name = "%s_%d"%(domain_name, 2)
                domainStruct[domain_name] = [mapToStruct[start_aln_val]+1,  mapToStruct[stop_aln_val]+1, numGaps]
                #domain_name may already exist, so add a new value to it
            except:
                #print("ERROR: could not map domains")
                return -1
    return domainStruct