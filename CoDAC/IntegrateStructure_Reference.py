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