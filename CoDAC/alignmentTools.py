
from cogent3 import make_aligned_seqs

from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def returnAlignment(seq1, seq2, name1, name2):
    """
    Use pairwise2 alignment and return alignment as a cogent3 alignment class

    Parameters
    ----------
    seq1: str
        first sequence
    seq2: str
        second sequence
    name1: str
        reference name for seq1 to be used in accessing the aln object
    name2: str
        reference name for seq1 to be used in accessing the aln object

    Returns
    -------
    aln: cogent3.make_aligned_seqs
        Has two sequences, first entry is name1, aln_seq_1


    """
    alignments = pairwise2.align.localms(seq1, seq2, 2, -1, -1, -0.5)
    if not alignments:
        #ERROR in alignment, return 0
            return 0
    
    aln_seq_1 = alignments[0][0]
    aln_seq_2= alignments[0][1]

    aln = make_aligned_seqs([[name1,aln_seq_1], [name2, aln_seq_2]], array_align=True) #cogent3 pairwise alignment object
    return aln

def findDifferencesBetweenPairs(aln, region_start, region_stop, offset, fromName, toName):
    """
    Given an alignment object, return the differences that exist in the region of fromName to toName for a specific region of the alignment
    
    """
    
    #mapToRef = aln.get_gapped_seq(fromName).gap_maps()[0]
    aln_tuples = aln.positions[region_start:region_stop]   
    diffList = []
    for i in range(0, len(aln_tuples)):
        pos = aln_tuples[i]
        pos_in_protein = i+offset 

        if pos[0] != pos[1]:
                #print("FOUND MuTATION")
                diffList.append("%s%d%s"%(pos[0], pos_in_protein, pos[1]))
    return diffList
