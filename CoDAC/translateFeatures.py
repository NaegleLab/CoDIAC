import cogent3
from Bio import pairwise2
from Bio import SeqIO
import numpy as np
import pandas as pd

'''
These tools allow you to translate features from an input fasta version of protein or protein coding regions onto a reference.
'''

def read_feature_file(feature_file):
    feature_file='lipid_features.txt'
    feature_dict={}
    type_dict = {}
    with open(feature_file, 'r') as file:
        for row in file:
            columns = row.split('\t')
            if len(columns) > 2:
                header = columns[1]
                feature = int(columns[4])
                if header not in feature_dict:
                    feature_dict[header] = []
                feature_dict[header].append(feature)
            elif len(columns) == 2:
                feature_type = columns[0]
                color = columns[1]
                type_dict[feature_type] = color
    return type_dict, feature_dict

def return_sequence_dict_from_fasta(fasta_file):
    """
    Returns a dictionary with keys equal to the fasta headers in a file and values equal to the protein sequence
    reports if redundant fasta header is found by returning that in a error_list

    Parameters
    ----------
    fasta : str
        Name of fasta file

    Returns
    --------
    seq_dict: dict
        Keys are the fasta header of the fasta file. Values or protein sequences (strings)
    error_list: list
        A list of headers that had more than one entry 
    """
    seqIO_obj = list(SeqIO.parse(open(fasta_file), "fasta"))
    seq_dict = {}
    error_list = []
    for record in seqIO_obj:
        if record.id not in seq_dict: 
    #for those cases where more than one SH2 domain exists, need to reference that according to start position in protein
            seq_dict[record.id] = str(record.seq)
        else:
            error_list.append(record.id)
    return seq_dict, error_list




def returnAlignment(seq1, seq2):
    """
    Given two protein sequences and the names for reference, create an alignment using BioSeq pairwise alignment
    
    Parameters
    ----------
    seq1 : str
        first sequence
    seq2: str
        second sequence
    name1: str
        name reference for hanlde to seq1
    name2: str
        name reference for handle to seq2
    
    Returns
    --------
    aln: cogent3 pairwise alignment object
        The alignment object of seq1, name1 with seq2, name2    
    
    """
    alignments = pairwise2.align.localms(seq1, seq2, 2, -1, -1, -0.5)
    if not alignments:
        #ERROR in alignment, return 0
            return 0
    aln_seq_1 = alignments[0][0]
    aln_seq_2= alignments[0][1]
    aln = cogent3.make_aligned_seqs([['Input',aln_seq_1], ['Reference', aln_seq_2]], array_align=True) #cogent3 pairwise alignment object
    return aln


def return_best_match(input_header, input_sequence, reference_dict):
    """
    Given an input header and sequence, find the best match inside a reference_dict (keys are headers and values are sequences from a reference fasta)
    Will try quick search first, assumes the same ID is used at the front of the headers. Otherwise, performs a brute-force all by all pairwise alignment

    Parameters
    ----------
    input_header: str
        Input header, string separated by | for fields. Uses the [0] element in quick search by split('|')
    input_sequence: str
        input sequence
    reference_dict: dict
        Reference fasta converted to dictionary, keys are fasta headers and values are sequences
    
    Returns
    -------
    best_header: str
        The header of the best match in reference_dict
    best_aln: cogent3 alignment object
        alignment of 'Input' to 'Reference' of best match
    aln_score: float
        score of alignment
    percent_matched:
        percent of input sequence that matches the reference (excludes gaps)
    percent_gapped:
        percent of input sequence with gaps to reference
    
    """

    header_vals = input_header.split('|')
    match_headers = return_headers_with_matching_ID(header_vals[0], reference_dict)
    aln_dict = generate_alignment_dict(match_headers, input_sequence, reference_dict)
    best_header, best_aln, aln_score = return_best_match_from_set(aln_dict)
    percent_matched, percent_gapped = return_percent_matched(best_aln)
    return best_header, best_aln, aln_score, percent_matched, percent_gapped



def return_headers_with_matching_ID(ID, reference_dict):
    """
    For a string ID (e.g. an accession number like 'P00533'), find all sequences in seq_dict that have that common ID at the start of the headers
    Requires an exact match of ID to the first value of fasta header, where fasta headers are separated by '|'

    Parameters
    ----------
    ID: string
        The common ID we are looking for at start of fasta headers
    seq_dict: dict
        Keys are the fasta header of the fasta file. Values or protein sequences (strings)
    
    Returns
    -------
    match_list: list
        List of headers in seq_dict that have common ID
    """
    match_list = []
    for fasta_header in reference_dict:
        header_items = fasta_header.split('|')
        if header_items[0]==ID:
            match_list.append(fasta_header)
    return match_list


def generate_alignment_dict(match_list, input_sequence, reference_dict):
    """
    If there are matches in the match_list for the given input_sequence, this generates an alignment dict of those 
    sequences with the input_sequence. If there are no matches identified, it generates an alignment_dict of all of the 
    reference sequences aligned with the input_sequence.
    
    Returns
    -------
    aln_dict: dictionary of alignments between the input_sequence and either the sequences in the match_list, or all 
    sequences in the reference_dict.
    """

    aln_dict={}
    if len(match_list):
        for ref_header in match_list:
            aln_dict[ref_header]=returnAlignment(input_sequence,reference_dict[ref_header])
    else:
        for ref_header in reference_dict:
            aln_dict[ref_header]=returnAlignment(input_sequence,reference_dict[ref_header])
    return aln_dict



def match_and_align(input_header, input_sequence, reference_dict):
    """
    For a string ID (e.g. an accession number like 'P00533'), find all sequences in seq_dict that have that common ID at the start of the headers
    Requires an exact match of ID to the first value of fasta header, where fasta headers are separated by '|'
    
    Then, If there are matches in the match_list for the given input_sequence, this generates an alignment dict of those 
    sequences with the input_sequence. If there are no matches identified, it generates an alignment_dict of all of the 
    reference sequences aligned with the input_sequence.

    If the fast match cannot be made, then will return 
    
    Parameters
    ----------
    ID: string
        The common ID we are looking for at start of fasta headers
    input_sequence: string
        the protein sequence that will be used to align to proteins identified from the match_list
    reference_dict: dict
        Keys are the fasta header of the fasta file. Values or protein sequences (strings)
    
    Returns
    -------
    match_list: list
        List of headers in seq_dict that have common ID
        
    matched_dict: dict
        alignments between the input_sequence and either the sequences in the match_list, or all 
    sequences in the reference_dict.
    """
    for header in feature_dict:
        match_list = []
        matched_dict={}
        values = header.split('|')
        ID = values[0]
        for fasta_header in reference_dict:
            header_items = fasta_header.split('|')
            if header_items[0]==ID:
                match_list.append(fasta_header)
        #align based on match_list results
        if len(match_list):
            for ref_header in match_list:
                matched_dict[ref_header]=returnAlignment(input_sequence,reference_dict[ref_header])
        else:
            for ref_header in reference_dict:
                matched_dict[ref_header]=returnAlignment(input_sequence,reference_dict[ref_header])
    return(matched_dict,match_list)


def return_best_match_from_set(matched_dict):
    
    """
    Given an input header, alignment dictionary, and the positions of the features in the given input protein, this will systematically         align the input sequence with all options identified previously and return the header, alignment object, and alignment quality that         best matched. This will also give the percent similarity for the input header and the corresponding reference header. Lastly, this will     give the positions of any positions that are not able to be mapped as a feature due to mutation or a gap in the sequence.
    
    Parameters
    -------
    alignment_dict: dict
        dictionary of all of the alignment objects corresponding to your input header
    gaps: arr
        array of which positions correspond to gaps. This is important so that these gaps do not affect the sequence alignment 
    
    Returns
    -------
    best_header: string
        The best matched header based on alignment for the given input_sequence
    aln_score: int
        The alignment score for the best aligned pair
    percent_match: int
        The percent similarity value between 0-100
    """
    best_header=''
    aln_score=0
    for header in matched_dict:
        aln_score_temp=matched_dict[header].alignment_quality()
        if aln_score_temp > aln_score:
            aln_score = aln_score_temp
            best_header=header
            best_aln=matched_dict[best_header]

    #Return percent matched here as well. 



    return best_header, best_aln, aln_score


def return_percent_matched(aln):
    """
    Given a cogent3 alignment object, calculate the percent matched (exact amino acid matches of the input sequence to the comparator)
    and percent gapped

    Parameters
    ----------
    aln: cogent3 alignment object
        assumes two sequences, 'Input' and 'Reference'
    
    Returns
    -------
    percent_matched: float
        percent matched, excluding gaps
    percent_gapped: float
        percent of input that is gapped, relative to reference
    

    """
    total_gaps = aln.count_gaps_per_seq()['Input']
    input_length = aln.get_lengths()['Input']
    seq_input = aln.get_gapped_seq('Input')
    seq_reference = aln.get_gapped_seq('Reference')
    matches = sum(a == b and a!='-' for a, b in zip(seq_input, seq_reference))
    percent_matched = 100*matches/input_length
    percent_gaps = 100*total_gaps/input_length
    return(percent_matched, percent_gaps)

def return_translation_dict_for_all_feature_seqs(feature_dict, input_dict, reference_dict):
    """
    Given a set of unique features that we need to translate from an input feature file (see ead_feature_file)
    Return a dictionary with input fasta headers and a dictionary type that includes the matching output header, alignment and quality scores

    Parameters
    ----------
    feature_dict: dict
        keys are fasta headers, values are the list of features to translate
    input_dict: dict
        input dictionary, keys fasta headers and values are sequences
    reference_dict: dict
        reference dictionary, keys fasta headers, and values are sequences
    
    Returns
    -------
    feature_trans_dict: dict
        Dict of dicts, outer key is the header and inner keys are the following special keys
            'header': header to output
            'aln': alignment object of best header
            'percent_matched': percent of input sequence that matches, not including gaps
            'percent_gapped': percent of input sequence that is gapped 
            'aln_score': alignment score, calculated by cogent3
    
    """
    feature_trans_dict = {}
    for fasta_header in feature_dict:
        header, aln, aln_score, percent_matched, percent_gapped = return_best_match(fasta_header,  input_dict[fasta_header], reference_dict)
        temp_dict = {}
        temp_dict['header'] = header
        temp_dict['aln'] = aln
        temp_dict['percent_matched'] = percent_matched
        temp_dict['percent_gapped'] = percent_gapped
        temp_dict['aln_score'] = aln_score
        feature_trans_dict[fasta_header] = temp_dict
    return feature_trans_dict

def translate_features(matched_alns, feature_pos):
    """
    translating the features
    matched_alns = matched_dict[fasta_header]
    feature_pos = feature_dict[fasta_header]
    """
    validated_features=[]
    gaps=[]
    mutations=[]
    y=matched_dict[fasta_header]
    map_to_ref= y.get_gapped_seq('Input').gap_maps()[0]
    aa=list(y)
    for feature in feature_dict[fasta_header]:
        if aa[feature-1][0]==aa[feature-1][1]: #true features
            validated_features.append(feature)
            #print('Validated features:', validated_features)
        if aa[feature-1][0] != aa[feature-1][1] and aa[feature-1][0]=='-': 
            aa[feature-1][0]= 'Z'
    for x in range(0, len(feature_dict[fasta_header])-1):
        if aa[feature-1][0]=='Z' or aa[feature-1][0]=='-':
            gaps.append(feature_dict[fasta_header][x])
            print('Gap identified at position:', gaps[x])
        if aa[feature-1][0] != aa[feature-1][1] and aa[feature-1][0]!='-' and aa[feature-1][0]!='Z': #mutation
            mutations.append(feature)
            print('Mutation identified at position:', mutations[x])
    
    #percent match code 
    for gap in gaps:
        if aa[gap-1][0] != aa[gap-1][1] and aa[gap-1][0]=='-': 
            aa[gap-1][0]= 'Z'
    input_seq = [pair[0] for pair in aa]
    ref_seq=[pair[1] for pair in aa]
    seq_length = min(len(input_seq), len(ref_seq))
    match=0
    for i in range(0, seq_length):
        if input_seq[i] == ref_seq[i]:
            match=match+1
        i=i+1
    percent_match=(match/seq_length) *100

    if percent_match < 95:
        temp_dict = {}
        temp_dict['header'] = best_header
        temp_dict['aln'] = best_aln
        matched_dict[fasta_header] = temp_dict
        percent_match= 'ERROR: BELOW SIMILARITY THRESHOLD'
        del feature_dict[fasta_header][feature]
   
    print('The alignment score is:', aln_score)
    print('The validated features are at positions:', validated_features)
 
    return mutations,gaps, validated_features







        
          