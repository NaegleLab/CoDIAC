
from Bio import SeqIO

from cogent3 import load_aligned_seqs

import pandas as pd
import numpy as np

import sys


"""
CoDIAC.featureTools 
Assumes fasta and feature files, such as used by Jalview definitions.
Includes Functions to integrate feature files, create annotation tracks, and translate features between possible different sources
"""

COLORS = ['332288', '117733', '44AA99', '88CCEE', 'DDCC77', 'CC6677', 'AA4499', '882255'] #Paul Tol library for color blindness

def return_features_from_file(file):
    """
    Given a feature file, this will parse the feature set types and return feature_colors and seq_dict. If a combined feature is encountered,
    indicated by a name_feature1:name_feature2 (i.e. colon separated files), then we add these individual features to the list of 
    features observed

    Parameters
    ----------
    file : str
        Name of feature file

    
    Returns
    --------
    feature_colors: dict
        keys are the feature names and values are the color values from the feature file
    seq_dict: dict
        Dict of dicts, with parent key the fasta header seq id, then inner dict with keys equal to the amino acid position
        and values equal to a list of features set for that feature.

    """ 
    feature_colors = {}
    seq_dict = {}
    lines = open(file, 'r')
    for line in lines:
        line = line.strip()
        #check if the line is a preamble for a feature or a feature line
        line_arr = line.split('\t')
        if len(line_arr) == 2:
            feature_colors[line_arr[0]] = line_arr[1]
        elif len(line_arr) == 6:
            feature_desc, seq_id, offset, pos, pos, feature = line_arr
            if seq_id not in seq_dict:
                seq_dict[seq_id] = {}
            #now check if that position already has an entry 
            if pos not in seq_dict[seq_id]:
                #if feature is a colon separated, then add all the features of that type to seq_id
                feature_set = feature.split(':')
                seq_dict[seq_id][pos] = [feature]
                if len(feature_set) > 1:
                    for feat in feature_set:
                        seq_dict[seq_id][pos].append(feat)
        else:
            print("Skipping %s since it does not match a color preamble or a feature line"%(line))

    return feature_colors, seq_dict

def combine_feature_sets(feature_dict_1, feature_dict_2):
    """
    Given feature file sets and the colors used, combine these into a new feature set,
    and creating new features if needed. This sets a combined feature from two files for the same
    residue as the longest combination of features. 
    e.g. of K15 in feature_dict_1 has 'N6-acetyllysine:Ubiquitination' and in feature_dict_2 it is 'Methylation',
    the output feature will be 'N6-acetyllysine:Ubiquitination:Methylation' and it will have all three individual features, 
    so that annotation tracks can be created 

    Parameters
    ----------
    feature_dict_1 : dict
        Dict in the form returned by return_features_from_file(file). This is a dictionary of one feature file. 
        Dict of dicts, with parent key the fasta header seq id, then inner dict with keys equal to the amino acid position
        and values equal to a list of features set for that feature.
    feature_dict_2: dict
        Same as feature_dict_1, but generated from a second feature file

    
    Returns
    --------
    feature_combined_dict: dict
        Same output as input as dict of dictionaries, where the two features have been combined. 


    """
    # feature_lists = list(feature_dict_1.keys()).append(list(feature_dict_2.keys()))
    # unique_features = set(feature_lists)
    # feature_colors = {}
    # for feature in unique_features:
    #   feature_colors{feature} = '' #for now, let's just note that we have a color we need to set.

    feature_combined_dict = {}

    for seq_id in feature_dict_1.keys():
        # walk through each in feature_dict_1 and if it also exists in dict_2, then combine those feature types into new categories.
        #once that is taken care of, do not walk through it in dict_2 separately
        if seq_id not in feature_combined_dict:
            feature_combined_dict[seq_id] = {} #add all features from dict_1
            for pos in feature_dict_1[seq_id]:
                feature_combined_dict[seq_id][pos] = feature_dict_1[seq_id][pos]
                if seq_id in feature_dict_2:
                    if pos in feature_dict_2[seq_id]:
                        feature_combined_dict[seq_id][pos].extend(feature_dict_2[seq_id][pos])
                        #also need to make a new category which is the largest combined names of the two feature dictionaries
                        longest_feature_1 = return_longest_feature(feature_dict_1[seq_id][pos])
                        longest_feature_2 = return_longest_feature(feature_dict_2[seq_id][pos])
                        feature_combined_dict[seq_id][pos].append(longest_feature_1+':'+longest_feature_2)

    #now add features in 2 that did not overlap with positions in 1
    for seq_id in feature_dict_2.keys():
        if seq_id not in feature_combined_dict:
            feature_combined_dict[seq_id] = {}
        for pos in feature_dict_2[seq_id]:
            if pos not in feature_combined_dict[seq_id]:
                feature_combined_dict[seq_id][pos] = feature_dict_2[seq_id][pos]
                    



    return feature_combined_dict


def print_jalview_feature_file(feature_dict, feature_color_dict, feature_file, append=False):
    """ 
    Print all features in a feature_dict to a feature_file, based on the color mapping given in the feature_color_dict
    This ensures that the longest combined feature is printed for a residue that has multiple annotations.

    Parameters
    ----------
    feature_dict : dict
        Dict in the form returned by return_features_from_file(file) or combine_feature_sets(feature_dict_1, feature_dict_2). 
        Dict of dicts, with parent key the fasta header seq id, then inner dict with keys equal to the amino acid position
        and values equal to a list of features set for that feature.
    color_dict: dict
        Keys are the features in feature_dict and values are the hex code values to designate in the feature file.
    feature_file: str
        File name of string, this will be overwritten by default
    append: bool
        Default false, meaning the feature file will be overwritten. Otherwise, set to true to append new lines to an existing
        feature file.



    """

    if append:
        ff = open(feature_file, "a")
    else:
        ff = open(feature_file, "w")

    #first write the preamble lines of the colors 
    for feature in feature_color_dict:
        ff.write("%s\t%s\n"%(feature, feature_color_dict[feature]))

    #now write the lines for every feature (have to check if a color has not been given, randomly select from end of color wheel)
    for seq_id in feature_dict:
        for pos in feature_dict[seq_id]:
            feature = return_longest_feature(feature_dict[seq_id][pos])
            if feature not in feature_color_dict:
                print("ERROR: did not find color for %s"%(feature))
            ff.write("%s\t%s\t-1\t%s\t%s\t%s\n"%(feature, seq_id, pos, pos, feature))
    ff.close()


def return_unique_features(feature_dict):
    """
    Given a feature_dict, return all unique features.

    Parameters
    ----------
    feature_dict : dict
        Dict in the form returned by return_features_from_file(file) or combine_feature_sets(feature_dict_1, feature_dict_2). 
        Dict of dicts, with parent key the fasta header seq id, then inner dict with keys equal to the amino acid position
        and values equal to a list of features set for that feature.

    Returns
    -------
    features: list
        List of unique features
    feature_numbers: dict
        A dict of features as keys and number of features of that type observed
    """
        #set the colors
    feature_unique = {}

    for seq_id in feature_dict:
        for pos in feature_dict[seq_id]:
            for feature in feature_dict[seq_id][pos]:
                if feature not in feature_unique:
                    feature_unique[feature] = 0
                feature_unique[feature] +=1

    return list(feature_unique.keys()), feature_unique

def return_longest_feature(feature_list):
    """
    returns the longest feature that occurs in a list, according to the most number of : found

    Parameters
    ----------
    feature_list: list
        list of strings of features
    Returns
    -------
    longest_feature: str
        The longest feature found in feature_list, based on having concatenations with :
    """


    longest_feature = feature_list[0]


    #return the feature with the most : and that is the longest feature

    if len(feature_list)> 1:
        for feature in feature_list:
            if feature.count(':') > longest_feature.count(':'):
                longest_feature = feature
    return longest_feature

def return_feature_ann_dict(aln, feature_dict):
    """
    Given an aligned seqs object and a feature dictionary, this will return a dictionary of numpy dataframes with the number of features
    in each sequence, column position. You can look at this for the properties of a feature in a column

    Parameters
    ----------
    aln: Cogent3.load_aligned_seqs
        Cogent3 sequence alignment object, such as createdg by aln = load_aligned_seqs(alignment_file, moltype='protein')

    feature_dict : dict
        Dict in the form returned by return_features_from_file(file)
        Dict of dicts, with parent key the fasta header seq id, then inner dict with keys equal to the amino acid position
        and values equal to a list of features set for that feature.

    Returns
    -------
    feature_annotations: dict of dataframes
        keys are the unique features found in feature_dict and the values are dataframes based on numpy arrays of num_sequences x num_positions_alignment

    """
    feature_list, feature_num = return_unique_features(feature_dict)

    seq_translation_dict = {}
    for seq_id in feature_dict.keys():
        seq_translation_dict[seq_id] = aln.get_gapped_seq(seq_id).gap_maps()[0] #first argument is the key of the original sequence, value is the gapped position
    #remember this gap_map is 0-based counting, but feature_dict are 1-based counting

    num_positions = len(aln)
    num_sequences = len(aln.names)
    #instantiate the numpy arrays
    #let's make it a pandas dataframe
    feature_annotations = {}
    for feature in feature_list:
        df = pd.DataFrame(data=np.zeros([num_sequences, num_positions]), index = aln.names, columns=range(0,len(aln)))
        feature_annotations[feature] = df
    #now walk through the features and add to the position in the correct dataframe and position.
    for seq_id in feature_dict.keys():
        for pos in feature_dict[seq_id].keys():
            for feature_type in feature_dict[seq_id][pos]:
                try:
                    aln_pos = seq_translation_dict[seq_id][int(pos)-1]
                    feature_annotations[feature_type].loc[seq_id,aln_pos]=1
                except: #in case a feature exists outside the boundary of the alignment, which has happened on occasion
                    print("ERROR: cannot map %s to translation dict for %s"%(pos, seq_id))


    return feature_annotations

def print_ann_file(feature_annotations, feature_colors, annotation_file):
    """
    Given feature annotations, such as created by return_feature_ann_dict, print a Jalview Annotation file with bargraphs. 

    Parameters
    ----------
    feature_annotations: dict of dataframes
        keys are the unique features found in feature_dict and the values are dataframes based on numpy arrays of num_sequences x num_positions_alignment
    feature_colors: dict
        This is the color scheme read from the feature file (e.g. feature_colors, features = jf.return_features_from_file(feature_file) )
        so that the bar graph colors match the feature colors from the feature file
    annotation_file: str
        This is the name of the file to output the annotation tracks to, it will be overwritten

    Returns
    -------
    """

    af = open(annotation_file, "w")
    af.write("JALVIEW_ANNOTATION\n")
    #use the same colors as loaded from the feature_file
    for feature in feature_annotations:
        af.write("BAR_GRAPH\t%s\t%s\t"%(feature, feature))
        binValues = feature_annotations[feature].sum().values
        line = ''
        for i in range(0, len(binValues)):
            char = '-,'
            if binValues[i]>10:
                line += "%d,%s,%d|"%(int(binValues[i]), '*', int(binValues[i]))
            elif binValues[i]>0: 
                line += "%d,%s,%d|"%(int(binValues[i]), '*', int(binValues[i]))
            else:
                line += "%d,%s,|"%(int(binValues[i]), '-')
                #char='*,%d'%(int(binValues[i]))
            
        af.write(line)
        af.write("\n")
        af.write("COLOUR\t%s\t%s\n"%(feature,feature_colors[feature]))

    af.close()
    
def return_position_mapping_dict(alignment_file):
    """
    Given an alignment file (fasta format), return a dictionary with the mapping of the alignment position to the 
    original sequence position to the new positions. 

    Parameters
    ----------
    alignment_file: str
        Name of the alignment file in fasta format
    
    Returns
    -------
    seq_translation_dict: dict
        Dict with keys as the sequence ids and values a dictionary. Subdictionary has keys int positions (zero-based)
         of the protein sequence and values the alignment position (also int and zero-based)

    """
    aln = load_aligned_seqs(alignment_file,  moltype='protein')
    seq_translation_dict = {}
    for seq_id in aln.names:
        seq_translation_dict[seq_id] = aln.get_gapped_seq(seq_id).gap_maps()[0] #first argument is the key of the original sequence, value is the gapped position
    #remember this gap_map is 0-based counting, but feature_dict are 1-based counting
    return seq_translation_dict, aln

def return_intersection(feature_set_1, feature_set_2):
    """
    Given two feature set dictionaries, return the intersection of the second set with the first
    Returns a dictionary outer key the header, inner key the position, and values a list of the features that 
    directly overlapped. 
    Assumes that both feature sets are on the same positioning system (i.e. translated to the same reference sequences)

    Parameters
    ----------
    feature_set_1:
        dict in the form returned by return_features_from_file(file)
        Dict of dicts, with parent key the fasta header seq id, then inner dict with keys equal to the amino acid position
        and values equal to a list of features set for that feature.
    feature_set_2:
        dict in the form returned by return_features_from_file(file)
        Dict of dicts, with parent key the fasta header seq id, then inner dict with keys equal to the amino acid position
        and values equal to a list of features set for that feature.
    Returns
    -------
    intersection_dict:
        A dictionary containing the headers of feature_set_1 for sequences where at least one position overlapped 
        with feature_set_2. The inner dictionary key are the positions and these point to value arrays of 
        feautres in set 2 that occurred at the same position. 
    

    """

    headers_not_overlapping = []
    feature_intersection = {}
    for header in feature_set_1:
        feature_intersection[header] = {}
        for site in feature_set_1[header]:
            feature_intersection[header][site] = []

            if header in feature_set_2:
                if site in feature_set_2[header]:
                    feature_intersection[header][site] += feature_set_2[header][site]
                #print(header, site, domain_int_dict[header][site])
    else:
        headers_not_overlapping.append(header)
    return feature_intersection, headers_not_overlapping