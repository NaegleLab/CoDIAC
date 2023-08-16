import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from six.moves import urllib
from sys import exit
import re
from CoDIAC import contactMap as cm

#Here are helpers for orienting a ligand chain with a pTyr for analysis of SH2 binding pockets

def return_pos_of_interest(resNums, seq, pos, n_term_num=2, c_term_num=4, PTR_value = 'y'):
	"""
	Given the resNums and a transDict, return the resNums starting region and ending region to that captures a 
	pTyr ligand of interest. Include the new labels as well for easy heatmap plotting. Will return a list of tuples, so if
	there is more than one pTyr, then there is a set of regions returned. 
	
	Parameters
	----------
	resNums: list
		list of residue numbers that comes from chainMap class 
	seq: string
		sequence string to use that matches resNums for assembly of the tick_labels
	pos: int
		position of the residue of interest (i.e. a PTR)
	n_term_num: int
		number to the n-terminal side for flanking region
	c_term_num: 
		number of c-terminal flanking residues. 
	PTR_value: str
		char or string value to be printed as the xtick-label at position of interest. Default 'y'

	Returns
	-------

	"""

	minRes = resNums[0]
	maxRes = resNums[-1]
	num_pad_c_term = 0
	num_pad_n_term = 0
	if pos not in resNums:
		raise NameError("Transdict residue %d not in available resNums, be sure you are matching the correct object"%(pos))

	if len(PTR_value) > 1:
		raise NameError("PTR_value must be a single character")

	res_start = pos - n_term_num
	if res_start < minRes:
		num_pad_n_term = minRes - res_start
		res_start = minRes

	res_end = pos + c_term_num #have to offset for one's positioning and that we want the central position and n_term beyond
	if maxRes - res_end < 0:
		num_pad_c_term = res_end - maxRes
		res_end = maxRes

	#make the aligned string, first change the string to indicate the character change

	seq_changed = seq[0:pos-minRes] + PTR_value + seq[pos+1-minRes:maxRes]
	#print("DEBUG: seq_changed %s"%(seq_changed))

	aligned_str = ''
	for i in range(0, num_pad_n_term):
		aligned_str = aligned_str+'-'
	aligned_str+= seq_changed[res_start-minRes:res_end-minRes+1]
	for i in range(0, num_pad_c_term):
		aligned_str = aligned_str + '-'


	#to build the tick_labels, I want to combine normal resNums and sequence, except we'll rename the resNums around the central position
	resNums_copy = resNums.copy()
	#central position
	for i in range(res_start, res_end+1):
		value = i-pos
		ind = i-minRes
		resNums_copy[ind] = value 
	tick_labels = cm.return_aa_pos_list(seq_changed, resNums_copy)

	return res_start, res_end, aligned_str, tick_labels



