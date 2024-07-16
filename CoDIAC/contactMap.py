import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from six.moves import urllib
from sys import exit
import re
#import globals
#globals.initialize()
from CoDIAC.PDBHelper import PDBEntitiesClass

class globals:

  def initialize(): 
    global PTM_CONTACT_DICT 
    #to look up a ligand code replace the <> with intention in this URL
    # https://www.rcsb.org/ligand/<ALY>
    
    #PTM_CONTACT_DICT created using AdjacencyFiles.makePTM_dict
    PTM_CONTACT_DICT = {'PTR': 'Y',
                         'CAS': 'C',
                         'SEP': 'S',
                         '1PA': 'F',
                         'PM3': 'F',
                         'CME': 'C',
                         'FTY': 'Y',
                         'PTH': 'Y',
                         'CSO': 'C',
                         'ALY': 'K',
                         'MSE': 'M',
                         'SLZ': 'K',
                         'PTM': 'Y',
                         '1AC': 'A',
                         '02K': 'A',
                         'CCS': 'C',
                         'TPO': 'T'}
    #FYQ, AYI, YEN are synthetic pTYR sequences that act as inhibitors, have to decide what to do about that.
    global MOLECULES_TO_REMOVE
    MOLECULES_TO_REMOVE = ['HOH', 'SO4', 'NBS', 'CSO', 'MN', 'MG', 'ZN', 'MYR', 'P16', 'GOL', 'QUE', 'CA', 'ANP', 'EDO', 'DVT', 'CL', 'PO4', 'FMT', 'ACE', 'CAT', '1N1', 'VSH', 'PB']
    synthetic_inhibitors = ['FYQ', 'AYI', 'YEN', 'AYQ']
    MOLECULES_TO_REMOVE = MOLECULES_TO_REMOVE+synthetic_inhibitors


class chainMap:

  def __init__(self, PDB_ID, entity):
    """
    Returns an empty chainMap object for filling by translation or by read of a path file (use construct)

    """

    self.PDB_ID = PDB_ID
    self.entity = entity
    self.aaRes_dict = {}
    self.arr = []

    self.transDict = {}
    self.structSeq = ''
    self.resNums = []
    self.adjacencyDict = {}
    self.unmodeled_list = []
    globals.initialize()


  def construct(self, PATH):
    """
    Given a dataframe from a contact map, return an object with structure sequence extracted, unmodeled regions, modified amion acids, and the intrachain 
    contact map. 

    Parameters
    ----------
    PDB_ID: str
        The PDB_ID of interest, this assumes that the contactMap of this file exists in the target directory 
    PATH: string
      location of PDB contact map files (directory root where folders of PDB ID are under)
    chain: char
        chain to create information about 
    ang_cutoff: float
      Angstrom cutoff to use for considering something a contact (default <= 4Angstroms)
    noc: int
      Number of minimum contacts required
    Attributes
    ----------
    PDB_ID: str
      storing the PBD ID of this object
    chain: char
      stores the chain 
    adjacencyDict: dict
      dict of dict, with outer key being the residue and inner keys being the residue that outer key interacts with and value is number of contacts
    aaRes: list
      lists of the amino acids, concatenated with the residue number
    arr: array
      A zeros-based index adjancency matrix of contact numbers between all residues extracted from the file, square as it is chain-chain contacts
    structSeq: string
      the amino acid string extracted directly from the contact map, using - for missing residues as the value
    resNums: list
      The list of amino acids extracted directly from the contact map. 
    transDict: dict
      Dictionary of amino acid positions and the modification found (using globals.PTM_CONTACT_DICT to translate)
    unmodeled_dict: list
      resolved list of unmodeled residues with start, stop sets in the list

    """
    #self.chainMap(PDB_ID, chain)

    PDB_ID = self.PDB_ID
    entity = self.entity
    file = PATH+PDB_ID + '/' + PDB_ID + '_BF.txt'

    try: 
      df = pd.read_csv(file, sep='\t')
    except:
      exit(PDB_ID+" does not exist at " +file)

    for molecule in MOLECULES_TO_REMOVE:
      df = df.loc[lambda df: df['Res1'] != molecule]
      df = df.loc[lambda df: df['Res2'] != molecule]
    df = df[df['Res1'].notnull()]
    df = df[df['Res2'].notnull()]
   # df = df[df['Distance'] <= ang_cutoff]
   # df = df[df['Number of atomic contacts'] >= noc]

    self.structSeq, self.resNums, self.transDict, self.adjacencyDict = self.return_struct_sequence_info(df, entity)
    self.unmodeled_list = self.return_unmodeled_regions()
    try:
      self.arr = self.return_symmetric_single_chain_arr()
    except:
      print("error constructing matrix")

    #return the array

  

  def return_struct_sequence_info(self, df, entity):
    """
    This iterates through the dataframe, df, for a single chain to construct the structure sequence and intrachain array

    Returns
    -------
    structSeq: str
      The structure sequence, using '-' to indicate unmodeled locations
    resNums: list
      residue numbers extracted in sequential order including the missing residues
    transDict: dict
      Keys are the residue numbers, in structure sequence numbering, and value is the PTM (that was converted from)
    adjacencyDict: dict
      dict of dict, with outer key being the residue and inner keys being the residue that outer key interacts with and value is number of contacts  

    """
    df_single_chain = df[df['Entity1']==entity]
    df_single_chain = df_single_chain[df_single_chain['Entity2']==entity]
    
    aa_dict = {}
    adjList = {}
    transDict = {}
    for name, group in df_single_chain.groupby(['ResNum1', 'ResNum2']):
      res1_num = group.ResNum1.values[0]
      res1_aa = group.Res1.values[0]
      res2_num = group.ResNum2.values[0]
      res2_aa = group.Res2.values[0]
      if res1_num not in aa_dict:
        aa_dict[res1_num] = res1_aa
      if res2_num not in aa_dict:
          aa_dict[res2_num] = res2_aa
      innerDict = {}

      innerDict[res1_num] = len(group)
      if res2_num not in adjList:
        #if res2_num not in adjList:
        adjList[res2_num] = innerDict
      else:
        dict_temp = adjList[res2_num] #get a pointer to be able to append to the dictionary of dictionary
        if res1_num not in dict_temp:
          dict_temp[res1_num] = group.Binary_Feature.values[0]


    resNums_temp = list(aa_dict.keys())
    resNums_temp.sort()
    #print(resNums_temp)

    resNums = list(range(int(resNums_temp[0]), int(resNums_temp[-1])+1)) #not sure why, but the list with range ends one number prematurely
    #print(resNums)

    structSeq = '-'*len(resNums)
    structSeqList = list(structSeq)

    minRes = resNums[0]
    maxRes = resNums[-1]
    for res in resNums:
      if res in aa_dict:
        aa = aa_dict[res]
        if len(aa)>1:
            if aa in PTM_CONTACT_DICT:
                transDict[res] = aa
                aa = PTM_CONTACT_DICT[aa]
            else:
                print("ERROR: encountered new amino acid code %s"%(aa))
                aa = '-'
        pos = res - minRes        
        structSeqList[pos] = aa #update position, else it's left as -
          #print("updating %d with %s"%(pos, aa))
    structSeq = ''.join(structSeqList)
    return structSeq, resNums, transDict, adjList

  def return_min_residue(self):
    return self.resNums[0]
  def return_max_residue(self):
    return self.resNums[-1]

  def return_unmodeled_regions(self):
    """
    Uses the structure sequence with the - to indicate there was missing data
    unmodeled list is returned as 0-based counting 
    """
    unmodeled_list = []

    #add all missing contacts in the middle of the string
    result = [_.start() for _ in re.finditer('-', self.structSeq)]

    if result: 
      #walk through and find the first and end of a sequence
      start = result[0]
      end = result[0]
      i = 0
      while i < len(result)-1:
        i+=1
        #looking forward 1 to compare to current end. If it is a continuous number, then update end
        if result[i] == end + 1:
          end = result[i]
        else:
          unmodeled_list.append(start+self.resNums[0])
          unmodeled_list.append(end + self.resNums[0])
          start = result[i]
          end = result[i]
      #for the last one need to add the final start, end since we ran out of results
      unmodeled_list.append(start+self.resNums[0])
      unmodeled_list.append(end + self.resNums[0])

    return unmodeled_list

  def return_symmetric_single_chain_arr(self):
    """
    Given an adjacencyDict for a single chain, create the symmetric adjacecny matrix
    """
    arr = []

    
    minResidue =self.return_min_residue()
    maxResidue = self.return_max_residue()
    #print("minResidue is %d, maxresidue is %d"%(minResidue, maxResidue))
    size = maxResidue - minResidue+1
    arr = np.zeros((size, size))
    #print("array size is: "+str(size))
    #for i in range(minResidue+1, maxResidue):
    for i in self.adjacencyDict:
      adjacencyDict_inner = self.adjacencyDict[i]
      for j in adjacencyDict_inner:
        arr[i-minResidue][j-minResidue] = adjacencyDict_inner[j]
        arr[j-minResidue][i-minResidue] = adjacencyDict_inner[j]
  
      
    return arr


  def print_fasta_feature_files(self, featureStart, N_offset1,featureEnd, C_offset1,contactFromStart, N_offset2,contactFromEnd, C_offset2,fastaHeader, contactLabel, outputFileBase, threshold = 1, append = True, color = '117733', use_ref_seq_aligned=True):
    """
    Create a feature file for the ROI_1 that has contacts to ROI_2 of protein. 

    efined as ROI_1 and ROI_2 that are between (xStart, xEnd)
    and (yStart, yEnd)

    featureStart: int
        start position of ROI_1
    featureEnd: int
        end position of ROI_1
    contactFromStart: int
        start position of ROI_2
    contactFromEnd: int
        end position of ROI_2
    fastaHeader: str
        fastaHeader to be used to reference in jalview 
    contactLabel: str
        label of feature
    outputFileBase: str
        file name base to write feature and fasta to (it will be <outputFileBase>.fasta and <outpufileBase>.fea)
    threshold: float
        number of atomistic contacts to consider as a contact
    append: bool
        Whether to append or overwrite the feature file
    color: 
        color number for javliew to use when loading feature.
    use_ref_seq_aligned: bool
      True if you want to use the ref_seq_aligned where unmodeled regions have been imputed during translation, otherwise uses structseq

    """

    minRes = self.return_min_residue()  
    if use_ref_seq_aligned:
      seq = self.refseq
    else:
      seq = self.structSeq #right now there might be unmodeled regions 

    print_fasta_feature_files(self.arr, seq, featureStart, N_offset1, featureEnd,C_offset1, minRes, contactFromStart,N_offset2,contactFromEnd, C_offset2,minRes, fastaHeader, contactLabel, outputFileBase, threshold, append, color)
    

  def generateAnnotatedHeatMap(self, xStart, xEnd, yStart, yEnd, remove_no_contacts=True, text_annotate = 'on', use_ref_seq_aligned=True):
    '''
    Create a heatmap between two different regions of the protein for a single chain. Defined as ROI_1 and ROI_2 that are between (xStart, xEnd)
    and (yStart, yEnd). The heatmap will plot the first region as the rows and second as columns in the array.

    xStart: int
        start position of ROI_1, assumes it is in the crystal structure offset number
    xEnd: int
        end position of ROI_1, assumes it is in the crystal structure offset number
    yStart: int
        start position of ROI_2, assumes it is in the crystal structure offset number
    yEnd: int
        end position of ROI_1, assumes it is in the crystal structure offset number
    text_annotate: str
        'on' or any other string sets to not on
    use_ref_seq_aligned: bool
      True if you want to use the ref_seq_aligned where unmodeled regions have been imputed during translation, otherwise uses structseq

    '''

    offset = self.return_min_residue()

    if use_ref_seq_aligned:
      seqToUse = self.refseq
    else:
      seqToUse = self.structSeq

    rowTickLabels = return_aa_pos_list(seqToUse, self.resNums)
    generateAnnotatedHeatMap(self.arr, xStart, xEnd, offset, yStart, yEnd, offset, rowTickLabels, rowTickLabels, remove_no_contacts, text_annotate)


def return_single_chain_dict(PDB_descriptive_df, PDB_ID, PATH, entity):
  """
  This will create a dictionary of a single chain contact map and aligned contact map for a chain, entity  
  of a PDB_ID. Keys of this dictionary include:
    'chain'
    entity'
    'pdb_class', 'cm', 'cm_aligned'
  Parameters
  ----------
  PDB_description_df: pandas dataframe
    dataframe made by reading the descriptive file of PDBs and information
  PDB_ID: str
    PDB_ID
  PATH: str
    Path to where contact maps are
  chain: char
    chain 
  entity: int
    entity number - should match a chain.

  Returns
  -------
  cm_dict: dict
    Dictionary with controlled keys that assemble the contact maps (unaligned and aligned) for a chain and entity. Keys
    'chain' - chain (char)
    'entity' - entity (int)
    'pdb_class' - pdb_class object for the master file
    'cm' - chain map object
    'cm_aligned' - aligned chain map


  """
  cm_dict = {}
  #cm_dict['chain'] = chain
  cm_dict['entity'] = entity
  cm_dict['PDB_ID'] = PDB_ID

  entities = PDBEntitiesClass(PDB_descriptive_df, PDB_ID)
  if entity not in entities.pdb_dict:
    raise NameError("Entity %d not found in main file"%(entity))

  cm_dict['pdb_class'] = entities.pdb_dict[entity]

  #if chain not in entities.pdb_dict[entity].chain_list:
   # raise NameError("Chain %s not found for entity %d"%(chain, entity))

  cm_dict['cm'] = chainMap(PDB_ID, entity)
  cm_dict['cm'].construct(PATH)

  cm_dict['cm_aligned'] = translate_chainMap_to_RefSeq(cm_dict['cm'], cm_dict['pdb_class'])

  return cm_dict



def return_interChain_adj(PATH, from_dict, to_dict):
    """
    Given two single_chain_dict objects of two chains, return the adjList and the array 
    of contacts between chains. This uses the from_dict as the outer keys of adjList (row values of matrix)
    and to_dict as the inner keys of adjList (col values of matrix)
    
    Parameters
    ----------
    PATH: str
      Path to where contact maps are
    from_dict: dict
      The dict created by return_single_chain_dict that will serve as mapping FROM from_dict to to_dict
    to_dict: dict
      The dict created by return_single_chain_dict that will serve as mapping to_dict 
    ang_cutoff: float
      Angstrom cutoff to use for considering something a contact (default <= 4Angstroms)
    noc: int
      Number of minimum contacts required

    Returns
    -------
    adjDict: dict
      adjacency dict with outer keys as the from_dict residue numbers, inner keys the to_dict residue numbers
      and value equal to the number of contacts made between pair of residues (between chains) of ang_cutoff and noc
    arr: array
      Array with from_dict x to_dict size (numresidues) and values equal to the number of contacts from adjDict
    """
    
    PDB_ID = from_dict['PDB_ID']
    if PDB_ID != to_dict['PDB_ID']:
        raise NameError("The PDB_IDs of the two dictionaries are not the same")
    file = PATH+PDB_ID + '/' + PDB_ID + '_BF.txt'
    df = pd.read_csv(file, sep='\t')
    
 #   df = df[df['Distance'] <= ang_cutoff]
 #   df = df[df['Number of atomic contacts'] >= noc]
    df = df[df['Res1'].notnull()]
    df = df[df['Res2'].notnull()]
    
    
    chains = [from_dict['entity'], to_dict['entity']]
    df_a_to_b = df[(df['Entity1']==from_dict['entity']) & (df['Entity2']==to_dict['entity'])]
    df_b_to_a = df[(df['Entity2']==from_dict['entity']) & (df['Entity1']==to_dict['entity'])]
    # df_multichain = df_a_to_b.append(df_b_to_a)
    df_multichain = pd.concat([df_a_to_b, df_b_to_a])
    
    
    
    
    #this dict will allow me to access how to address the from versus to relationship by addressing directly 
    dict_combined = {}
    dict_combined[from_dict['entity']] = from_dict 
    dict_combined[to_dict['entity']] = to_dict
    adjList = {}
    for name, group in df_multichain.groupby(['ResNum1', 'ResNum2']):
        chain1 = group['Entity1'].iloc[0]
        chain2 = group['Entity2'].iloc[0]
        res1 = group['ResNum1'].iloc[0]
        res2 = group['ResNum2'].iloc[0]

        res1_offset = dict_combined[chain1]['cm_aligned'].offset
        res2_offset = dict_combined[chain2]['cm_aligned'].offset

        if res1-res1_offset in dict_combined[chain1]['cm_aligned'].resNums:
            if res2 - res2_offset in dict_combined[chain2]['cm_aligned'].resNums:

                #have to figure out if the key is the outer dict or inner dict key based on from/to
                if chain1 == from_dict['entity']:
                    outer_key = res1-res1_offset
                    inner_key = res2-res2_offset
                else: 
                    outer_key = res2-res2_offset
                    inner_key = res1-res1_offset


                innerDict = {}
                innerDict[inner_key] = len(group)

                if outer_key not in adjList:
            #if res2_num not in adjList:
                    adjList[outer_key] = innerDict
                else:
                    dict_temp = adjList[outer_key] #get a pointer to be able to append to the dictionary of dictionary
                    if inner_key not in dict_temp:
                        dict_temp[inner_key] = len(group)
                        
                        
        
    arr = np.zeros( [len(from_dict['cm_aligned'].resNums), len(to_dict['cm_aligned'].resNums)])
    
    from_start = from_dict['cm_aligned'].return_min_residue()
    to_start = to_dict['cm_aligned'].return_min_residue()
    for from_key in adjList.keys():
        from_index = from_key - from_start
        for to_key in adjList[from_key].keys():
            to_index = to_key-to_start
            arr[from_index, to_index] = adjList[from_key][to_key]
            
    return adjList, arr

def copy_chainMap(chain):
  """
  Returns a new copy of chain map

  """
  chain_copy = chainMap(chain.PDB_ID, chain.entity)
  chain_copy.arr = chain.arr.copy()
  chain_copy.transDict = chain.transDict.copy()
  chain_copy.structSeq = chain.structSeq
  chain_copy.resNums = chain.resNums.copy()
  chain_copy.adjacencyDict = chain.adjacencyDict.copy()
  chain_copy.unmodeled_list = chain.unmodeled_list.copy()
  return chain_copy



def translate_chainMap_to_RefSeq(entity, pdbClass):
  """
  Given a chainMap class, which includes the defined chain, use the pdbClass information to return
  a new chainMap class that is aligned to the reference sequence positions and full sequence expected.

  Parameters
  ----------
  chain: chainMap
    a chainMap created from a contact file
  pdbClass: pdbClass
    a pdbClass created from the annotation file for the matching PDB ID of chainMap

  Returns
  -------
  chainMap_aligned: chainMap
    A new chainmap that updates the residue numbers and array. This object appends the following things
    offset: int
      The offset required to map from refseq to the pdb contact map file if returning to that structure file
    ERROR_CODE: int
      0 if no errors, otherwise following codes are indicated. 1:no match found; 2:issue between struct and reference sequences after alignment; 3:struct and reference were different lengths


  """

  chainMap_aligned = copy_chainMap(entity)

  match, offset = find_offset(entity, pdbClass)
  n_term_gaps = 0
  c_term_gaps = 0
  c_term_excess = 0
  n_term_excass = 0
  structSeq_new = ''
  chainMap_aligned.match = match
  chainMap_aligned.offset = offset
  chainMap_aligned.ERROR_CODE = 0


  #cases to handle:
  # no match
  # match and we have to remove sequence at the end (struct Seq is > refseq)
  # it matches and we have to pad either the start or the end of the sequence and array

  if match:
    
    #here's where we translate old res nums to ref sequence resNums.
    newResNums = []
    for res in chainMap_aligned.resNums:
      newResNums.append(res-offset)
    chainMap_aligned.resNums = newResNums
    newUnmodeled_list = []
    for val in chainMap_aligned.unmodeled_list:
      newUnmodeled_list.append(val-offset)
    chainMap_aligned.unmodeled_list = newUnmodeled_list
    if chainMap_aligned.transDict:
      tempDict = {}
      for key in chainMap_aligned.transDict:
        tempDict[key-offset] = chainMap_aligned.transDict[key]
      chainMap_aligned.transDict = tempDict


    n_term_gaps = newResNums[0] - pdbClass.ref_seq_positions[0]
    n_term_excess = pdbClass.ref_seq_positions[0] - newResNums[0] #we have to cut stuff

    updateVals = 0
    if n_term_gaps > 0: # the contact map starts after pdbClass
      print("Adding %d n_term positions"%n_term_gaps)
      structSeq_new = pdbClass.ref_seq_mutated[0:n_term_gaps] + chainMap_aligned.structSeq
      # the refseq positions starting at the first and going forward need to be added too
      for i in range(0, n_term_gaps):
        newResNums.append(pdbClass.ref_seq_positions[0]+i)
      newResNums.sort()
      new_arr = np.insert(chainMap_aligned.arr, 0, np.zeros([n_term_gaps, len(chainMap_aligned.arr)]), axis=0)
      new_arr = np.insert(new_arr, 0, np.zeros([n_term_gaps, len(new_arr)]), axis=1)
      updateVals = 1
      chainMap_aligned.unmodeled_list.insert(0, newResNums[0]+n_term_gaps-1) #insert the end of the n-term gap and then insert the beginning
      chainMap_aligned.unmodeled_list.insert(0, newResNums[0])

      

    elif n_term_excess > 0:
      print("Deleting %d amino acids"%(n_term_excess))
      structSeq_new = chainMap_aligned.structSeq[n_term_excess:-1]
      newResNums = newResNums[n_term_excess:-1]
      new_arr = chainMap_aligned.arr
      for i in range(0, n_term_excess): #delete the number of rows and columns of excess from beginning of array
        new_arr = np.delete(new_arr, 0, axis=0)
        new_arr = np.delete(new_arr, 0, axis=1)

      updateVals = 1
      #Remove unmodeled from umodeled list -  is that best just to return the unmodeled from the new structSeq? 
     

    if updateVals:
      chainMap_aligned.structSeq = structSeq_new
      newResNums.sort()
      chainMap_aligned.resNums = newResNums
      #print("Debug: resnums after update")
      chainMap_aligned.arr = new_arr



    c_term_gaps = pdbClass.ref_seq_positions[-1] - chainMap_aligned.return_max_residue()
    #print("Debug: chainMap last residue currently is %d"%(chainMap_aligned.return_max_residue()))
    c_term_excess = -c_term_gaps

    if c_term_gaps > 0:
      # need to pad the c-terminal tail
      print("Adding %d c-terminal gaps"%c_term_gaps)
      len_ref = len(pdbClass.ref_seq_mutated)
      c_term_ext_str = pdbClass.ref_seq_mutated[len_ref-c_term_gaps:len_ref]
      structSeq_new = chainMap_aligned.structSeq + c_term_ext_str
      maxRes = chainMap_aligned.return_max_residue()
      for i in range(0, len(c_term_ext_str)):
        ind_to_add = maxRes+i+1
        #print("DEBUG: extending list of positions")
        chainMap_aligned.resNums.append(ind_to_add)
        #add a new column and row of zero values to the array
      pos_of_end_of_arr = len(chainMap_aligned.arr)
      new_arr = np.insert(chainMap_aligned.arr, pos_of_end_of_arr, np.zeros([c_term_gaps, len(chainMap_aligned.arr)]), axis=0)
      new_arr = np.insert(new_arr, pos_of_end_of_arr, np.zeros([c_term_gaps, len(new_arr)]), axis=1) 
      chainMap_aligned.arr = new_arr
      chainMap_aligned.structSeq = structSeq_new

      #add the new unmodeled region at the ned
      chainMap_aligned.unmodeled_list.append(maxRes)
      chainMap_aligned.unmodeled_list.append(maxRes+c_term_gaps)



    elif c_term_excess > 0: 
      #need to remove some of the c-terminal
      print("Removing an excess of %d"%(c_term_excess))
      #remove the last resNums, update the structure sequence, drop columns and rows from the array and update the unmodeled list.
      len_ref = len(pdbClass.ref_seq_mutated)
      chainMap_aligned.structSeq = chainMap_aligned.structSeq[0:len_ref]
      chainMap_aligned.resNums = chainMap_aligned.resNums[0:len_ref]
      chainMap_aligned.arr = chainMap_aligned.arr[0:len_ref, 0:len_ref]
      unmodeled_list_temp = []
      count = 0
      while count < len(chainMap_aligned.unmodeled_list)-1:
        #go through in pairs, so we can repeat the start or stop if needed 
        start = chainMap_aligned.unmodeled_list[count]
        stop = chainMap_aligned.unmodeled_list[count+1]
        count+=2
        if start in chainMap_aligned.resNums:
          unmodeled_list_temp.append(start)
          if stop in chainMap_aligned.resNums: 
            unmodeled_list_temp.append(stop)
          else: #we've hit the end of the protein and we should append the end
            unmodeled_list_temp.append(chainMap_aligned.return_max_residue())
        elif stop in chainMap_aligned.resNums: #the start of an unmodeled region is not in it, but the stop is
          unmodeled_list_temp.append(chainMap_aligned.return_min_residue())
          unmodeled_list_temp.append(stop)
      chainMap_aligned.unmodeled_list = unmodeled_list_temp

      #Remove unmodeled from umodeled list -  is that best just to return the unmodeled from the new structSeq? 

    #In adjacency dict, remove offset and keep only residues in list at the end
    adjacencyDict_temp = {}
    for res in chainMap_aligned.adjacencyDict:
      new_res = res - chainMap_aligned.offset

      if new_res in chainMap_aligned.resNums:
        adjacencyDict_temp[new_res] = {}
        innerDict = {}
        for res_2 in chainMap_aligned.adjacencyDict[res]:
          new_res2 = res_2 - chainMap_aligned.offset
          if new_res2 in chainMap_aligned.resNums:
            innerDict[new_res2] = chainMap_aligned.adjacencyDict[res][res_2]
        adjacencyDict_temp[new_res] = innerDict
    chainMap_aligned.adjacencyDict = adjacencyDict_temp

    #now add an object that is the refseq sequence for the residues of the structure. 
    #since sometimes we have seen mutations in refseq that were not denoted, we should simply replace unmodeled regions of structSeq with refSeq

    if len(pdbClass.ref_seq_mutated) == len(chainMap_aligned.structSeq):
      chainMap_aligned.refseq = build_string_by_replacing_unmodeled(chainMap_aligned.structSeq, pdbClass.ref_seq_mutated, chainMap_aligned.unmodeled_list, chainMap_aligned.return_min_residue(), chainMap_aligned.return_max_residue())

      
      if chainMap_aligned.refseq != pdbClass.ref_seq_mutated:
        print("ERROR: When replacing unmodled regions, found that struct sequence and ref_seq_mutated are not the same")
        print(pdbClass.ref_seq_mutated)
        print(chainMap_aligned.refseq)
        chainMap_aligned.ERROR_CODE = 2
    else:
      print("ERROR: refseq is not the same size as structSeq. Did not assign a refseq")
      chainMap_aligned.ERROR_CODE = 3

  else: #else match did not work
    print("FATAL ERROR, could not find match")
    chainMap_aligned.ERROR_CODE = 1

  return chainMap_aligned 

def build_string_by_replacing_unmodeled(structSeq, refSeq, unmodeled_list, struct_start, struct_end):
  """
  Given the structure sequence with unmodeled_list of unmodeled regions
  Build a string that uses refseq from the unmodeled list to replace '-'
  """

  if len(structSeq) != len(refSeq):
    print("ERROR: cannot use refseq to replace unmodeled list")
    return '-1'
  new_seq = ''
  #change the unmodeled_list into a list of positions that need to come from refseq
  unmodeled_pos = []
  ind = 0
  while ind < len(unmodeled_list)-1:
    start = unmodeled_list[ind]-struct_start
    stop = unmodeled_list[ind+1]-struct_start+1
    vals = list(range(start, stop))
    unmodeled_pos.extend(vals)
    ind+=2
 # print("DEBUG: unmodeled pos")
  # print(unmodeled_pos)
  for i in range(0, struct_end-struct_start+1):
    if i in unmodeled_pos:
      new_seq = new_seq + refSeq[i]
    else:
      new_seq  = new_seq  + structSeq[i]

  return new_seq



def find_offset(chainClass, pdbClass):
  """
  Given a structure sequence, which may have unmodeled residues indicated by '-'
  find a match within ref_seq of a string of struct_seq and identify the position where
  struct_seq reliably begins in ref_seq. This uses the longest continuous chain of modeled residues
  found in struct_seq s

  Parameters:
  -----------
  chainClass: chainClass 
    chainClass built from a contact map
  pdbClass: pdbClass
    pdbClass that matches the chain for chainClass

  Returns:
  --------
  match: bool
    True if match was found, false if not
  offset: int
    offset, where the start of the structure sequence was found in the reference_sequence
  """

  pos_start, length_str, sr = return_longest_string(chainClass.resNums, chainClass.unmodeled_list)
  start_in_structSeq = pos_start-chainClass.return_min_residue()+1
  substring = chainClass.structSeq[start_in_structSeq:start_in_structSeq+length_str-1]
  search = re.search(substring, pdbClass.ref_seq_mutated)
  match = False
  offset = 0
  if search:
    start, stop = search.span()
    ref_seq_pos = pdbClass.ref_seq_positions[0]+start#position in the substring + 
    print(ref_seq_pos)
    offset = pos_start +1 - ref_seq_pos 
    match = True

  else: # if it didn't work to look for a substring (e.g. b/c there was a fusion component), let's look for refseq in structure.
    totalLength = len(pdbClass.ref_seq_mutated)
    numToTry = 5
    strSize = min(int(totalLength/numToTry), 10)


    if(totalLength) < 20:
      numToTry = 3
      strSize = min(int(totalLength/numToTry), 3) #if the string is too small, like in peptides, we can shorten match size
      #print("DEBUG: shortening to 3")

    iterNum = 0
    search = re.search(pdbClass.ref_seq_mutated, chainClass.structSeq)
    subStart = 0
    #if a full match did not occur of ref_seq_mutated in structSeq (due to unomdeled residues, then we'll take different segments )
    while not search and iterNum<numToTry:

      subStart = iterNum*strSize
      #print("DEBUG: Trying indexes %d and %d"%(subStart, min(subStart+strSize, totalLength)))
      substring_val = pdbClass.ref_seq_mutated[subStart:subStart+min(subStart+strSize, totalLength)]
      search = re.search(substring_val, chainClass.structSeq)
      #print("DEBUG: trying a substring match with %s"%(substring_val))
      iterNum+=1
    if search:
      start, stop = search.span() #in this case, start is the string position where ref_seq matched the crystal structure.
      #print("DEBUG: start of refseq in crystal struct is %d, moved by %d of match"%(start, subStart))
      pos_in_ref_seq_start = pdbClass.ref_seq_positions[0] + subStart
      pos_in_struct_where_match = chainClass.resNums[0] + start
      offset = pos_in_struct_where_match - pos_in_ref_seq_start
      #print("Found offset")
      #therefore, that position in the crystal structure is the offset from ref_seq
      # if chainClass.resNums[start] > pdbClass.ref_seq_positions[0]:
      #   offset = pdbClass.ref_seq_positions[0] - chainClass.resNums[start] #offset is positive, have to subtract the offset
      # else:
      #   offset = chainClass.resNums[start] - pdbClass.ref_seq_positions[0]  

      match = True
  return(match, offset)




def return_longest_string(resNums, unmodeled_list):
  """
  Give the residue numbers and the unomdeled list, return the position (in resNum space)
  of the longest substring of a sequence.

  Parameters
  ----------
  resNums: list
    residue number list for understanding unmodeled list

  unmodeled_list: list
    paired list of unmodeled residues

  """
  stringSizes = {}
  if not unmodeled_list:
    stringSizes[resNums[0]] = len(resNums)
    pos_start_longest_str = resNums[0]
    length_longest_str  = len(resNums)
  else:
    stringSizes[resNums[0]] = unmodeled_list[0]-resNums[0]
    i = 1
    while i<len(unmodeled_list)-1:
        stringSizes[unmodeled_list[i]] = unmodeled_list[i+1]-unmodeled_list[i]
        i+=2
    stringSizes[unmodeled_list[-1]] = resNums[-1] - unmodeled_list[-1]
    pos_start_longest_str = max(stringSizes, key=stringSizes.get)
    length_longest_str = stringSizes[pos_start_longest_str]
  return pos_start_longest_str, length_longest_str, stringSizes

def print_fasta_feature_files(contact_arr, seq, featureStart, N_offset1, featureEnd, C_offset1, feature_minRes, contactFromStart, N_offset2, contactFromEnd, C_offset2,contact_minRes, fastaHeader, contactLabel, outputFileBase, threshold = 1, append = True, color = '117733', use_ref_seq_aligned=True):
  """
  Create a feature file for the ROI_1 (Region of Interest) that has contacts to ROI_2 of an array of interest. This assumes that the region of interest is represented in the rows of 
  the contact_arr.

  Parameters
  ----------
  contact_arr: array
    Has rows of the [featureStart, featureEnd] of interest and cols of the [contactFromStart, contactFromEnd] ROI. Cols can be the same or different species
  seq: str
    Protein sequence used to generate fasta file
  featureStart: int
      start position of ROI_1, assumes it is in the crystal structure minRes number
  featureEnd: int
      end position of ROI_1, assumes it is in the crystal structure minRes number
  feature_minRes: int
    minimum residue of the feature 
  contactFromStart: int
      start position of ROI_2, assumes it is in the crystal structure minRes number
  contactFromEnd: int
      end position of ROI_1, assumes it is in the crystal structure minRes number
  contact_minRes: int
    minimum residue of the contact component
  fastaHeader: str
      fastaHeader to be used to reference in jalview 
  contactLabel: str
      label of feature
  outputFileBase: str
      file name base to write feature and fasta to (it will be <outputFileBase>.fasta and <outpufileBase>.fea)
  threshold: float
      number of atomistic contacts to consider as a contact
  append: bool
      Whether to append or overwrite the feature file
  color: 
      color number for javliew to use when loading feature.
  use_ref_seq_aligned: bool
    True if you want to use the ref_seq_aligned where unmodeled regions have been imputed during translation, otherwise uses structseq

  """

  feature_file_name = outputFileBase+'.fea'
  fasta_file_name = outputFileBase+'.fasta'
  if not append:
      feature_file = open(feature_file_name, "w")
      feature_file.close()

      fasta_file = open(fasta_file_name, "w")
      fasta_file.close()

  # account for offset by subtracting featureStart by it; also account for python beginning at index 0 by subtracting 1
  featureStart = featureStart-feature_minRes
  featureEnd = featureEnd-feature_minRes
  contactFromStart = contactFromStart-contact_minRes
  contactFromEnd = contactFromEnd-contact_minRes
  
  rows, cols = contact_arr.shape
  
  if featureStart < 0:
      print("Invalid value for featureStart: value reset to 0")
      featureStart = 0
  elif contactFromStart < 0:
      print("Invalid value for contactFromStart: value reset to 0")
      contactFromStart = 0
  elif featureEnd > rows:
      print("Invalid value for featureEnd: value reset to max array index")
      featureEnd = rows
  elif contactFromEnd > cols:
      print("Invalid value for contactFromEnd: value reset to max array index")
      contactFromEnd = cols
  else:
      pass
  
  if append:
      feature_file = open(feature_file_name,"a")
  else:
      feature_file = open(feature_file_name,"a")
      feature_file.write("%s\t %s \n"%(contactLabel, color)) #this really needs to be added to the top of the file

  
  arr = contact_arr[featureStart+N_offset1:featureEnd+1+C_offset1, contactFromStart+N_offset2:contactFromEnd+1+C_offset2]
  arr = arr == threshold
  arr_sum = arr.sum(axis=1)
  featureLength = featureEnd +N_offset1- featureStart +1+C_offset1
  # print(featureLength,len(arr_sum))
  #contactLength = contactFromEnd - contactFromStart 
  for i in range(featureLength):
    if arr_sum[i]: #if we wanted, could require arr_sum[i] >= NumResiduesContacted
      feature_file.write("%s\t%s\t-1\t%s\t%s\t%s\n"%(contactLabel, fastaHeader, str(i+1), str(i+1), contactLabel))
      # print(i)
      #for j in range(contactLength):
          #if arr[i][j] > threshold:
           #   feature_file.write("%s\t%s\t-1\t%s\t%s\t%s\n"%(contactLabel, fastaHeader, str(i+1), str(i+1), contactLabel))
            #  break
  feature_file.close() 

  #write the fasta file
  fasta_file = open(fasta_file_name, "a")
  fasta_file.write('>'+fastaHeader+'\n')
  #if use_ref_seq_aligned:
  #  seq = self.refseq[featureStart:featureEnd]
  #else:
  seq_sub = seq[featureStart+N_offset1:featureEnd+1+C_offset1] #right now there might be unmodeled regions 
  fasta_file.write(seq_sub+'\n')
  fasta_file.close()


def return_aa_pos_list(seq, resNums):
  """
  Given a sequence and list of residue numbers, create an array of strings with the aminoacid and residue numbers concatenated
  """
  list_AA = []
  if len(resNums) != len(seq):
    print("ERROR: resNums and struct Seq do not match in length")
  else:
    for i in range(0, len(resNums)):
      list_AA.append(seq[i] + str(resNums[i]))
  return list_AA
 # print(list_AA)  

def generateAnnotatedHeatMap(contact_arr, rowStart, rowEnd, rowMinRes, colStart, colEnd, colMinRes, rowTickLabels, colTickLabels, remove_no_contacts=True, text_annotate = 'on'):
  '''
  Create a heatmap of arr for a subset of the array from ROI_1 defined by [rowStart,rowEnd] and columns from [colStart, colEnd]
  using the minimum residues to adjust for moving from reference sequence into the array

  rowStart: int
      start position of ROI_1 (the "from" region)
  rowEnd: int
      end position of ROI_1
  rowMinRes: int
    minimum residue number of the sequence in rows
  colStart: int
      start position of ROI_2 (the "to" region)
  colEnd: int
      end position of ROI_2
  colMinRes: int
    minimum residue number of the sequence in cols
  rowTickLabels: list
    list of strings that will be used as labels for rows
  colTickLabels: list
    list of strings that will be used for column labels 
  remove_no_contacts: bool
    Remove from the printed array any rows and columns that have no contacts.
  text_annotate: str
      'on' or any other string sets to not on
 
  '''

  # xStart = colStart
  # xEnd = colEnd
  # x_offset = colMinRes

  # yStart = rowStart
  # yEnd = rowEnd
  # y_offset = rowMinRes

  # #account for the offset
  # xStart = xStart-x_offset
  # xEnd = xEnd-x_offset
  # yStart = yStart-y_offset
  # yEnd = yEnd-y_offset

  # if xStart < 0:
  #   print("ERROR: xStart, when considering offset, is less than 0, setting new parameter")
  #   xStart = 0

  # if yStart < 0:
  #   print("ERROR: yStart, when considering offset, is less than 0, setting to 0")
  

  # rows, cols = contact_arr.shape
  # if xEnd > cols:
  #   print("ERROR: xEnd is past the array, setting to end")
  #   xEnd = cols
  # if yEnd > rows:
  #   print("ERROR: yEnd is past the array, setting to end")
  #   yEnd = rows

  # x_AA = colTickLabels[xStart:xEnd+1]
  # y_AA = rowTickLabels[yStart:yEnd+1]
  
  # arr = contact_arr[yStart:yEnd+1, xStart:xEnd+1]

  arr, x_AA, y_AA = return_arr_subset_by_ROI(contact_arr, rowStart, rowEnd, rowMinRes, colStart, colEnd, colMinRes, rowTickLabels, colTickLabels)

  if remove_no_contacts:
    row_vals_to_keep = np.where(arr.any(axis=1))[0]
    col_vals_to_keep = np.where(arr.any(axis=0))[0]
    subset = np.ix_(row_vals_to_keep, col_vals_to_keep)
    arr = arr[subset]
    x_AA_new = []
    y_AA_new = []

    for val in col_vals_to_keep:
      x_AA_new.append(x_AA[val])
    for val in row_vals_to_keep:
      y_AA_new.append(y_AA[val])
    x_AA = x_AA_new
    y_AA = y_AA_new


  #print("DEBUG: size of array to print is")
  #print(arr.shape)
  #print("DEBUG: supposed to be %d x %d"%(yEnd-yStart, xEnd-xStart))
  #print("DEBUG: object starting array is")
  #print(self.arr.shape)

  
  fig, ax = plt.subplots(figsize = (20,10))
  im = ax.imshow(arr, cmap='YlGn')

  # We want to show all ticks...
  ax.set_xticks(np.arange(len(x_AA)))
  ax.set_yticks(np.arange(len(y_AA)))
  # ... and label them with the respective list entries
  ax.set_xticklabels(x_AA)
  ax.set_yticklabels(y_AA)

  # Rotate the tick labels and set their alignment.
  plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
          rotation_mode="anchor")

  # Loop over data dimensions and create text annotations.
  if text_annotate == 'on':
    for i in range(len(y_AA)):
        for j in range(len(x_AA)):
            if arr[i,j] == 0:
                pass
            else:
                text = ax.text(j, i, arr[i, j], ha="center", va="center", color="k")

  ax.set_title("Annotated Contact Map")
  fig.tight_layout()
  plt.show() 


def return_arr_subset_by_ROI(contact_arr, rowStart, rowEnd, rowMinRes, colStart, colEnd, colMinRes, rowTickLabels=[], colTickLabels=[]):
  '''
  Return a subset of the array from ROI_1 defined by [rowStart,rowEnd] and columns from [colStart, colEnd]
  using the minimum residues to adjust for moving from reference sequence into the array

  rowStart: int
      start position of ROI_1 (the "from" region)
  rowEnd: int
      end position of ROI_1
  rowMinRes: int
    minimum residue number of the sequence in rows
  colStart: int
      start position of ROI_2 (the "to" region)
  colEnd: int
      end position of ROI_2
  colMinRes: int
    minimum residue number of the sequence in cols
  rowTickLabels: list
    list of strings that will be used as labels for rows
  colTickLabels: list
    list of strings that will be used for column labels 

  '''
  xStart = colStart
  xEnd = colEnd
  x_offset = colMinRes

  yStart = rowStart
  yEnd = rowEnd
  y_offset = rowMinRes

  #account for the offset
  xStart = xStart-x_offset
  xEnd = xEnd-x_offset
  yStart = yStart-y_offset
  yEnd = yEnd-y_offset

  if xStart < 0:
    print("ERROR: xStart, when considering offset, is less than 0, setting new parameter")
    xStart = 0

  if yStart < 0:
    print("ERROR: yStart, when considering offset, is less than 0, setting to 0")
  

  rows, cols = contact_arr.shape
  if xEnd > cols:
    print("ERROR: xEnd is past the array, setting to end")
    xEnd = cols
  if yEnd > rows:
    print("ERROR: yEnd is past the array, setting to end")
    yEnd = rows

  x_AA = []
  y_AA = []
  if len(rowTickLabels):
  
    y_AA = rowTickLabels[yStart:yEnd+1]
  if len(colTickLabels):
    x_AA = colTickLabels[xStart:xEnd+1]

  
  arr = contact_arr[yStart:yEnd+1, xStart:xEnd+1]

  return arr, x_AA, y_AA







