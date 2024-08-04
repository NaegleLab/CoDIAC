import os
import requests
import numpy as np
import pandas as pd
import logging
import ast
import json
import urllib.request

from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

AA_dict = {"ALA":'A',
                  "ARG":'R',
                  "ASN":'N',
                  "ASP":'D',
                  "CYS":'C',
                  "GLU":'E',
                  "GLN":'Q',
                  "GLY":'G',
                  "HIS":'H',
                  "ILE":'I',
                  "LEU":'L',
                  "LYS":'K',
                  "MET":'M',
                  "PHE":'F',
                  "PRO":'P',
                  "SER":'S',
                  "THR":'T',
                  "TRP":'W',
                  "TYR":'Y',
                  "VAL":'V'}

def ProcessArpeggio(input_file_path, outfile_path, mmCIF_file, small_molecule = False):
        ''' This function generates a text file of contacts that is a processed version of contact information extracted using Arpeggio (json file) 
        
        Parameters
        ----------
            input_file_path : str
                Input a .json file generated by running Arpeggio
            outfile_path : str 
                Path to save the output .txt contactmap file
            mmCIF_file : str
                Path to fetch the mmCIF file
            small_molecule : bool
                Whether to include or exclude protein interactions with any small molecules that may be present in the structure complex
                    
        Returns
        -------
            Adjacency file '{PDB_ID}.txt' to the specified path in outfile variable'''

       # filepath = self.PATH + '/' + self.input_json
        with open(input_file_path,'r') as file:
            data = json.load(file)
            
        pdb_info = MMCIF2Dict(mmCIF_file)
        
        parse_cifpath = mmCIF_file.split('/')
        for i in (parse_cifpath):
            if '.cif' in i:
                CIF_PDB_ID, ext = i.split('.')
                print(CIF_PDB_ID)
                
        if small_molecule == False:

            HETATMS_dict = {}

            if '_pdbx_nonpoly_scheme.pdb_mon_id' in pdb_info.keys():
                het_res = pdb_info['_pdbx_nonpoly_scheme.pdb_mon_id']
                het_res_id = pdb_info['_pdbx_nonpoly_scheme.pdb_seq_num']

                for key in het_res_id:
                    for value in het_res:
                        HETATMS_dict[int(key)] = value
                        het_res.remove(value)
                        break

        if small_molecule == True:
            HETATMS_dict = {}
            
        contact_types=['aromatic','carbonyl','hbond','hydrophobic','ionic','polar',
                       'vdw','vdw_clash','weak_hbond','weak_polar','xbond']
        
        parse_path = input_file_path.split('/')
        for i in (parse_path):
            if '.json' in i:
                PDB_ID, ext = i.split('.')
                print(PDB_ID)
                
        outfile = outfile_path+'/'+PDB_ID+'.txt'

        if CIF_PDB_ID == PDB_ID :

            with open(outfile,'w') as f:
                f.write('PDB'+'\t'+'Chain1'+'\t'+'Chain2'+'\t'+'Res1'+'\t'+'ResNum1'+'\t'+
                       'Res2'+'\t'+'ResNum2'+'\t'+'Atoms'+'\t'+'Distance'+'\t'+'Contact_type'+'\n')

                for i in range(len(data)):
                    chain_1 = ((data[i]['bgn']['auth_asym_id']))
                    atom_1 =((data[i]['bgn']['auth_atom_id']))
                    resnum_1 = (int(data[i]['bgn']['auth_seq_id']))
                    res_1 = ((data[i]['bgn']['label_comp_id']))
                    chain_2 = ((data[i]['end']['auth_asym_id']))
                    atom_2 = ((data[i]['end']['auth_atom_id']))
                    resnum_2 = (int(data[i]['end']['auth_seq_id']))
                    res_2 = ((data[i]['end']['label_comp_id']))
                    distance = float((data[i]['distance']))
                    interaction = (data[i]['interacting_entities'])
                    contactlist = data[i]['contact'] #stores data in a list 

                    if res_1 not in HETATMS_dict.values():               
                        if res_2 not in HETATMS_dict.values():                
                            if interaction == 'INTRA_SELECTION':
                                    if distance < 5:

                                        for contact in contactlist:

                                            if contact in contact_types:
                                                if resnum_1 > resnum_2:
                                                    atom_pair = atom_1 +'-'+atom_2

                                                    f.write(str(PDB_ID)+'\t'+str(chain_1)+'\t'+str(chain_2)+'\t'+
                                                            str(res_1)+'\t'+str(resnum_1)+'\t'+str(res_2)+'\t'+str(resnum_2)+'\t'+
                                                           str(atom_pair)+'\t'+str(distance)+'\t'+str(contactlist)+'\n')

                                                else:
                                                    atom_pair = atom_2 +'-'+atom_1
                                                    f.write(str(PDB_ID)+'\t'+str(chain_2)+'\t'+str(chain_1)+'\t'+
                                                            str(res_2)+'\t'+str(resnum_2)+'\t'+str(res_1)+'\t'+str(resnum_1)+'\t'+
                                                           str(atom_pair)+'\t'+str(distance)+'\t'+str(contactlist)+'\n')

                                                break

        df = pd.read_csv(outfile,  sep='\t')

        df.replace({"Res1": AA_dict},inplace=True)
        df.replace({"Res2": AA_dict},inplace=True)

        df.sort_values(["Chain1","Chain2","ResNum1","ResNum2"],
                            axis=0,
                            ascending=[True, True, True, True],
                            inplace=True)

        df.to_csv(outfile,sep='\t',encoding='utf-8',index=False)
        print('Adjacency File generated for {PDB}'.format(PDB = PDB_ID))
        

def BinaryFeatures(PDB_ID, PATH, translate_resid=False):
    ''' Produces an adjacency file for intra (domain-domain) and inter (liagnd) protein contacts that are represented as binary features. Contacts merged across and within protein entities and chains using defined thresholds. 
    
    Parameters
    ----------
        PDB_ID : str
            Input a PDB ID
        PATH : str
            Path to write the output file 
        translate_resid : bool
            chose whether to update structure sequence residue positions to the ones in the Uniprot Reference
            
    Returns
    -------
        OUTFILE : str
            binarized version of the adjacency txt file is created
    '''
    
    OUTFILE = PATH+'/'+PDB_ID+'_BF.txt'
    interpro_Adj = Interprotein_AdjFile(PDB_ID, PATH, translate_resid)
    intrapro_Adj = Intraprotein_AdjFile(PDB_ID, PATH, translate_resid)
    df = pd.concat([intrapro_Adj, interpro_Adj]) 
    df.to_csv(OUTFILE, index=None, sep='\t', mode='w') 
    print('Adjacency File with Binary features generated for', PDB_ID)
    
    
def Intraprotein_AdjFile(PDB_ID, PATH, translate_resid=False):
    ''' Produces adjacency file for only Intraprotein contacts along with binary features
    
    Parameters
    ----------
        PDB_ID : str
            Input a PDB ID
        PATH : str
            Path to find the arpeggio generated adjacency file and produce dataframe with binarized features for intraprotein 
        translate_resid : bool
            whether to update structure sequence residue positions to the ones in the Uniprot Reference

    '''
    
    INPUTFILE = PATH+'/'+PDB_ID+'.txt'
    #OUTFILE = PATH+'/'+PDB_ID+'/'+PDB_ID+'_intra.txt'
    
    df = pd.read_csv(INPUTFILE, sep='\t')
    
    #duplicate columns for residue numbers
    df['ResNum1_upd'] = df.loc[:, 'ResNum1']
    df['ResNum2_upd'] = df.loc[:, 'ResNum2']
    
    #join chain and residue number column values to use as identifier for updating the residue numbers if required
    df['map1'] = df['Chain1'].astype(str)+(df['ResNum1']).astype(str)
    df['map2'] = df['Chain2'].astype(str)+(df['ResNum2']).astype(str)
        
    #update residue numbers if needed
    if translate_resid == True:
        if ResID_from_DBREF(PDB_ID)[1].values():
            dict_resnums = (ResID_from_DBREF(PDB_ID)[1])
    if translate_resid == False:
        if ResID_from_DBREF(PDB_ID)[2].values():
            dict_resnums = (ResID_from_DBREF(PDB_ID)[2])
         
    df['ResNum1_upd'] = df['map1'].map(dict_resnums)
    df['ResNum2_upd'] = df['map2'].map(dict_resnums)
    
    #handle residues that are part of SEQADV record in PDB header (expression tags, variants, conflicts, etc.)
    for i in range(len(df)):
        if str(df.loc[i,'ResNum2_upd']) == 'nan':
            df.loc[i,'ResNum2_upd'] = df.loc[i, 'ResNum2']
        if str(df.loc[i,'ResNum1_upd']) == 'nan':
            df.loc[i,'ResNum1_upd'] = df.loc[i, 'ResNum1']
        
    #pair the chain and residue number for each row (representing a feature)
    df['ResNum pair'] = df['ResNum1_upd'].astype(str)+'-'+(df['ResNum2_upd']).astype(str)
    df['Chain Pair'] = df['Chain1'].astype(str)+'-'+(df['Chain2']).astype(str)

    entity_dict = entity_for_chain(PDB_ID)[0]
    chain_count = len(set(entity_for_chain(PDB_ID)[1]))
    unique_entities = set(entity_dict.values())
    #map entity values based on chain and add columns 
    df['Entity1'] = df['Chain1'].map(entity_dict)
    df['Entity2'] = df['Chain2'].map(entity_dict)

    #for intraprotein we filter features within same chain and repeat this for all entities present
    df2 = df[df['Chain1']==df['Chain2']]
    df3 = df2[df2['Entity1']==df2['Entity2']]
    df3 = df3.reset_index(drop=True)
    df3.insert(11,'Binary_Feature','')

    for entity in unique_entities:
        threshold = intraprotein_threshold(entity_dict, entity)

        for i in range(len(df3)):
            check_entity = df['Entity1'][i]
            respair = df3['ResNum pair'][i]
            if check_entity == entity:

                #group chain pairs based on the respair to find the occurrence of the feature across all chains
                fea_in_chains = df3.loc[df3['ResNum pair'] == respair, 'Chain Pair'].values.tolist()

                #modify threshold if needed
#                 if chain_count == 2: #for two chains separate condition
#                     if len(set(fea_in_chains)) > threshold:
#                         BF = 1
#                     else:
#                         BF = 0
#                 else:

                if len(set(fea_in_chains)) >= threshold:
                    BF = 1
                else:
                    BF = 0

                df3.at[i,'Binary_Feature'] = BF

    df_edit = df3.drop(['ResNum1', 'ResNum2', 'ResNum pair', 'Chain Pair', 'Atoms', 'Contact_type', 'Distance','map1','map2'], axis=1)
    df3 = df_edit.drop_duplicates()
    df3 = df3.drop(['Chain1','Chain2'],axis=1)
    df4 = df3.drop_duplicates()  
    df4 = df4.rename({'ResNum1_upd': 'ResNum1', 'ResNum2_upd': 'ResNum2'}, axis=1)
    new_cols = ['PDB', 'ResNum1', 'Res1', 'Entity1', 'ResNum2', 'Res2', 'Entity2', 'Binary_Feature']
    df4 = df4.reindex(columns=new_cols)
    df4 = df4.reset_index(drop=True) 
    df4.sort_values(["ResNum1","ResNum2"], 
                                    axis=0,
                                    ascending=[True, True], 
                                    inplace=True)
    df4 = df4.drop_duplicates()
    df4 = df4.reset_index(drop=True)
    df4['ResNum1'] = df4['ResNum1'].astype('int')
    df4['ResNum2'] = df4['ResNum2'].astype('int')
    df4['Binary_Feature'] = df4['Binary_Feature'].astype('int')
    #df4.to_csv(OUTFILE, index=None, sep='\t', mode='w')
    return(df4)

def Interprotein_AdjFile(PDB_ID,PATH, translate_resid=False):
    ''' Produces adjacency file for only Interprotein contacts along with binary features
    Parameters
    ----------
        PDB_ID : str
            Input a PDB ID
        PATH : str
            Path to find the arpeggio generated adjacency file and produce dataframe with binarized features for interprotein 
        translate_resid : bool
            whether to update structure sequence residue positions to the ones in the Uniprot Reference '''
    
    INPUTFILE = PATH+'/'+PDB_ID+'.txt'
    #OUTFILE = PATH+'/'+PDB_ID+'/'+PDB_ID+'_inter.txt'

    df = pd.read_csv(INPUTFILE, sep='\t')

    #duplicate columns for residue numbers
    df['ResNum1_upd'] = df.loc[:, 'ResNum1']
    df['ResNum2_upd'] = df.loc[:, 'ResNum2']
    
    #join chain and residue number column values to use as identifier for updating the residue numbers if required
    df['map1'] = df['Chain1'].astype(str)+(df['ResNum1']).astype(str)
    df['map2'] = df['Chain2'].astype(str)+(df['ResNum2']).astype(str)
    
    #update residue numbers if needed
    if translate_resid == True:
        if ResID_from_DBREF(PDB_ID)[1].values():
            dict_resnums = (ResID_from_DBREF(PDB_ID)[1])
    if translate_resid == False:
        if ResID_from_DBREF(PDB_ID)[2].values():
            dict_resnums = (ResID_from_DBREF(PDB_ID)[2])
         
    df['ResNum1_upd'] = df['map1'].map(dict_resnums)
    df['ResNum2_upd'] = df['map2'].map(dict_resnums)

    #handle residues that are part of SEQADV record in PDB header (expression tags, variants, conflicts, etc.)
    for i in range(len(df)):
        if str(df.loc[i,'ResNum2_upd']) == 'nan':
            df.loc[i,'ResNum2_upd'] = df.loc[i, 'ResNum2']
        if str(df.loc[i,'ResNum1_upd']) == 'nan':
            df.loc[i,'ResNum1_upd'] = df.loc[i, 'ResNum1']
    
    #pair the chain and residue number for each row (representing a feature)
    df['ResNum Pair'] = df['ResNum1_upd'].astype(str)+'-'+(df['ResNum2_upd']).astype(str)
    df['Chain Pair'] = df['Chain1'].astype(str)+'-'+(df['Chain2']).astype(str)

    entity_dict = entity_for_chain(PDB_ID)[0]
    threshold = interprotein_threshold(PDB_ID)
    
    #map entity values based on chain and add columns 
    df['Entity1'] = df['Chain1'].map(entity_dict)
    df['Entity2'] = df['Chain2'].map(entity_dict)
    df['Entity Pair'] = df['Entity1'].astype(str)+'-'+(df['Entity2']).astype(str)
    
    #for interprotein we filter features that are present between entities
    df2 = df[df['Entity1']!=df['Entity2']]
    df2 = df2.reset_index(drop=True)
    df2.insert(11,'Binary_Feature','')
    
    for i in range(len(df2)):
        respair = df2['ResNum Pair'][i] 
        
        #group chain pairs based on the respair to find the occurrence of the feature across all chains
        fea_in_chains = df2.loc[df2['ResNum Pair'] == respair, 'Chain Pair'].values.tolist()

        #modify threshold if needed
        if len(set(fea_in_chains)) >= threshold:
            BF = 1
        else:
            BF = 0
        df2.at[i,'Binary_Feature'] = BF
    
    df_edit = df2.drop(['ResNum1', 'ResNum2', 'ResNum Pair', 'Chain Pair','Entity Pair', 'Distance','Atoms','Contact_type','map1', 'map2'], axis=1)
    df3 = df_edit.drop_duplicates()
    df3 = df3.drop(['Chain1','Chain2'],axis=1)
    df4 = df3.drop_duplicates()
    df4.reset_index(drop=True)   
    df4[['PDB', 'ResNum1_upd', 'Res1', 'Entity1', 'ResNum2_upd', 'Res2', 'Entity2', 'Binary_Feature']]
    df4 = df4.rename({'ResNum1_upd': 'ResNum1', 'ResNum2_upd': 'ResNum2'}, axis=1)
    new_cols = ['PDB', 'ResNum1', 'Res1', 'Entity1', 'ResNum2', 'Res2', 'Entity2', 'Binary_Feature']
    df4 = df4.reindex(columns=new_cols)
    df4 = df4.reset_index(drop=True)
    
    #Flip the Residue number, residue name and entity values to align the features sequentially in the file
    for i in range(len(df4)):
        ent1 = int(df4['Entity1'][i])
        ent2 = df4['Entity2'][i]
        rn1 = df4['ResNum1'][i]
        rn2 = df4['ResNum2'][i]
        res1 = df4['Res1'][i]
        res2 = df4['Res2'][i]
        if ent1 > ent2:
            df4.at[i,'ResNum1'] = rn2
            df4.at[i,'ResNum2'] = rn1
            df4.at[i,'Res1'] = res2
            df4.at[i,'Res2'] = res1
            df4.at[i,'Entity1'] = ent2
            df4.at[i, 'Entity2'] = ent1
    
    df4.sort_values(["ResNum1","ResNum2"], 
                                    axis=0,
                                    ascending=[True, True], 
                                    inplace=True)
    df4 = df4.drop_duplicates()
    df4 = df4.reset_index(drop=True)
    df4['ResNum1'] = df4['ResNum1'].astype('int')
    df4['ResNum2'] = df4['ResNum2'].astype('int')
    df4['Binary_Feature'] = df4['Binary_Feature'].astype('int')
    #df4.to_csv(OUTFILE,index=None, sep='\t', mode='w')
    return(df4)
    
    
def ResID_from_DBREF(PDB_ID):
    ''' This function uses the DBREF record from PDB header to find differences in residue numbers between PDB and reference (such as Uniprot) If PDB legacy formats are unavailable, we make use of the mmCif headers and extract the same information 
    
    Parameters
    ----------
        PDB_ID : str
            Input a PDB ID
        
    Returns
    -------
        dbref_dict : dict
            Dictionary with chains as keys and list of sequence start and end residue numbers for both databases as values
        dict_upd_id : dict
            Dictionary whose keys are chain and residue number of the residues and its values are updated residue number. 
            Upon finding differences in residue numbers, the PDB sequence numbers will be replaced by the Uniprot residue number values stored in this dict 
        dict_id : dict
            Dictionary whose keys are chain and residue number of the residues and its values are residue number in the structure sequence(PDB). 
            '''
    
    dbref_dict={}

    res = requests.get("https://files.rcsb.org/header/"+ \
                            PDB_ID + ".pdb")
    if res.status_code != 404:
        response = urllib.request.urlopen("https://files.rcsb.org/header/"+ \
                        PDB_ID + ".pdb")
        for line in response:
            list_attributes = []
            line = (line.strip())
            if (line.startswith(b'DBREF ')):
                if isinstance(line, bytes):
                    data =(line.decode())
                    chain = data[12]
                    begin_pdb = int(data[14:18])
                    end_pdb = int(data[20:24])
                    begin = int(data[55:60])
                    end = int(data[62:67])
                    list_attributes.append(begin_pdb)
                    list_attributes.append(end_pdb)
                    list_attributes.append(begin)
                    list_attributes.append(end)
                    dbref_dict[chain] = list_attributes
    else:
        i = PDBList()
        file = i.retrieve_pdb_file(PDB_ID, pdir ='Data/AdjacencyFiles/',file_format='mmCif')
        p = MMCIF2Dict("Data/AdjacencyFiles/"+\
                    PDB_ID+'.cif')

        df = pd.DataFrame()
        df['chain'] = p['_struct_ref_seq.pdbx_strand_id']
        df['pdb'] = p['_struct_ref_seq.pdbx_PDB_id_code']
        df['pdb_start'] = p['_struct_ref_seq.pdbx_auth_seq_align_beg']
        df['pdb_end'] = p['_struct_ref_seq.pdbx_auth_seq_align_end']
        df['db_start'] = p['_struct_ref_seq.db_align_beg']
        df['db_end'] = p['_struct_ref_seq.db_align_end']

        for i in range(len(df)):
            l=[]
            l.append(int(df.iloc[i,2]))
            l.append(int(df.iloc[i,3]))
            l.append(int(df.iloc[i,4]))
            l.append(int(df.iloc[i,5]))
            dbref_dict[df.iloc[i,0]] = l 
    
    dict_id = {}
    for k,v in dbref_dict.items():
        sequence_length = v[3] - v[2] +1
        index = 0
        for i in range(sequence_length):
            value = v[0]+index
            key = k + str(value)
            dict_id[key] = value
            index +=1
    
    dict_upd_id = {}
    for k, v in dbref_dict.items():

        if v[0] == v[2] and v[1] == v[3]:
            sequence_length = v[3] - v[2] +1
#             print('No difference in residue numbers in chain',k) 
            index = 0
            for i in range(sequence_length):
                value = v[0]+index
                key = k + str(value)
                dict_upd_id[key] = v[2]+index
                index +=1
        else:
#             print('Difference in residue numbers in chain',k)
            sequence_length = v[3] - v[2] +1
            index = 0
            for i in range(sequence_length):
                value = v[0]+index
                key = k + str(value)
                dict_upd_id[key] = v[2]+index
                index +=1
                    
    return(dbref_dict, dict_upd_id, dict_id)


def intraprotein_threshold(chain_entity_dict, entity):
    '''Calculates threshold to collapse features within chains from same entity
    
    Parameters
    ----------
        chain_entity_dict : dict
            Dict obtained from entity_for_chain function
        entity : int
            given an entity identifier (1, 2, etc.) 
        
    Returns
    -------
        threshold : int
            threshold value to be used for intraprotein features
            '''
    index = 0
    for i in chain_entity_dict.values():
        if i == entity:
            index +=1
    threshold = (index*25)/100
    return(int(threshold))         
        
    
def interprotein_threshold(PDB_ID):
    '''Calculates threshold to collapse features between chains from different entities
    
    Parameters
    ----------
        PDB_ID : str
            Input a PDB ID 
        
    Returns
    -------
        threshold : int
            threshold value to be used for interprotein features
            '''
    
    d = entity_for_chain(PDB_ID)[0]
    unique_entities =[]
    for entity_val in d.values(): 
        if entity_val in unique_entities: 
            continue 
        else:
            unique_entities.append(entity_val)

    len_of_entity = []
    for e in unique_entities:
        count = 0
        for x in enumerate(d.items()):
            if x[1][1] == e:
                count +=1
        len_of_entity.append(count)
        
    #here using 50% as threshold (can be modified) 
    #len_of_entity list stores values with the number of chains for each entity and the threshold calculated by using the maximum value of this list 
    threshold = ((max(len_of_entity))*50)/100
    return(int(threshold))
    

def entity_for_chain(PDB_ID):
    '''Gives us chain and entity information for a PDB structure
    
    Parameters
    ----------
        PDB_ID : str
            Input a PDB ID 
        
    Returns
    -------
        chain_entity_dict : dict
            Dictionary with entity (values) that are identified for each chain (key) in the structure
        chainlist : list 
            List of chains present in the structure '''
    
    response1 = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}').json()
    entity = response1['rcsb_entry_container_identifiers']['polymer_entity_ids']
    chain_entity_dict= {}
    chainlist = []

    for entity_id in entity:
        polymer_entity_url = "https://data.rcsb.org/rest/v1/core/polymer_entity/" + \
                            PDB_ID + "/" + entity_id
        response2 = requests.get(polymer_entity_url).json()
        chains = response2['entity_poly']['pdbx_strand_id']
        splitchains = chains.split(',')

        for c in splitchains:
            if isinstance(c, str):
                chain_entity_dict[c] = int(entity_id) #removed E
                chainlist.append(c)

    return(chain_entity_dict, chainlist)


def makePTM_dict(pdb_ref_file, cif_file_path):
    '''Creates a dictionary of PTMs that are present in the structures of the PDB Reference file. This dictionary is used as PTM_CONTACT_DICT in contactMap globals class. 
    
    Parameters
    ----------
        pdb_ref_file : str
            Path to input the PDB reference file
        cif_file_path : str
            Path to find the cif files
            
    Returns
    -------
        ptm_dict : dict
            dictionary with PTMs and their one letter code of amino acid
            '''
        
    ptm_dict = {}
    df = pd.read_csv(pdb_ref_file)
    for i in range(len(df)):
        PDB_ID = df.iloc[i,0]
        cif_file = cif_file_path+PDB_ID+'.cif'
        cif_info = MMCIF2Dict(cif_file)
#         print(PDB_ID)

        if '_pdbx_struct_mod_residue.auth_comp_id' in cif_info:
            length_of_dict = len(cif_info['_pdbx_struct_mod_residue.auth_comp_id'])
            for j in range(length_of_dict):
                mod_res = cif_info['_pdbx_struct_mod_residue.auth_comp_id'][j]
                mod_res_parent = cif_info['_pdbx_struct_mod_residue.parent_comp_id'][j]
                ptm_dict[mod_res] = mod_res_parent
                tmp_dict[mod_res] = mod_res_parent
       
        
    for mod_res, mod_res_parent in ptm_dict.items():
        if mod_res_parent in AA_dict.keys():
            ptm_dict[mod_res] = AA_dict[mod_res_parent]
            
    return ptm_dict