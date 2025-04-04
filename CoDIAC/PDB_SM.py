import requests
import pandas as pd
import os
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

COLUMNS = ['PDB_ID', 'ENTITY_ID', 'CHAIN_ID', 'pdbx_description',  'rcsb_gene_name', 'pdbx_gene_src_scientific_name',
            'pdbx_seq_one_letter_code', 'pdbx_seq_one_letter_code_can', 'rcsb_sample_sequence_length',
            'modifications locations','modifications',  'rcsb_entity_polymer_type', 'macromolecular_type', 
            'database_name', 'database_accession', 
            'rcsb_uniprot_protein_sequence',  'entity_beg_seq_id', 'ref_beg_seq_id','aligned_regions_length',
            'mutations exist (Y/N)','rcsb_mutation_count','mutations locations', 'pdbx_mutation joined',
            'molecular_weight', 'experimental_method', 'resolution_combined', 
            'deposit_date', 'audit_author_name', 'title','pdbx_database_id_doi','pdbx_database_id_pub_med']

class PDB_interface:
    """
    A class that obtains and outputs relevant metadata/annotations for a PDB ID in the form of a dictionary
    and/or ability to write to a csv file.
    
    Attributes
    -----------
    PDB_list : list of strings
        list of PDB IDs (4-character alphanumeric identifier) that the user wants to gather information on
    
    Methods
    --------
    get_anno_dict
        Produces a dictionary full of the annotations for
         the PDB ID.
    print_output_csv
        Produces a csv file (to the location of the inputted file path) full of the annotations for all of the pdb ids
        from the inputted csv file.
    get_all
        Produces a dictionary, csv file, & Excel file (to the location of the input file path) full of the annotations
        for all of the pdb ids from the inputted csv file.
    
	
	"""
    def __init__(self, PDB_ID):
        """
        Constructor for the PDB_interface class.
        
        Parameters
        ----------
        PDB : list of strings
            PDB ID (4-character alphanumeric identifier) that the user wants to gather information on
        
        Returns
        -------
        None.

        """
        self.PDB_ID = PDB_ID

       
    def get_anno_dict(self):
        """
        Produces a dictionary full of the annotations for all of the PDB IDs from the PDB List of IDs used to instantiate the object

        Returns
        -------
        ERROR: bool
            0 if no error, 1 if error
        annotated_dict_list : list of dictionaries
            list of each entitiy dictionary 
            where the dictionary has PDB.COLUMNS keys and information for each entity
            list is empty if error
        """

        #IDs = self.PDB_list
        
        overall_dict = {}
        pdb = self.PDB_ID
#        for each_id in IDs:
        try:
            ERROR = 0
            metadata = self.PDB_metadata(self.PDB_ID)
            annotated_dict_list = metadata.set_PDB_API_annotation()
            # overall_dict[pdb] = annotated_dict_list
            # print(ERROR, annotated_dict_list)
            return ERROR, annotated_dict_list
        except:
            #print(self.PDB_ID + ' could not be fetched')
            ERROR = 1
            # print(ERROR, annotated_dict_list)
            return ERROR, []


    
    def print_output_csv(self, outputFile):
        """
        Produces a csv file (to the location of the inputted file path) full of the annotations for all of the pdb ids
        from the inputted csv file.

        Parameters
        ----------
        outputFile : string
            name of the csv file you want to generate

        Returns
        -------
        None.

        """
        ERROR, overall_dict = self.get_anno_dict()
        data_list = []
        for each_key in overall_dict.keys():
            for each_entity in overall_dict[each_key]:
                data_list.append(each_entity)
        df = pd.DataFrame(data_list, columns=COLUMNS)

        if outputFile[-4:] == '.csv':
            df.to_csv(outputFile, index=False)
        else:
            df.to_csv(outputFile + '.csv', index=False)
        
    def get_all(self):
        """
        Produces a dictionary and csv file  full of the annotations
        for all of the pdb ids from the inputted csv file.

        Returns
        -------
        overall_dict : dictionary
            dictionary full of the annotations for all of the pdb ids from the inputted file

        """
        csv_file_name = input('Enter a file name for the output csv file: ')
        ERROR, overall_dict = self.get_anno_dict()
        df = pd.DataFrame.from_dict(overall_dict)
        df.to_csv(csv_file_name + '.csv')
        return overall_dict

    class PDB_metadata:
        """
        A class that stores relevant metadata for individual PDB IDs.
        
        Attributes
        -----------
        meta_data_categories : list
            A defined list of relevant metadata categories from PDB that this class will store
        PDB_ID : string
            the PDB ID (4-character alphanumeric identifier) that the user wants to gather information on
        entry_dict : dictionary
            stores information from the entry service for the PDB_ID
        overall_polymer_entity_dict : dictionary
            stores polymer entity data for each entity id corresponding with the PDB_ID
                        
        Methods
        --------
        set_PDB_API_annotation
            Utilizes other class methods to generate dictionaries of information for a PDB ID, find the correct information
            ("value") for each attibute ("key") for that PDB ID, and store the key-value pairs in a dictionary containing
            all the desired information for that PDB ID.
        get_empty_anno_dict
            Creates a dictionary with the attributes contained in meta_data as keys with empty values.
        get_all_dicts
            Obtains information from the PDB API for a specific PDB ID and stores it in dictionaries.
        get_annotation
            Finds and returns the corresponding information for the specified attribute of the PDB ID.
        
        """

        meta_data_categories = COLUMNS

        def __init__(self, PDB_ID):
            """
            Constructor for the PDB_metadata class which takes in a PDB ID and generates global metadata dictionaries for it.
            
            Parameters
            ----------
            PDB_ID : string
                the PDB ID (4-character alphanumeric identifier) that the user wants to gather information on
            
            Returns
            -------
            None.

            """

            self.name = PDB_ID

            self.entry_dict, self.overall_polymer_entity_dict, self.overall_uniprot_dict, self.overall_nonpolymer_entity_dict, self.overall_cif_dict = self.get_all_dicts()
            self.mods = {}
            self.locs = {}

        def set_PDB_API_annotation(self):
            """
            Utilizes other class methods to generate dictionaries of information for a PDB ID, find the correct information
            ("value") for each attibute ("key") for that PDB ID, and store the key-value pairs in a dictionary containing
            all the desired information for that PDB ID.

            Returns
            -------
            anno_dict_list : list of dictionaries
                dictionary storing all information for the PDB_ID
                list is returned since multipel entities can be present in a single PDB file with different attributes

            """
            anno_dict_list = []
            polymer_entity_ids = self.entry_dict['rcsb_entry_container_identifiers']['polymer_entity_ids']
            for entity_id in polymer_entity_ids:
                anno_dict = self.return_anno_dict(entity_id, True)
                anno_dict_list.append(anno_dict)
                
            try:
                nonpolymer_entity_ids = self.entry_dict['rcsb_entry_container_identifiers']['non_polymer_entity_ids']
            except KeyError:
                nonpolymer_entity_ids = []
                
            for entity_id in nonpolymer_entity_ids:
                anno_dict = self.return_anno_dict(entity_id, False)
                anno_dict_list.append(anno_dict)
           
            # print(anno_dict_list)     
            return anno_dict_list

        def get_empty_anno_dict(self):
            """
            Creates a dictionary with the attributes contained in meta_data as keys with empty values. 

            Returns
            -------
            inner : dictionary
                placeholder dictionary created before the values are added.

            """
            inner = {}
            meta_data = self.meta_data_categories
            for meta in meta_data:
                inner[meta] = ""
            return inner
        
        def return_anno_dict(self, polymer_entity_id, is_polymer):
            """
            Creates a dictionary with the attributes contained in meta_data as keys with empty values. 

            Returns
            -------
            anno_dict : dictionary
                updated anno_dict

            """
            anno_dict = self.get_empty_anno_dict()
            if is_polymer:
                polymer_entity_dict = self.overall_polymer_entity_dict[polymer_entity_id]
            if not is_polymer:
                nonpolymer_entity_dict = self.overall_nonpolymer_entity_dict[polymer_entity_id]


            anno_dict['PDB_ID']= self.name
            anno_dict['ENTITY_ID'] = polymer_entity_id
            
            if is_polymer:
                #CHAIN ID
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                        CHAIN_ID = entity_poly['pdbx_strand_id']
                        anno_dict['CHAIN_ID'] = CHAIN_ID
                except:
                    anno_dict['CHAIN_ID'] = 'not found'
                
                #database name
                try:
                    if len(polymer_entity_dict) > 3:
                        rcsb_container_identifiers = polymer_entity_dict['rcsb_polymer_entity_container_identifiers']
                        reference_sequence_identifiers = rcsb_container_identifiers['reference_sequence_identifiers']
                        middle = reference_sequence_identifiers[0]
                    anno_dict['database_name'] = middle['database_name']
                except KeyError:
                    anno_dict['database_name'] = 'not found'
    
                #gene name
    
                #elif (attribute == 'rcsb_gene_name'):
                try:
                    rcsb_gene_nameS_list = []
                    if len(polymer_entity_dict) > 3:
                        rcsb_entity_source_organism = polymer_entity_dict['rcsb_entity_source_organism']
                        middle = rcsb_entity_source_organism[0]
                        rcsb_gene_name = middle['rcsb_gene_name']
                        for each_dict in rcsb_gene_name:
                            gene_name = each_dict['value']
                            if gene_name not in rcsb_gene_nameS_list:
                                rcsb_gene_nameS_list.append(gene_name)
                    anno_dict['rcsb_gene_name'] = '; '.join(rcsb_gene_nameS_list)
                except KeyError:
                    anno_dict['rcsb_gene_name'] = 'not found'
                   
    
                #elif (attribute == 'database_accession'):
                try:
                    if len(polymer_entity_dict) > 3:
                        rcsb_container_identifiers = polymer_entity_dict['rcsb_polymer_entity_container_identifiers']
                        reference_sequence_identifiers = rcsb_container_identifiers['reference_sequence_identifiers']
                        middle = reference_sequence_identifiers[0]
                    anno_dict['database_accession'] = middle['database_accession']
                except KeyError:
                    anno_dict['database_accession'] = 'not found'
    
                #elif (attribute == 'pdbx_seq_one_letter_code'):
                try:
                    pdbx_seq_one_letter_code = ''
                    if len(polymer_entity_dict)>3:
                        entity_poly = polymer_entity_dict['entity_poly']
                        pdbx_seq_one_letter_code = entity_poly['pdbx_seq_one_letter_code']
                    anno_dict['pdbx_seq_one_letter_code'] = pdbx_seq_one_letter_code
                except KeyError:
                    anno_dict['pdbx_seq_one_letter_code'] = 'not found'
                
                
                #elif (attribute == 'rcsb_uniprot_protein_sequence'):
                entity_id =  anno_dict['ENTITY_ID']
                uniprot_url = "https://data.rcsb.org/rest/v1/core/uniprot/" + \
                    anno_dict['PDB_ID'] + "/" + entity_id
                #print("DEBUG: %s"%(uniprot_url))
                resp = requests.get(uniprot_url)
                #print("DEBUG: %s"%(resp.status_code))
                if resp.status_code != 200:
                    #print('Failed to get %s from rcsb uniprot with code %s:'%(PDB_ID, resp.status_code))
                    anno_dict['rcsb_uniprot_protein_sequence'] = 'not found'
                else: #in case we got a response, but it does not conform to json, try, except here
                    try:
                        uniprot_entity_dict = resp.json()
                        SEQ = ''
    
                        if isinstance(uniprot_entity_dict, list): #len(uniprot_entity_dict) < 3:
                            #print(len(uniprot_entity_dict))
                            middle = uniprot_entity_dict[0]
                            if isinstance(middle, dict):
                                rcsb_uniprot_protein = middle['rcsb_uniprot_protein']
                                SEQ = rcsb_uniprot_protein['sequence']
                                #print(SEQ)
                                anno_dict['rcsb_uniprot_protein_sequence'] = SEQ
                            else:
                                anno_dict['rcsb_uniprot_protein_sequence'] = 'not found'
                    except KeyError:
                        anno_dict['rcsb_uniprot_protein_sequence'] = 'not found'
                
                #elif (attribute == 'rcsb_sample_sequence_length'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                        sequence_length = entity_poly['rcsb_sample_sequence_length']
                    anno_dict['rcsb_sample_sequence_length'] = sequence_length
                except KeyError:
                    anno_dict['rcsb_sample_sequence_length'] = 'not found'
                
                
                #elif (attribute == 'rcsb_entity_polymer_type'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                    anno_dict['rcsb_entity_polymer_type'] = entity_poly['rcsb_entity_polymer_type']
                except KeyError:
                    anno_dict['rcsb_entity_polymer_type'] = 'not found'
                
    
                #elif (attribute == 'macromolecular_type'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                    anno_dict['macromolecular_type'] = entity_poly['type']
                except KeyError:
                    anno_dict['macromolecular_type'] = 'not found'
                
                #elif (attribute == 'molecular_weight'):
                entry_info = self.entry_dict['rcsb_entry_info']
                anno_dict['molecular_weight'] = entry_info['molecular_weight']
    
                #elif (attribute == 'experimental_method'):
                entry_info = self.entry_dict['rcsb_entry_info']
                anno_dict['experimental_method'] = entry_info['experimental_method']
                
                #elif (attribute == 'resolution_combined'):
                try:
                    summary = self.entry_dict['rcsb_entry_info']
                    # need to check if more than one value is ever reported
                    anno_dict['resolution_combined'] = summary['resolution_combined'][0]
                except KeyError:
                    anno_dict['resolution_combined'] = 'not found'
    
                #elif (attribute == 'pdbx_seq_one_letter_code_can'):
                try:
                    SEQ = ''
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                        SEQ = entity_poly['pdbx_seq_one_letter_code_can']
                    anno_dict['pdbx_seq_one_letter_code_can']= SEQ
                except KeyError:
                    anno_dict['pdbx_seq_one_letter_code_can']= 'not found'
                
                #elif (attribute == 'entity_beg_seq_id'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_align = polymer_entity_dict['rcsb_polymer_entity_align']
                        middle = entity_align[0]
                        aligned_info = middle['aligned_regions']
                        middle2 = aligned_info[0]
                    anno_dict['entity_beg_seq_id'] = str(middle2['entity_beg_seq_id'])
                except KeyError:
                    anno_dict['entity_beg_seq_id'] = 'not found'
                
                #elif (attribute == 'ref_beg_seq_id'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_align = polymer_entity_dict['rcsb_polymer_entity_align']
                        middle = entity_align[0]
                        aligned_info = middle['aligned_regions']
                        middle2 = aligned_info[0]
                    anno_dict['ref_beg_seq_id'] = str(middle2['ref_beg_seq_id'])
                except KeyError:
                    anno_dict['ref_beg_seq_id'] = 'not found'
                
                #elif (attribute == 'aligned_regions_length'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_align = polymer_entity_dict['rcsb_polymer_entity_align']
                        middle = entity_align[0]
                        aligned_info = middle['aligned_regions']
                        middle2 = aligned_info[0]
                    anno_dict['aligned_regions_length'] = str(middle2['length'] )
                except KeyError:
                    anno_dict['aligned_regions_length'] = 'not found'
                
                #elif (attribute == 'pdbx_gene_src_scientific_name'):
                try:
                    if  len(polymer_entity_dict) > 3:
                        entity_src_gen = polymer_entity_dict['entity_src_gen']
                        middle = entity_src_gen[0]
                        species = middle['pdbx_gene_src_scientific_name']
                    anno_dict['pdbx_gene_src_scientific_name'] = species
                except KeyError:
                    anno_dict['pdbx_gene_src_scientific_name'] = 'not found'
    
                #elif (attribute == 'mutations exist (Y/N)'):
                mutation_count = 0
                if len(polymer_entity_dict) > 3:
                    entity_poly = polymer_entity_dict['entity_poly']
                    mutation_count = entity_poly['rcsb_mutation_count']
                if mutation_count != 0:
                    anno_dict['mutations exist (Y/N)'] = 'Y'
                else:
                    anno_dict['mutations exist (Y/N)'] = 'N'
                    
    
                #elif (attribute =='rcsb_mutation_count'):
                if len(polymer_entity_dict) > 3:
                    entity_poly = polymer_entity_dict['entity_poly']
                    anno_dict['rcsb_mutation_count'] = entity_poly['rcsb_mutation_count']
                #return mutation_count
                
                #elif (attribute =='modifications locations'):
                try:
                    if len(polymer_entity_dict) > 3:
                        locations = self.get_mod_locations(polymer_entity_id=polymer_entity_id)
                        if len(locations)>0:
                            anno_dict['modifications locations'] = '; '.join(locations)
                        else:
                            anno_dict['modifications locations'] = 'N/A'
                except KeyError:
                    anno_dict['modifications locations'] = 'N/A'
                
    
                #elif (attribute =='modifications'):
                try:
                    if len(polymer_entity_dict) > 3:
                        types = self.get_mod_types(polymer_entity_id=polymer_entity_id)
                        if len(types)>0:
                            anno_dict['modifications'] = '; '.join(types)
                        else:
                            anno_dict['modifications'] = 'N/A'
                except KeyError:
                    anno_dict['modifications'] = 'N/A'   
    
               # elif (attribute =='mutations locations'):
                try:
                    if len(polymer_entity_dict) > 3:
                        polymer_entity = polymer_entity_dict['rcsb_polymer_entity']
                        mutation = polymer_entity['pdbx_mutation'].split(', ')
                        locations = []
                        for item in mutation:
                            locations.append(item[1:(len(item)-1)])
                        anno_dict['mutations locations'] = '; '.join(locations)
                except KeyError:
                     anno_dict['mutations locations'] = 'N/A'
                
                #elif (attribute == 'pdbx_mutation joined'):
                try:
                    if len(polymer_entity_dict) > 3:
                        polymer_entity = polymer_entity_dict['rcsb_polymer_entity']
                        mutation = polymer_entity['pdbx_mutation']
                    anno_dict['pdbx_mutation joined'] = mutation
                except KeyError:
                    anno_dict['pdbx_mutation joined'] = 'N/A'
                
                #elif (attribute =='deposit_date'):
                summary = self.entry_dict['rcsb_accession_info']
                anno_dict['deposit_date'] = summary['deposit_date'][:10]            
    
                #elif (attribute == 'audit_author_name'): #fixed
                AUTHORS_LIST = []
                temp_authors = self.entry_dict['audit_author']
                for each in temp_authors:
                    name = each['name']
                    AUTHORS_LIST.append(name)
                anno_dict['audit_author_name'] = ', '.join(AUTHORS_LIST)
    
                #elif (attribute == 'title'):
                citation = self.entry_dict['rcsb_primary_citation']
                anno_dict['title'] = citation['title']
    
                #elif (attribute =='pdbx_database_id_doi'):
                try:
                    citation = self.entry_dict['rcsb_primary_citation']
                    DOI = citation['pdbx_database_id_doi']
                    anno_dict['pdbx_database_id_doi'] = DOI
                except KeyError:
                    try:
                        citation = self.entry_dict['citation']
                        middle = citation[0]
                        anno_dict['pdbx_database_id_doi'] = middle['pdbx_database_id_doi']
                    except KeyError:
                        anno_dict['pdbx_database_id_doi'] = 'not found'
                
                #elif (attribute =='pdbx_database_id_pub_med'):
                try:
                    citation = self.entry_dict['rcsb_primary_citation']
                    anno_dict['pdbx_database_id_pub_med'] = citation['pdbx_database_id_pub_med']
                except KeyError:
                    try:
                        citation = self.entry_dict['citation']
                        middle = citation[0]
                        anno_dict['pdbx_database_id_pub_med'] = middle['pdbx_database_id_PubMed']
                    except KeyError:
                        anno_dict['pdbx_database_id_pub_med'] = 'not found'
                
                #elif (attribute == "pdbx_description"):
                try:
                    if len(polymer_entity_dict) > 3:
                        rcsb_polymer_entity = polymer_entity_dict['rcsb_polymer_entity']
                    anno_dict['pdbx_description'] = rcsb_polymer_entity['pdbx_description']
                except KeyError:
                    anno_dict['pdbx_description'] =  'N/A'

            if not is_polymer:
                #CHAIN ID
                try:
                   # insert for loop here
                   each_index = 0
                   chains = []                       
                   for each_index in (self.overall_cif_dict).keys():
                       cur_entity_id = self.overall_cif_dict[each_index]['entity_id']
                       if cur_entity_id == polymer_entity_id:            
                           chains.append(self.overall_cif_dict[each_index]['chain']) # string of all chains for that entity id 
                   anno_dict['CHAIN_ID'] = (', '.join(set(chains)))
                except KeyError:
                   anno_dict['CHAIN_ID'] = 'NA'
                
                #database name
                anno_dict['database_name'] = 'not found'
    
                #gene name
                try:
                   nonpoly_name = nonpolymer_entity_dict['rcsb_nonpolymer_entity_container_identifiers']['chem_ref_def_id']
                   anno_dict['rcsb_gene_name'] = nonpoly_name
                except KeyError:
                   anno_dict['rcsb_gene_name'] = 'NA'
    
                   
    
                #elif (attribute == 'database_accession'):
                anno_dict['database_accession'] = 'not found'
    
                #elif (attribute == 'pdbx_seq_one_letter_code'):
                anno_dict['pdbx_seq_one_letter_code'] = 'not found'
                
                
                #elif (attribute == 'rcsb_uniprot_protein_sequence'):
                anno_dict['rcsb_uniprot_protein_sequence'] = 'not found'
                
                #elif (attribute == 'rcsb_sample_sequence_length'):
                try:
                   max_length = nonpolymer_entity_dict['rcsb_nonpolymer_entity_feature_summary'][0]['maximum_length']
                   print('max lenght')
                   anno_dict['rcsb_sample_sequence_length'] = max_length
                except KeyError:
                    anno_dict['rcsb_sample_sequence_length'] = 'not found'
                
                
                #elif (attribute == 'rcsb_entity_polymer_type'):
                anno_dict['rcsb_entity_polymer_type'] = 'non-polymer'
                
    
                #elif (attribute == 'macromolecular_type'):
                anno_dict['macromolecular_type'] = 'non-polymer'
                
                #elif (attribute == 'molecular_weight'):
                formula_weight = nonpolymer_entity_dict['rcsb_nonpolymer_entity']['formula_weight']
                anno_dict['molecular_weight'] = formula_weight
    
                #elif (attribute == 'experimental_method'):
                entry_info = self.entry_dict['rcsb_entry_info']
                anno_dict['experimental_method'] = entry_info['experimental_method']
                
                #elif (attribute == 'resolution_combined'):
                try:
                    summary = self.entry_dict['rcsb_entry_info']
                    # need to check if more than one value is ever reported
                    anno_dict['resolution_combined'] = summary['resolution_combined'][0]
                except KeyError:
                    anno_dict['resolution_combined'] = 'not found'
    
                #elif (attribute == 'pdbx_seq_one_letter_code_can'):
                anno_dict['pdbx_seq_one_letter_code_can']= 'not found'
                
                #elif (attribute == 'entity_beg_seq_id'):
                try:
                    # insert for loop here
                    each_index = 0
                    seq_ids = []
                    for each_index in (self.overall_cif_dict).keys():
                        #print(each_index)
                        cur_entity_id = self.overall_cif_dict[each_index]['entity_id']
                        #print(cur_entity_id)
                        if cur_entity_id == polymer_entity_id:
                            #print(cur_entity_id + ' from .cif, ' + polymer_entity_id + ' from og')
                            #seq_id = self.overall_cif_dict[each_index]['resid']
                            seq_ids.append(self.overall_cif_dict[each_index]['resid'])
                    seq_ids_int = [eval(i) for i in seq_ids]
                    min_val = min(seq_ids_int)
                    max_val = max(seq_ids_int)
                    anno_dict['entity_beg_seq_id'] = min_val
                    anno_dict['ref_beg_seq_id'] = min_val
                    if min_val == max_val:
                        anno_dict['aligned_regions_length'] = 0
                    else:
                        anno_dict['aligned_regions_length'] = max_val - min_val + 1
                    # anno_dict['entity_beg_seq_id'] = ', '.join(seq_ids)
                except KeyError:
                    anno_dict['entity_beg_seq_id'] = 'not found'
                    anno_dict['ref_beg_seq_id'] = 'NA'
                    anno_dict['aligned_regions_length'] = 'NA'
                    
                        
                    
                
                #elif (attribute == 'ref_beg_seq_id'):
                # anno_dict['ref_beg_seq_id'] = 'not found'
                
                #elif (attribute == 'aligned_regions_length'):
                # try:
                #     num_features = nonpolymer_entity_dict['rcsb_nonpolymer_entity_feature_summary'][0]['count']
                #     anno_dict['aligned_regions_length'] = num_features
                # except KeyError:
                #     anno_dict['aligned_regions_length'] = 0
                
                #elif (attribute == 'pdbx_gene_src_scientific_name'):
                # try:
                #     # insert for loop here
                #     each_index = 0
                #     names = []
                #     for each_index in (self.overall_cif_dict).keys():
                #         #print(each_index)
                #         cur_entity_id = self.overall_cif_dict[each_index]['entity_id']
                #         #print(cur_entity_id)
                #         if cur_entity_id == polymer_entity_id:
                #             #print(cur_entity_id + ' from .cif, ' + polymer_entity_id + ' from og')
                #             #name = self.overall_cif_dict[each_index]['resn']
                #             names.append(self.overall_cif_dict[each_index]['resn'])
                #     anno_dict['pdbx_gene_src_scientific_name'] = ', '.join(names)
                # except KeyError:
                #     anno_dict['pdbx_gene_src_scientific_name'] = 'not found'

                try:
                   nonpoly_name = nonpolymer_entity_dict['rcsb_nonpolymer_entity_container_identifiers']['chem_ref_def_id']
                   anno_dict['rcsb_gene_name'] = nonpoly_name
                except KeyError:
                   anno_dict['rcsb_gene_name'] = 'NA'
    
                #elif (attribute == 'mutations exist (Y/N)'):
                anno_dict['mutations exist (Y/N)'] = 'NA'
                    
    
                #elif (attribute =='rcsb_mutation_count'):
                anno_dict['rcsb_mutation_count'] = 'NA'
                
                #elif (attribute =='modifications locations'):
                anno_dict['modifications locations'] = 'N/A'
                
    
                #elif (attribute =='modifications'):
                anno_dict['modifications'] = 'N/A'   
    
               # elif (attribute =='mutations locations'):
                anno_dict['mutations locations'] = 'N/A'
                
                #elif (attribute == 'pdbx_mutation joined'):
                anno_dict['pdbx_mutation joined'] = 'N/A'
                
                #elif (attribute =='deposit_date'):
                summary = self.entry_dict['rcsb_accession_info']
                anno_dict['deposit_date'] = summary['deposit_date'][:10]            
    
                #elif (attribute == 'audit_author_name'): #fixed
                AUTHORS_LIST = []
                temp_authors = self.entry_dict['audit_author']
                for each in temp_authors:
                    name = each['name']
                    AUTHORS_LIST.append(name)
                anno_dict['audit_author_name'] = ', '.join(AUTHORS_LIST)
    
                #elif (attribute == 'title'):
                citation = self.entry_dict['rcsb_primary_citation']
                anno_dict['title'] = citation['title']
    
                #elif (attribute =='pdbx_database_id_doi'):
                try:
                    citation = self.entry_dict['rcsb_primary_citation']
                    DOI = citation['pdbx_database_id_doi']
                    anno_dict['pdbx_database_id_doi'] = DOI
                except KeyError:
                    try:
                        citation = self.entry_dict['citation']
                        middle = citation[0]
                        anno_dict['pdbx_database_id_doi'] = middle['pdbx_database_id_doi']
                    except KeyError:
                        anno_dict['pdbx_database_id_doi'] = 'not found'
                
                #elif (attribute =='pdbx_database_id_pub_med'):
                try:
                    citation = self.entry_dict['rcsb_primary_citation']
                    anno_dict['pdbx_database_id_pub_med'] = citation['pdbx_database_id_pub_med']
                except KeyError:
                    try:
                        citation = self.entry_dict['citation']
                        middle = citation[0]
                        anno_dict['pdbx_database_id_pub_med'] = middle['pdbx_database_id_PubMed']
                    except KeyError:
                        anno_dict['pdbx_database_id_pub_med'] = 'not found'
                
                #elif (attribute == "pdbx_description"):

                try:
                    description = nonpolymer_entity_dict['rcsb_nonpolymer_entity']['pdbx_description']
                    anno_dict['pdbx_description'] =  description
                except KeyError:
                    anno_dict['pdbx_description'] =  ''

            

            return anno_dict


           
    
        def find_mods(self, polymer_entity_id):
            polymer_entity_dict = self.overall_polymer_entity_dict[polymer_entity_id]
            pdb_seq = polymer_entity_dict['entity_poly']['pdbx_seq_one_letter_code']
            start = 0
            offset = 0
            mods = []
            locs = []
            while pdb_seq.find('(', start) != -1:
                index = pdb_seq.find('(', start)
                end = pdb_seq.find(')', index)
                mod = pdb_seq[index+1:end]
                loc = index + 1 - offset
                mod_loc = mod+'-'+str(loc)
                mods.append(mod_loc)
                locs.append(str(loc))
                offset = offset + len(mod) + 1
                start = end
            return mods, locs
    
        def get_mod_locations(self, polymer_entity_id):
            if polymer_entity_id in self.locs:
                locs = self.locs[polymer_entity_id]
            else:
                mods, locs = self.find_mods(polymer_entity_id)
                self.mods[polymer_entity_id] = mods
                self.locs[polymer_entity_id] = locs
            return locs
        
        def get_mod_types(self, polymer_entity_id):
            if polymer_entity_id in self.mods:
                mods = self.mods[polymer_entity_id]
            else:
                mods, locs = self.find_mods(polymer_entity_id)
                self.mods[polymer_entity_id] = mods
                self.locs[polymer_entity_id] = locs                
            return mods

        def get_all_dicts(self):
            """
            Obtains information from the PDB API for a specific PDB ID and stores it in dictionaries.

            Returns
            -------
            entry_dict : dictionary
                stores information from the entry service for the PDB_ID
            pubmed_dict : dictionary
                stores the pubmed annotations (data integrated from PubMed) for the PDB_ID
            schema_dict : dictionary
                stores information from the entry schema for the PDB_ID
            schema_uniprot_dict : dictionary
                stores information from the uniprot schema for the PDB_ID
            overall_polymer_entity_dict : dictionary
                stores polymer entity data for each entity id corresponding with the PDB_ID

            """

            PDB_ID = self.name

            ERROR = 0
            #check for repsonse issues, collecting errors and handling them
            entry_url = "https://data.rcsb.org/rest/v1/core/entry/" + PDB_ID
            resp = requests.get(entry_url)
            if resp.status_code != 200:
                #print('Failed to get %s from rcsb entry with code %s:'%(PDB_ID, resp.status_code))
                ERROR = 1
            entry_dict = resp.json() #THIS one is used!!
        
                  
            entry_id_dict = entry_dict['entry']
            ENTRY_ID = entry_id_dict['id']
            rcsb_entry_container_identifiers = entry_dict['rcsb_entry_container_identifiers']
            ENTITY_IDS = rcsb_entry_container_identifiers['entity_ids']
            polymer_entity_ids = rcsb_entry_container_identifiers['polymer_entity_ids']
            # nonpolymer_entity_ids = []
            # for each in ENTITY_IDS:
            #     if each not in polymer_entity_ids:
            #         nonpolymer_entity_ids.append(each)
            
            if 'non_polymer_entity_ids' in rcsb_entry_container_identifiers:
                nonpolymer_entity_ids = rcsb_entry_container_identifiers['non_polymer_entity_ids']
            else:
                nonpolymer_entity_ids = []
                
            overall_polymer_entity_dict = {}
            for each_entity_id in polymer_entity_ids:
                polymer_entity_url = "https://data.rcsb.org/rest/v1/core/polymer_entity/" + \
                    ENTRY_ID + "/" + each_entity_id
                resp = requests.get(polymer_entity_url)
                if resp.status_code != 200:
                    #print('Failed to get %s from rcsb polymer entity with code %s:'%(PDB_ID, resp.status_code))
                    ERROR = 1
                polymer_entity_dict = resp.json()
                overall_polymer_entity_dict[each_entity_id] = polymer_entity_dict


            overall_uniprot_dict = {}
            for each_entity_id in ENTITY_IDS:
                uniprot_url = "https://data.rcsb.org/rest/v1/core/uniprot/" + \
                    ENTRY_ID + "/" + each_entity_id
                resp = requests.get(uniprot_url)
                uniprot_dict = resp.json()
                overall_uniprot_dict[each_entity_id] = uniprot_dict

            overall_nonpolymer_entity_dict = {}
            for each_entity_id in nonpolymer_entity_ids:
                nonpolymer_entity_url = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/" + \
                    ENTRY_ID + "/" + each_entity_id
                resp = requests.get(nonpolymer_entity_url)
                nonpolymer_entity_dict = resp.json()
                overall_nonpolymer_entity_dict[each_entity_id] = nonpolymer_entity_dict

            cif_path = '/Users/adk9hq/Documents/Research_UVA/metabolism_codiac/met_ciffiles'
            mmCIF_file_path = cif_path + '/' + PDB_ID + '.cif'
            
            if not os.path.exists(mmCIF_file_path):
                PDB.download_cifFile([PDB_ID],cif_path)
            
            pdb_info = MMCIF2Dict(mmCIF_file_path)
            
            each_pdb_cif_dict = {}
            overall_cif_dict = {}
            if '_pdbx_nonpoly_scheme.auth_mon_id' in pdb_info.keys():
                res = pdb_info['_pdbx_nonpoly_scheme.auth_mon_id']
                for i in range(len(res)):
                    each_res_cif_dict = {}
                    if 'HOH' not in (pdb_info['_pdbx_nonpoly_scheme.pdb_mon_id'][i]):
                        # PDB residue number
                        each_res_cif_dict['resid'] = pdb_info['_pdbx_nonpoly_scheme.pdb_seq_num'][i]
                        
                        # PDB strand/chain id.
                        each_res_cif_dict['chain'] = pdb_info['_pdbx_nonpoly_scheme.pdb_strand_id'][i]
                        
                        # NDB/RCSB residue number
                        each_res_cif_dict['ndb_num'] = pdb_info['_pdbx_nonpoly_scheme.ndb_seq_num'][i]
                        
                        # PDB residue identifier
                        each_res_cif_dict['resn'] = pdb_info['_pdbx_nonpoly_scheme.pdb_mon_id'][i]
                        
                        # Pointer to _atom_site.label_entity_id.
                        # The value of _entity.id must uniquely identify a record in the ENTITY list.
                        # Note that this item need not be a number; it can be any unique identifier.
                        each_res_cif_dict['entity_id'] = pdb_info['_pdbx_nonpoly_scheme.entity_id'][i]
                        
                        each_pdb_cif_dict[i] = each_res_cif_dict
                
                #overall_cif_dict[PDB_ID] = each_pdb_cif_dict
                
            overall_cif_dict = each_pdb_cif_dict
        
        
            return entry_dict, overall_polymer_entity_dict, overall_uniprot_dict, overall_nonpolymer_entity_dict, overall_cif_dict

        
def generateStructureRefFile_fromUniprotFile(uniprotRefFile, outputFile):
    '''
    Creates a PDB Structure Reference File

    Parameters
    ----------
        uniprotRefFile: str
            Location of the uniprot referenct to use
        outputFile: str
            name of the output file
        
    Returns
    -------
        Writes to the outputFile location a CSV formatted, annotated structure file
    '''
    #check that a uniprotRefFile exists
    if not os.path.isfile(uniprotRefFile):
        print("ERROR: %s cannot be found"%(uniprotRefFile))
        return
    uniprot_df = pd.read_csv(uniprotRefFile)
    PDB_IDs = []
    for index, row in uniprot_df.iterrows():
        pdbs_col = row['PDB IDs']
   
        if isinstance(pdbs_col, str):
            pdbs = pdbs_col.split(';')
            for pdb_id in pdbs:
                PDB_IDs.append(pdb_id)
    #print("DEBUG: PDB IDs")
    #print(PDB_IDs)
    generateStructureRefFile(PDB_IDs, outputFile)    

def progress_bar(current, total, bar_length=20):
    fraction = current / total

    arrow = int(fraction * bar_length - 1) * '-' + '>'
    padding = int(bar_length - len(arrow)) * ' '

    ending = '\n' if current == total else '\r'

    print(f'Progress: [{arrow}{padding}] {int(fraction*100)}%', end=ending)

def generateStructureRefFile(PDB_IDs, outputFile):
    '''
    Creates a PDB Structure Reference File

    Parameters
    ----------
        PDB_IDs : a list of PDB IDs
        outputFile : name of the output file
        
    Returns
    -------
        .csv Structure reference file with relevant metadata
    '''
    dict_list = []
    bad_PDBs = []
    for PDB_ID in PDB_IDs:
        #This can have a lot of fetching time, we would like to print to let the user
        # know what the status is of the job fetch. 
        progress_bar(PDB_IDs.index(PDB_ID), len(PDB_IDs))
        interface = PDB_interface(PDB_ID)
        ERROR, annotations = interface.get_anno_dict()
        # reiteration attempt - in case of connectivity issue. 
        num_attempts = 0
        max_attempts = 5
        while ERROR and num_attempts < max_attempts:
            ERROR, annotations = interface.get_anno_dict()
            num_attempts += 1
        if not ERROR:
            for annotation in annotations: #multiple entities can come back and each should become something in a dataframe
                dict_list.append(annotation)
        else:
            bad_PDBs.append(PDB_ID)
    df = pd.DataFrame(dict_list, columns=COLUMNS)
    df.to_csv(outputFile, index=False)
    
    print('Structure Reference File successfully created!')
    if(len(bad_PDBs) == 0):
        print('All PDBs successfully fetched')
    else:
        print("Could not fetch the following PDBs, these encountered errors, despite retrying %d times:"%(max_attempts))
        print(bad_PDBs)
    #interface.print_output_csv(outputFile)
    #return dict_list, bad_PDBs
    
def download_cifFile(PDB_list, PATH):
    '''generates .cif files for PDB structures
    
    Parameters
    ----------
        PDB_list : list of PDB IDs 
        PATH : path to the store the .cif file
        
    Returns
    -------
        .cif files to the path specified'''
    
    for PDB_ID in PDB_list:
        
        response=requests.get(f'https://files.rcsb.org/view/{PDB_ID}.cif')
        cif_path = PATH+'/'+PDB_ID+'.cif'

        with open(cif_path, 'wb') as f:
            f.write(response.content)
        
    return cif_path