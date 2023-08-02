import requests
import pandas as pd
import csv
import os
import errno

class PDB_interface:
    """
    A class that obtains and outputs relevant metadata/annotations for individual PDB IDs in the form of a dictionary
    and/or a csv file.
    
    Attributes
    -----------
    PDB_list : list of strings
        list of PDB IDs (4-character alphanumeric identifier) that the user wants to gather information on
    
    Methods
    --------
    get_anno_dict
        Produces a dictionary full of the annotations for all of the PDB IDs from the inputted csv file.
    print_output_csv
        Produces a csv file (to the location of the inputted file path) full of the annotations for all of the pdb ids
        from the inputted csv file.
    get_all
        Produces a dictionary, csv file, & Excel file (to the location of the input file path) full of the annotations
        for all of the pdb ids from the inputted csv file.
    
	
	"""
    def __init__(self, PDB_list):
        """
        Constructor for the PDB_interface class.
        
        Parameters
        ----------
        PDB_list : list of strings
            list of PDB IDs (4-character alphanumeric identifier) that the user wants to gather information on
        
        Returns
        -------
        None.

        """
        self.PDB_list = PDB_list

       
    def get_anno_dict(self):
        """
        Produces a dictionary full of the annotations for all of the PDB IDs from the PDB List of IDs used to instantiate the object

        Returns
        -------
        overall_dict : dictionary
            dictionary full of the annotations for all of the pdb ids from the inputted file

        """

        IDs = self.PDB_list
        
        overall_dict = {}
        
        for each_id in IDs:
            try:
                metadata = self.PDB_metadata(each_id)
                annotated_dict = metadata.set_PDB_API_annotation()
                overall_dict[each_id] = annotated_dict
            except:
                print(each_id + ' is not a valid PDB ID')
        return overall_dict


    
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
        columns = ['PDB_ID', 'ENTITY_ID', 'ENTITY_DESCRIPTION', 'CHAIN_ID', 'DATABASE', 'GENE_NAME', 
                   'ACCESS', 'PDB_SEQ', 'UNIPROT_SEQ', 'CHAIN_LENGTH', 'POLYMER_TYPE', 'MACROMOLECULAR_TYPE', 
                   'MOLECULAR_WEIGHT', 'EXPERIMENT_TYPE', 'RESOLUTION', 'CANNONICAL_REF_SEQ', 
                   'PDB_SEQ_BEG_POSITION', 'CANNONICAL_SEQ_BEG_POSITION','REF_SEQ_LENGTH',
                   'SPECIES', 'MUTATIONS/MODS (Y/N)', 'MUTATIONS/MODS (#)', 'MUTATIONS (LOCATION)',
                   'MUTATIONS (TYPE)', 'MODIFICATIONS (LOCATION)', 'MODIFICATIONS (TYPE)', 'DEPOSITED (DATE)', 
                   'DEPOSITED (AUTHORS)', 'TITLE', 'DOI', 'PUBMED_ID']
        overall_dict = self.get_anno_dict()
        data_list = []
        for each_key in overall_dict.keys():
            for each_entity in overall_dict[each_key]:
                data_list.append(each_entity)
        df = pd.DataFrame(data_list, columns=columns)

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
        overall_dict = self.get_anno_dict()
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
        pubmed_dict : dictionary
            stores the pubmed annotations (data integrated from PubMed) for the PDB_ID
        schema_dict : dictionary
            stores information from the entry schema for the PDB_ID
        schema_uniprot_dict : dictionary
            stores information from the uniprot schema for the PDB_ID
        overall_polymer_entity_dict : dictionary
            stores polymer entity data for each entity id corresponding with the PDB_ID
        overall_uniprot_dict : dictionary
            stores UniProt annotations for a given macromolecular entity (each entity id corresponding with the PDB_ID)
                        
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

        meta_data_categories = ['PDB_ID', 'ENTITY_ID', 'ENTITY_DESCRIPTION', 'CHAIN_ID', 'DATABASE', 'GENE_NAME', 'ACCESS', 
                                'PDB_SEQ', 'UNIPROT_SEQ', 'CHAIN_LENGTH', 'POLYMER_TYPE', 'MACROMOLECULAR_TYPE', 'MOLECULAR_WEIGHT',
                                'EXPERIMENT_TYPE', 'RESOLUTION', 'CANNONICAL_REF_SEQ', 'PDB_SEQ_BEG_POSITION', 
                                'CANNONICAL_SEQ_BEG_POSITION', 'REF_SEQ_LENGTH', 'SPECIES', 'MUTATIONS/MODS (Y/N)', 'MUTATIONS/MODS (#)', 
                                'MUTATIONS (LOCATION)', 'MUTATIONS (TYPE)','MODIFICATIONS (LOCATION)', 'MODIFICATIONS (TYPE)', 
                                'DEPOSITED (DATE)', 'DEPOSITED (AUTHORS)', 'TITLE', 'DOI', 'PUBMED_ID']

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

            self.entry_dict, self.pubmed_dict, self.schema_dict, self.schema_uniprot_dict, self.overall_polymer_entity_dict, self.overall_uniprot_dict = self.get_all_dicts()
            self.mods = {}
            self.locs = {}

        def set_PDB_API_annotation(self):
            """
            Utilizes other class methods to generate dictionaries of information for a PDB ID, find the correct information
            ("value") for each attibute ("key") for that PDB ID, and store the key-value pairs in a dictionary containing
            all the desired information for that PDB ID.

            Returns
            -------
            anno_dict : dictionary
                dictionary storing all information for the PDB_ID

            """
            anno_dict_list = []
            for entity_id in self.entry_dict['rcsb_entry_container_identifiers']['polymer_entity_ids']:
                anno_dict = self.get_empty_anno_dict()
                
                for attribute in anno_dict:
                    anno_dict[attribute] = self.get_annotation(attribute, entity_id)

                anno_dict_list.append(anno_dict)
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
            overall_uniprot_dict : dictionary
                stores UniProt annotations for a given macromolecular entity (each entity id corresponding with the PDB_ID)

            """

            PDB_ID = self.name

            entry_url = "https://data.rcsb.org/rest/v1/core/entry/" + PDB_ID
            resp = requests.get(entry_url)
            entry_dict = resp.json()
        
            pubmed_url = "https://data.rcsb.org/rest/v1/core/pubmed/" + PDB_ID
            resp = requests.get(pubmed_url)
            pubmed_dict = resp.json()
        
            schema_url = "https://data.rcsb.org/rest/v1/schema/entry" + PDB_ID
            resp = requests.get(schema_url)
            schema_dict = resp.json()
        
            schema_uniprot_url = "https://data.rcsb.org/rest/v1/schema/uniprot" + PDB_ID
            resp = requests.get(schema_uniprot_url)
            schema_uniprot_dict = resp.json()
        
            entry_id_dict = entry_dict['entry']
            ENTRY_ID = entry_id_dict['id']
            rcsb_entry_container_identifiers = entry_dict['rcsb_entry_container_identifiers']
            ENTITY_IDS = rcsb_entry_container_identifiers['entity_ids']
            polymer_entity_ids = rcsb_entry_container_identifiers['polymer_entity_ids']
            nonpolymer_entity_ids = []
            for each in ENTITY_IDS:
                if each not in polymer_entity_ids:
                    nonpolymer_entity_ids.append(each)
        
            overall_polymer_entity_dict = {}
            for each_entity_id in polymer_entity_ids:
                polymer_entity_url = "https://data.rcsb.org/rest/v1/core/polymer_entity/" + \
                    ENTRY_ID + "/" + each_entity_id
                resp = requests.get(polymer_entity_url)
                polymer_entity_dict = resp.json()
                overall_polymer_entity_dict[each_entity_id] = polymer_entity_dict
        
            overall_uniprot_dict = {}
            for each_entity_id in ENTITY_IDS:
                uniprot_url = "https://data.rcsb.org/rest/v1/core/uniprot/" + \
                    ENTRY_ID + "/" + each_entity_id
                resp = requests.get(uniprot_url)
                uniprot_dict = resp.json()
                overall_uniprot_dict[each_entity_id] = uniprot_dict
        
            return entry_dict, pubmed_dict, schema_dict, schema_uniprot_dict, overall_polymer_entity_dict, overall_uniprot_dict

        def get_annotation(self, attribute, polymer_entity_id):
            """
            Finds and returns the corresponding information for the specified attribute of the PDB ID.

            Parameters
            ----------
            attribute : string
                each attribute is a specific piece of information wanted for the PDB_ID (attributes match the "keys" of the
                dictionaries created for each PDB_ID)

            Returns
            -------
            string
                string containing the corresponding information for that attribute.
                or string containing an error or N/A message when the desired information cannot be found/does not exist
            int
                int containing the corresponding information for that attribute.
                
            Raises
            --------
            KeyError if the desired information cannot be accessed, meaning a different pathway to access that information must be tried
            or that the desired information cannot be found/does not exist

            """
            
            PDB_ID = self.name
            polymer_entity_dict = self.overall_polymer_entity_dict[polymer_entity_id]
            if (attribute == 'PDB_ID'):
                return PDB_ID
            elif (attribute == 'ENTITY_ID'):
                return polymer_entity_id
            elif (attribute == 'CHAIN_ID'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                        CHAIN_ID = entity_poly['pdbx_strand_id']
                    return CHAIN_ID
                except KeyError:
                    return 'not found'
            elif (attribute == 'DATABASE'):
                try:
                    if len(polymer_entity_dict) > 3:
                        rcsb_container_identifiers = polymer_entity_dict['rcsb_polymer_entity_container_identifiers']
                        reference_sequence_identifiers = rcsb_container_identifiers['reference_sequence_identifiers']
                        middle = reference_sequence_identifiers[0]
                    return middle['database_name']
                except KeyError:
                    return 'not found'
            elif (attribute == 'GENE_NAME'):
                try:
                    GENE_NAMES_list = []
                    if len(polymer_entity_dict) > 3:
                        rcsb_entity_source_organism = polymer_entity_dict['rcsb_entity_source_organism']
                        middle = rcsb_entity_source_organism[0]
                        rcsb_gene_name = middle['rcsb_gene_name']
                        for each_dict in rcsb_gene_name:
                            gene_name = each_dict['value']
                            if gene_name not in GENE_NAMES_list:
                                GENE_NAMES_list.append(gene_name)
                    return '; '.join(GENE_NAMES_list)
                except KeyError:
                    return 'not found'
                    # try:
                    #     GENE_NAMES_list = []
                    #     for each_dict in self.overall_polymer_entity_dict.values():
                    #         if len(each_dict) > 3:
                    #             entity_src_gen = each_dict['entity_src_gen']
                    #             middle = entity_src_gen[0]
                    #             all_names = (middle['pdbx_gene_src_gene']).split(', ')
                    #             for each in all_names:
                    #                 GENE_NAMES_list.append(each)
                    #     return ','.join(GENE_NAMES_list)
                    # except KeyError:
                    #     try:
                    #         for each_dict in self.overall_polymer_entity_dict.values():
                    #             if len(each_dict) > 3:
                    #                 rcsb_entity_source_organism = each_dict['rcsb_entity_source_organism']
                    #                 middle = rcsb_entity_source_organism[0]
                    #         return middle['source_type']
                    #     except KeyError:
                    #         return 'not found'
            elif (attribute == 'ACCESS'):
                try:
                    if len(polymer_entity_dict) > 3:
                        rcsb_container_identifiers = polymer_entity_dict['rcsb_polymer_entity_container_identifiers']
                        reference_sequence_identifiers = rcsb_container_identifiers['reference_sequence_identifiers']
                        middle = reference_sequence_identifiers[0]
                    return middle['database_accession']
                except KeyError:
                    return 'not found'
            elif (attribute == 'PDB_SEQ'):
                try:
                    PDB_SEQ = ''
                    if len(polymer_entity_dict)>3:
                        entity_poly = polymer_entity_dict['entity_poly']
                        PDB_SEQ = entity_poly['pdbx_seq_one_letter_code']
                    return PDB_SEQ
                except KeyError:
                    return 'not found'
            elif (attribute == 'UNIPROT_SEQ'):
                try:
                    uniprot_entity_dict = self.overall_uniprot_dict[polymer_entity_id]
                    SEQ = ''
                    if len(uniprot_entity_dict) < 3:
                        middle = uniprot_entity_dict[0]
                        rcsb_uniprot_protein = middle['rcsb_uniprot_protein']
                        SEQ = rcsb_uniprot_protein['sequence']
                    return SEQ
                except KeyError:
                    return 'not found'
            elif (attribute == 'CHAIN_LENGTH'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                        sequence_length = entity_poly['rcsb_sample_sequence_length']
                    return sequence_length
                except KeyError:
                    return 'not found'
            elif (attribute == 'POLYMER_TYPE'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                    return entity_poly['rcsb_entity_polymer_type']
                except KeyError:
                    return 'not found'
            elif (attribute == 'MACROMOLECULAR_TYPE'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                    return entity_poly['type']
                except KeyError:
                    return 'not found'
            elif (attribute == 'MOLECULAR_WEIGHT'):
                entry_info = self.entry_dict['rcsb_entry_info']
                return entry_info['molecular_weight']
            elif (attribute == 'EXPERIMENT_TYPE'):
                entry_info = self.entry_dict['rcsb_entry_info']
                return entry_info['experimental_method']
            elif (attribute == 'RESOLUTION'):
                try:
                    summary = self.entry_dict['rcsb_entry_info']
                    # need to check if more than one value is ever reported
                    return summary['resolution_combined'][0]
                except KeyError:
                    return 'not found'
            # elif (attribute == 'PDB_REF_SEQ'):
            #     try:

            #         for each_dict in self.overall_uniprot_dict.values():
            #             if len(polymer_entity_dict) < 3:
            #                 middle = each_dict[0]
            #                 rcsb_uniprot_protein = middle['rcsb_uniprot_protein']
            #                 SEQ = rcsb_uniprot_protein['sequence']
            #         REF_SEQ_POSITIONS_DICT = {'Beginning Position on PDB Sequence': '', 'Ref Sequence Length': ''}
            #         for each_dict in self.overall_polymer_entity_dict.values():
            #             if len(each_dict) > 3:
            #                 entity_align = each_dict['rcsb_polymer_entity_align']
            #                 middle = entity_align[0]
            #                 aligned_info = middle['aligned_regions']
            #                 middle2 = aligned_info[0]
            #                 REF_SEQ_POSITIONS_DICT['Beginning Position on PDB Sequence'] = middle2['entity_beg_seq_id']
            #                 REF_SEQ_POSITIONS_DICT['Ref Sequence Length'] = middle2['length']
            #         beg_number = REF_SEQ_POSITIONS_DICT['Beginning Position on PDB Sequence']-1
            #         end_number = beg_number + REF_SEQ_POSITIONS_DICT['Ref Sequence Length']
            #         return SEQ[beg_number:end_number]
            #     except KeyError:
            #         return 'not found'
            elif (attribute == 'CANNONICAL_REF_SEQ'):
                try:
                    SEQ = ''
                    if len(polymer_entity_dict) > 3:
                        entity_poly = polymer_entity_dict['entity_poly']
                        SEQ = entity_poly['pdbx_seq_one_letter_code_can']
                    return SEQ
                    # uniprot_entity_dict = self.overall_uniprot_dict[polymer_entity_id]
                    # if len(uniprot_entity_dict) < 3:
                    #     middle = each_dict[0]
                    #     rcsb_uniprot_protein = middle['rcsb_uniprot_protein']
                    #     SEQ = rcsb_uniprot_protein['sequence']
                    # REF_SEQ_POSITIONS_DICT = {'Beginning Position on Cannonical Sequence': '', 'Ref Sequence Length': ''}
                    # if len(polymer_entity_dict) > 3:
                    #     entity_align = polymer_entity_dict['rcsb_polymer_entity_align']
                    #     middle = entity_align[0]
                    #     aligned_info = middle['aligned_regions']
                    #     middle2 = aligned_info[0]
                    #     REF_SEQ_POSITIONS_DICT['Beginning Position on Cannonical Sequence'] = middle2['ref_beg_seq_id']
                    #     REF_SEQ_POSITIONS_DICT['Ref Sequence Length'] = middle2['length']
                    # beg_number = REF_SEQ_POSITIONS_DICT['Beginning Position on Cannonical Sequence']-1
                    # end_number = beg_number + REF_SEQ_POSITIONS_DICT['Ref Sequence Length']
                    # return SEQ[beg_number:end_number]
                except KeyError:
                    return 'not found'
            elif (attribute == 'PDB_SEQ_BEG_POSITION'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_align = polymer_entity_dict['rcsb_polymer_entity_align']
                        middle = entity_align[0]
                        aligned_info = middle['aligned_regions']
                        middle2 = aligned_info[0]
                    return str(middle2['entity_beg_seq_id'])
                except KeyError:
                    return 'not found'
            elif (attribute == 'CANNONICAL_SEQ_BEG_POSITION'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_align = polymer_entity_dict['rcsb_polymer_entity_align']
                        middle = entity_align[0]
                        aligned_info = middle['aligned_regions']
                        middle2 = aligned_info[0]
                    return str(middle2['ref_beg_seq_id'])
                except KeyError:
                    return 'not found'
            elif (attribute == 'REF_SEQ_LENGTH'):
                try:
                    if len(polymer_entity_dict) > 3:
                        entity_align = polymer_entity_dict['rcsb_polymer_entity_align']
                        middle = entity_align[0]
                        aligned_info = middle['aligned_regions']
                        middle2 = aligned_info[0]
                    return str(middle2['length'] )
                except KeyError:
                    return 'not found'
            elif (attribute == 'SPECIES'):
                try:
                    if  len(polymer_entity_dict) > 3:
                        entity_src_gen = polymer_entity_dict['entity_src_gen']
                        middle = entity_src_gen[0]
                        species = middle['pdbx_gene_src_scientific_name']
                    return species
                except KeyError:
                    return 'not found'
            elif (attribute == 'MUTATIONS/MODS (Y/N)'):
                mutation_count = 0
                if len(polymer_entity_dict) > 3:
                    entity_poly = polymer_entity_dict['entity_poly']
                    mutation_count = entity_poly['rcsb_mutation_count']
                if mutation_count != 0:
                    return 'Y'
                else:
                    return 'N'
            elif (attribute == 'MUTATIONS/MODS (#)'):
                if len(polymer_entity_dict) > 3:
                    entity_poly = polymer_entity_dict['entity_poly']
                    mutation_count = entity_poly['rcsb_mutation_count']
                return mutation_count
            elif (attribute == 'MODIFICATIONS (LOCATION)'):
                try:
                    if len(polymer_entity_dict) > 3:
                        locations = self.get_mod_locations(polymer_entity_id=polymer_entity_id)
                        if len(locations)>0:
                            return '; '.join(locations)
                        else:
                            return 'N/A'
                except KeyError:
                    return 'N/A'
            elif (attribute == 'MODIFICATIONS (TYPE)'):
                try:
                    if len(polymer_entity_dict) > 3:
                        types = self.get_mod_types(polymer_entity_id=polymer_entity_id)
                        if len(types)>0:
                            return '; '.join(types)
                        else:
                            return 'N/A'
                except KeyError:
                    return 'N/A'                
            elif (attribute == 'MUTATIONS (LOCATION)'):
                try:
                    if len(polymer_entity_dict) > 3:
                        polymer_entity = polymer_entity_dict['rcsb_polymer_entity']
                        mutation = polymer_entity['pdbx_mutation'].split(', ')
                        locations = []
                        for item in mutation:
                            locations.append(item[1:(len(item)-1)])
                        location = '; '.join(locations)
                    return location
                except KeyError:
                    return 'N/A'
            elif (attribute == 'MUTATIONS (TYPE)'):
                try:
                    if len(polymer_entity_dict) > 3:
                        polymer_entity = polymer_entity_dict['rcsb_polymer_entity']
                        mutation = polymer_entity['pdbx_mutation']
                    return mutation
                except KeyError:
                    return 'N/A'
            elif (attribute == 'DEPOSITED (DATE)'):
                summary = self.entry_dict['rcsb_accession_info']
                DEPOSITED_DATE = summary['deposit_date'][:10]
                return DEPOSITED_DATE
            elif (attribute == 'DEPOSITED (AUTHORS)'): #fixed
                AUTHORS_LIST = []
                temp_authors = self.entry_dict['audit_author']
                for each in temp_authors:
                    name = each['name']
                    AUTHORS_LIST.append(name)
                final_string = ', '.join(AUTHORS_LIST)
                return final_string
            elif (attribute == 'TITLE'):
                citation = self.entry_dict['rcsb_primary_citation']
                TITLE = citation['title']
                return TITLE
            elif (attribute == 'DOI'):
                try:
                    citation = self.entry_dict['rcsb_primary_citation']
                    DOI = citation['pdbx_database_id_doi']
                    return DOI
                except KeyError:
                    try:
                        citation = self.entry_dict['citation']
                        middle = citation[0]
                        DOI = middle['pdbx_database_id_doi']
                    except KeyError:
                        return 'not found'
            elif (attribute == 'PUBMED_ID'):
                try:
                    citation = self.entry_dict['rcsb_primary_citation']
                    PUBMED_ID = citation['pdbx_database_id_pub_med']
                    return PUBMED_ID
                except KeyError:
                    try:
                        citation = self.entry_dict['citation']
                        middle = citation[0]
                        DOI = middle['pdbx_database_id_PubMed']
                    except KeyError:
                        return 'not found'
            elif (attribute == "ENTITY_DESCRIPTION"):
                try:
                    if len(polymer_entity_dict) > 3:
                        rcsb_polymer_entity = polymer_entity_dict['rcsb_polymer_entity']
                    return rcsb_polymer_entity['pdbx_description']
                except KeyError:
                    return 'N/A'
                    

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

    interface = PDB_interface(PDB_IDs)
    interface.print_output_csv(outputFile)
    #interface.get_output_xlsx(outputFile)
    print('Structure Reference File successfully created!')
    
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