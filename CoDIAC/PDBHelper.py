import pandas as pd
import numpy as np
from string import digits


# This is a set of helpers to parse the file of information that annotates PDB IDs with sequences 
class PDBEntitiesClass:
    
    def __init__(self, df, PDB_ID):
        """

        Parameters
        ----------
        df: pandas DataFrame
            dataframe of annotated PDB files from contact mapping project

        PDB_ID: str
            PDB_ID of interest to pull data for from df

        Returns
        -------
        PDBEntitiesClass
            initiated PDB info class 

        Attributes
        ----------
        PDB_ID: str
            storing the PBD ID of this object
        entities: list
            list of ints of the entities that are contained in a single
        entities_dict: dict
            dictionary of PDBClass objects for each of the entities in a PDB 

        """
        #for each line in the dataframe for a PDB_ID, return a PDBClass object of information 

        self.pdb_dict = {}

        sub_df = df[df['PDB_ID']==PDB_ID]
        if sub_df.empty:
            print("ERROR: PDB ID %s not found in dataframe"%(sub_df))
            return -1


        for index, row in sub_df.iterrows():
            entity_id = int(row['ENTITY_ID'])
            resolution = row['resolution_combined']
            method = row['experimental_method']
            pdbclass = PDBClass(row)
            self.pdb_dict[entity_id] = pdbclass

        self.resolution = resolution
        self.method = method







class PDBClass:
    def __init__(self, row):
        """

        Given a row (i.e. a specific entity) of data from the dataframe describing PDBs and their annotations, parse the row into a class that is useful

        Parameters
        ----------
        row: pandas row
            row of annotated PDB files from contact mapping project

        Returns
        -------
        PDBClass
            initiated PDB info class 

        Attributes
        ----------
        PDB_ID: str
            storing the PBD ID of this object
        database_name: str
            accession database
        acc: str
            accession, usually uniprot, but stored in database_name
        chain_list: list
            list of chain characters found for that row (i.e. entity)
        structure_seq: str
            Sequence of the structure as defined in PDB
        entity_macroType: str
            type of macromolecule as defined in PDB
        species: str
            species as defined in PDB
        ref_seq: str
            the reference sequence that was the focus of protein, as extracted (and only as reference, does not contain mutations)
            may have less than structure_sequence as it will not include tags, etc.
        ref_seq_mutated: str
            the reference sequence as it was crystallized in the PDB, including mutations, should allow for alignment of positions in contact map to ref_seq positions
        ref_seq_positions: list
            list of integers defining the positions within the acc and ref_sequence 
        gene_name: str
            name of gene/protein of reference
        domain_dict: dict
            Has keys that lists the domain num, this is for uniqueness, and is a dict of dict
            second key is domain name, points to a list of [start, stop, num_gaps, num_muts]
        domains_of_interest: dict
            will be empty if no domains_of_interest are noted, otherwise will be a dictionary of 
            headers and sequences for printing to fasta.
        struct_domain_arch: str
            domain architecture caught fully within the crystal structure
        protein_domain_arch: str
            the full protein domain architectures
        """
        self.database_name = row['database_name']
        self.acc = row['database_accession']
        self.chain_list = self.return_chain_list(row)
        self.structure_seq = row['ref:struct/ref sequence']
        self.entity_macroType = row['macromolecular_type']
        self.species = row['pdbx_gene_src_scientific_name']
        self.PDB_ID = row['PDB_ID']
        self.Entity_ID = row['ENTITY_ID']
        self.ref_seq = row['rcsb_uniprot_protein_sequence']
        self.gene_name = row['ref:gene name']
        ref_sequence_positions = row['ref:reference range']
        
        if ref_sequence_positions != 'not found':
            ref_pos_list = ref_sequence_positions.split('-')
            self.ref_seq_positions = [int(ref_pos_list[0]), int(ref_pos_list[1])]

        domains = row['ref:domains']
        domain_dict = {}
        if domains != '-1':
            if isinstance(domains, str):
                domain_list = domains.split(';')

                domain_num = 0
                for domain in domain_list:
                    name, IPR, other = domain.split(':')
                    start, stop, gaps, muts = other.split(',')
                    small_dict = {}
                    small_dict[name] = [int(start), int(stop), int(gaps), int(muts)]
                    domain_dict[domain_num] = small_dict
                    domain_num += 1
        self.domains = domain_dict #will be empty dict if no domains


#         domains_of_interest = row['SH2 domain sequence from struct']
#         if isinstance(domains_of_interest, str):
#             headers = row['reference']
#             DOI_list = domains_of_interest.split(';')
#             DOI_headers = headers.split(';')
#             DOI_dict = {}
#             if len(DOI_headers) != len(DOI_list):
#                 print("ERROR: domains of interest have issue, different number of headers and sequences")
#                 self.domains_of_interest = -1
#             else:
#                 temp_dict = {}
#                 for index in range(0, len(DOI_headers)):
#                     temp_dict[DOI_headers[index]] = DOI_list[index]
#                 self.domains_of_interest = temp_dict
#         else:
#             self.domains_of_interest = {}

        self.struct_domain_arch = row['ref:struct domain architecture']

        self.protein_domain_arch = row['ref:protein domain architecture']
        self.ref_seq_mutated = row['ref:struct/ref sequence']

#         # BELOW HERE IS if we are working with an annotation file that has been processed for contact mapping
#         if 'struct seq ext' in row:
#             self.struct_seq_ext = row['struct seq ext']
#             self.ERROR_CODE = row['ERROR_CODE']
#             self.offset = row['offset']
        PTM_str = row['modifications']
        if isinstance(PTM_str, str):
            self.transDict = return_PTM_dict(PTM_str)
        else:
            self.transDict = {}





    def return_chain_list(self, row):
        """
        Return the chains of a PDB_ID as a list of characters
        if there are two entries for the PDB_ID, return as a list of lists, with Uniprot IDs to help sort which is which
        """
        chains_list = []

        chains = row['CHAIN_ID']
        chains= chains.replace("'", '')
        chains = chains.replace(" ", '')
        remove_digits = str.maketrans('', '', digits)
        chains= chains.translate(remove_digits)
        chains = chains.replace(",", '')
        chain_vals = []
        for ind in range(0, len(chains)):
            chain_vals.append(chains[ind])
        return chain_vals

def return_PTM_dict(PTM_str):
    """
    GIven a string from the table of PTMs, return a dictionary as it was setup in contact map checking
    """
    transDict = {}
    list_vals = PTM_str.split(';')
    for val in list_vals:
        if '-' in val:
            PTM, res = val.split('-')
            transDict[int(res)] = PTM

    return transDict
