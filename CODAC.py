import os
import requests
import numpy as np
import pandas as pd
import logging
import ast
import json

class Interpro:
    '''Given an Interpro Accession ID, this class returns a list of Uniprot Accession IDs associated to this Interpro Acc ID'''
    
    def __init__(self, IPR_ID, specie='Homo sapiens'):
        self.IPR_ID = IPR_ID
        self.specie = specie
        
    def fetch_uniprotid(self):
        
        fetch_all_results = []
        
        base_url = f"https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/InterPro/{self.IPR_ID}/?page_size=200"

        logging.basicConfig(filename='error.log',
                           format='%(asctime)s %(message)s',
                           encoding='utf-8',
                           level=logging.WARNING)
        try:
            response = requests.get(base_url)
            data = response.json()

            fetch_all_results = fetch_all_results + data['results']

            while data['next'] is not None:

                print("Found next page and downloading", data['next'])
                response = requests.get(data['next'])
                data = response.json()
                fetch_all_results = fetch_all_results + data['results']

        except requests.exceptions.RequestException as err:
            print('Found request exceptions...')
            logging.warning(err)
            
        UNIPROT_ID_LIST = []
        
        for entry in fetch_all_results:
            Uniprot_ACC = (entry['metadata']['accession'])
            Scientific_name = entry['metadata']['source_organism']['scientificName']
            fullname = (entry['metadata']['source_organism']['fullName'] )

            if Scientific_name == self.specie:

                UNIPROT_ID_LIST.append(Uniprot_ACC)
        print('Found {count} Uniprot IDs linked to {IPR_ID}'.format(count=len(UNIPROT_ID_LIST), IPR_ID=self.IPR_ID))        
        return(UNIPROT_ID_LIST)
    
class Uniprot:
    '''Given a Uniprot ID, this class can be used to fetch reference structure attributes such as Domains, Domain Architecture, Gene name, Reference sequence from Uniprot for specie of interest and returns a Domain reference csv file'''
    
    def __init__(self, UPROT_ID):
        self.UPROT_ID = UPROT_ID
        
    def Domains(self):
        res = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{self.UPROT_ID}').json() 
        gene = res['id']
        get_dom = []
        if 'features' in res.keys():
            for i in range(len(res['features'])):
                s = res['features'][i]
                for k, v in s.items():
                    if k == 'type':
                        if v == 'DOMAIN':
                            start = s['begin']
                            end = s['end']
                            name = s['description']
                            header = self.UPROT_ID+':'+gene+':'+name+':'+start+':'+end
                            #print(header)
                            get_dom.append(header)                    
        else:
            get_dom.append('None')
        return(get_dom)
    
    def DomainArchitecture(self):
        res = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{self.UPROT_ID}').json() 
        gene = res['id']
        get_dom = []
        domdict = {}
        if 'features' in res.keys():
            for i in range(len(res['features'])):
                s = res['features'][i]
                for k, v in s.items():
                    if k == 'type':
                        if v == 'DOMAIN':
                            start = s['begin']
                            end = s['end']
                            name = s['description']
                            header = self.UPROT_ID+':'+gene+':'+name+':'+start+':'+end
                            get_dom.append(header)  
                            domdict[start] = end, name
                    
        sorted_dict = dict(sorted(domdict.items(),reverse=True))
        domain_arch = []
        for key, value in sorted_dict.items():
            sort_start = key
            sort_end = value[0]
            domain = value[1]
            domain_arch.append(domain)
            
        if len(domain_arch) > 1:
            final_domarch = '|'.join(domain_arch)   
        else:
            final_domarch = domain
        #print('{count} number of domains with Architecture {Arch}'.format(count=len(domain_arch),Arch = final_domarch))  
        return(final_domarch)
    
    def Gene(self):
        res = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{self.UPROT_ID}').json()
        gene = res['id']
        return(gene)
    
    def RefSeq(self):
        res = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{self.UPROT_ID}').json()
        if 'sequence' in res.keys():
            seq = (res['sequence']['sequence'])
        else:
            seq = 'None'                 
        return(seq)
        
    def makeRefFile(list_of_UID, IPR_ID, Specie='Human'):
        '''Makes a Domain Reference file
        
        Parameters
        ----------
            list_of_UID : a list of Uniprot Accession IDs
            IPR_ID : Interpro Accession 
            Specie : Scientific name of organism of interest 
            
        Attributes
        ----------
            UPROT_Accession : Uniprot Accession ID
            Gene : gene name
            Domain Boundaries : reference domain boundary ranges with its start and end positions
            Ref Sequence : Reference sequence
            Domain Architecture : Domains found within the protein sequence arranged from N ter to C ter
            
        Returns
        -------
            .csv Domain reference file with all attributes'''
            
        Domain_list = []
        Gene_list = []
        Refseq_list = []
        Domain_Arch = []
        print('Fetch successful for Uniprot IDs :\n')
        for ID in list_of_UID:
            instance = Uniprot(ID)
            DOMAIN = Uniprot.Domains(instance)
            GENE = Uniprot.Gene(instance)
            REFSEQ = Uniprot.RefSeq(instance)
            DOMAIN_ARCH = Uniprot.DomainArchitecture(instance)
            
            Domain_list.append(DOMAIN)
            Gene_list.append(GENE)
            Refseq_list.append(REFSEQ)
            Domain_Arch.append(DOMAIN_ARCH)
            print(ID)
            
        df = pd.DataFrame(columns=['UPROT Accession','Gene','Domain Boundaries','Ref Sequence','Domain Architecture'])
        df['UPROT Accession'] = list_of_UID
        df['Gene'] = Gene_list
        df['Domain Boundaries'] = Domain_list
        df['Ref Sequence'] = Refseq_list
        df['Domain Architecture'] = Domain_Arch
        df.index = np.arange(1, len(df) + 1)
        outputfile = 'Domain_referencer_'+IPR_ID+'_'+Specie+'.csv'
        df.to_csv(outputfile)
        print('Domain Reference File successfully created!')
            
    def make_fastafile(DOMAIN, input_RefFile, IPR_ID, specie='HUMAN'):

        '''Makes fasta file for domain of interest

        Parameters
        ----------
            DOMAIN : str ('SH2', 'SH3', 'PH', etc.) This function searches for this domain string keyword in the domain name 
            input_RefFile : DomainReference file generated using write_referencer_csv
            interpro_Acc : str 

        Returns
        -------
            outfile_fasta : FASTA file with domain reference sequences
            Returns number of reference domains found '''

        df = pd.read_csv(input_RefFile)
        outfile_fasta = 'RefFasta_'+DOMAIN+'_'+IPR_ID+'_'+specie+'.fasta'
        COUNT = 0

        with open(outfile_fasta,'w') as f:

            for i in range(len(df)):

                DomBound = df['Domain Boundaries'][i]
                DomBoundList = ast.literal_eval(DomBound)
                DomainCount = len(DomBoundList)
                RefSeq = df['Ref Sequence'][i]

                for entry in DomBoundList: 

                    if entry == 'None':
                        print('Domain boundary information not found for Acc ID {UniprotID}! '
                                  .format(UniprotID = df['UPROT Accession'][i]))

                    else:
                        UniprotID, gene, domain_name, ref_start, ref_end = entry.split(':') 

                        if DOMAIN in domain_name:

                            domain_identifier = []

                            HEADER = UniprotID+'|'+gene+'|'+ref_start+'|'+ref_end
                            Domain_RefSeq = RefSeq[int(ref_start)-1:int(ref_end)]
                            f.write('>'+str(HEADER)+'\n'+str(Domain_RefSeq)+'\n')
                            COUNT += 1

        return('{count} number of {DOMAIN_NAME} domains found'.format(count=COUNT,DOMAIN_NAME = DOMAIN))

class fromArpeggio:
    
    def __init__(self, PATH, input_json, OUTPATH):
        self.PATH = PATH
        self.input_json = (input_json)
        self.OUTPATH = OUTPATH
        self.PDB_ID = (self.input_json).split('.')[0]
          
    
    def Create_ContactMap(self):
        ''' This function generates a contactmap text file.
        
        Parameters
        ----------
            input_json : input a .json file generated from running Arpeggio
            PATH : PATH to find the .json input file
            OUTPATH : PAth to save the output .txt contactmap file
            
        Returns
        -------
            Contactmap {PDB_ID}.txt to specified OUTPATH'''
        
        filepath = self.PATH + '/' + self.input_json
        with open(filepath,'r') as file:
            data = json.load(file)

        contact_types=['aromatic','carbonyl','hbond','hydrophobic','ionic','polar','proximal',
                       'vdw','vdw_clash','weak_hbond','weak_polar','xbond']
        
        outfile = self.OUTPATH+'/'+self.PDB_ID+'.txt'
        
        if os.path.isfile(filepath):
        
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

                    if interaction == 'INTRA_SELECTION':

                            if distance < 5:

                                for contact in contactlist:
                                    if contact in contact_types:

                                        if resnum_1 > resnum_2:

                                            atom_pair = atom_1 +'-'+atom_2 

                                            f.write(str(self.PDB_ID)+'\t'+str(chain_1)+'\t'+str(chain_2)+'\t'+
                                                    str(res_1)+'\t'+str(resnum_1)+'\t'+str(res_2)+'\t'+str(resnum_2)+'\t'+
                                                   str(atom_pair)+'\t'+str(distance)+'\t'+str(contactlist)+'\n')

                                        else:
                                            atom_pair = atom_2 +'-'+atom_1
                                            f.write(str(self.PDB_ID)+'\t'+str(chain_2)+'\t'+str(chain_1)+'\t'+
                                                    str(res_2)+'\t'+str(resnum_2)+'\t'+str(res_1)+'\t'+str(resnum_1)+'\t'+
                                                   str(atom_pair)+'\t'+str(distance)+'\t'+str(contactlist)+'\n')

                                        break

            df = pd.read_csv(outfile,  sep='\t')

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

            df.replace({"Res1": AA_dict},inplace=True) 
            df.replace({"Res2": AA_dict},inplace=True) 

            df.sort_values(["Chain1","Chain2","ResNum1","ResNum2"], 
                                axis=0,
                                ascending=[True, True, True, True], 
                                inplace=True)

            df.to_csv(outfile,sep='\t',encoding='utf-8',index=False) 
            print('CM generated for {PDB}'.format(PDB= self.PDB_ID))
            
        else:
            print('File Not Found!')



