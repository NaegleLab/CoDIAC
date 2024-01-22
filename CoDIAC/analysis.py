from CoDIAC import PDBHelper, contactMap
from CoDIAC import pTyrLigand_helpers as pTyr_helpers
import pandas as pd
import numpy as np
from Bio import SeqIO
import os

def CanonicalFeatures(pdb_ann_file, ADJFILES_PATH, reference_fastafile, error_structures_list, append_refSeq = True, PTM='PTR', mutation=False, 
                      domain_of_interest='SH2', SH2_file='SH2_C', PTM_file='pTyr_C'):
    '''Generates contact features that are present across canonical interfaces (between a domain and it sligand partner)
    
    Parameters
    ----------
        pdb_ann_file : str
            PDB reference file with all PDB structures annotated and filtered based on the domain of interest
        ADJFILES_PATH : str
            path to fetch adjacency files
        reference_fastafile : str
            fasta file with reference sequences of domain of interest obtained from the Uniprot reference csv file
        error_structures_list : list
            list of PDB structures that are present in the PDB reference file but not useful for contactmap analysis due to issues in the PDB structure (discontinuous chains), error generate adjacency files, unable to assign a reference sequence, etc.
        PTM : str
            PTM that binds to our domain of interest
        append_refSeq : boolean
            appending reference sequences to the fasta output file
        mutation : boolean
            fetches native/mutant structures. Default set to retrieve native structures
        domain_of_interest : str
            the domain of interest
        SH2_file : str
            name of fasta and feature files for the domain of interest
        PTM_file : str
            name of fasta and feature files for the ligand entities with PTMs on it
            
    Returns
    -------
        Fasta and feature files with canonical interface contact features'''
    
    main = pd.read_csv(pdb_ann_file)
    list_of_uniprotids=[]
    for name, group in main.groupby('PDB_ID'):
        PDB_ID = name
        print(PDB_ID)

        if PDB_ID not in error_structures_list:

            for index, row in group.iterrows():

                if isinstance(row['modifications'], str):
                    
                    transDict = PDBHelper.return_PTM_dict(row['modifications'])
                    for res in transDict:
                        if PTM in transDict[res]:
                            lig_entity = row['ENTITY_ID']                             

                            entities = PDBHelper.PDBEntitiesClass(main, PDB_ID)
                            for entity in entities.pdb_dict.keys():
                                domains = entities.pdb_dict[entity].domains
                                for domain_num in domains:
                                    if domain_of_interest in domains[domain_num]:
                                        SH2_entity = entity 
                                        get_mutation = (pd.isnull(main.loc[(main['ENTITY_ID'] == SH2_entity) & (main['PDB_ID'] == PDB_ID ), 'ref:variants']))
                                        check_mutation = get_mutation.values.tolist()
                                        df2 = main.loc[(main['ENTITY_ID'] == SH2_entity) & (main['PDB_ID'] == PDB_ID )]
                                        uniprot_id = df2['database_accession'].values.tolist()


                    transDict_stripped = []
                    for i in transDict.values():
                        i = i.strip()
                        transDict_stripped.append(i)
                    
                    if PTM in transDict.values():
                        if check_mutation != mutation:
                            list_of_uniprotids.append(uniprot_id[0])

                            pdbClass = entities.pdb_dict[lig_entity]
                            dict_of_lig = contactMap.return_single_chain_dict(main, PDB_ID, ADJFILES_PATH, lig_entity)
                            dict_of_SH2 = contactMap.return_single_chain_dict(main, PDB_ID, ADJFILES_PATH, SH2_entity)

                            cm_aligned = dict_of_lig['cm_aligned']
                            if hasattr(cm_aligned, 'refseq'):
                                value = True
                            else:
                                value = False

                            for res in cm_aligned.transDict:
                                if res in cm_aligned.resNums:
                                    if PTM in cm_aligned.transDict[res]: #print the aligned sequence
                                        res_start, res_end, aligned_str, tick_labels = pTyr_helpers.return_pos_of_interest(
                                            cm_aligned.resNums, cm_aligned.structSeq, res, n_term_num=5, c_term_num=5, PTR_value = 'y')

                                        from_dict = dict_of_lig
                                        to_dict = dict_of_SH2
                                        adjList, arr = contactMap.return_interChain_adj(ADJFILES_PATH, from_dict, to_dict)
                                        adjList_alt, arr_alt = contactMap.return_interChain_adj(ADJFILES_PATH, to_dict, from_dict)

                                        domains = dict_of_SH2['pdb_class'].domains
                                        for domain_num in domains:
                                            if domain_of_interest in str(domains[domain_num]):

                                                dom_header = list(domains[domain_num].keys())[0]
                                                SH2_start, SH2_stop, muts, gaps = domains[domain_num][dom_header]


                                                arr_sub, list_aa_from_sub, list_to_aa_sub = contactMap.return_arr_subset_by_ROI(arr, 
                                                                 res_start, res_end, from_dict['cm_aligned'].return_min_residue(), 
                                                                 SH2_start, SH2_stop, to_dict['cm_aligned'].return_min_residue())

                #                                             
                                                fasta_header = makeHeader(PDB_ID, SH2_entity,int(SH2_start), int(SH2_stop),domain_of_interest,pdb_ann_file, reference_fastafile)+'|lig_'+str(res)+'|'+PDB_ID


                                                if lig_entity == SH2_entity:


                                                    cm_aligned.print_fasta_feature_files(SH2_start, SH2_stop, res_start, res_end,
                                                             fasta_header, 'pTyr',PTM_file, append=True, 
                                                                use_ref_seq_aligned=value)

                                                    cm_aligned.print_fasta_feature_files(res_start, res_end, SH2_start, SH2_stop,
                                                             fasta_header, 'SH2',SH2_file, append=True, 
                                                                use_ref_seq_aligned=value)


                                                if lig_entity != SH2_entity:
                                                    if hasattr(to_dict['cm_aligned'], 'refseq'):
                                                        contactMap.print_fasta_feature_files(arr_alt, to_dict['cm_aligned'].refseq, 
                                                            SH2_start, SH2_stop, to_dict['cm_aligned'].return_min_residue(), 
                                                            res_start, res_end, from_dict['cm_aligned'].return_min_residue(),
                                                            fasta_header,'pTyr', PTM_file, threshold=1, append=True)
                                                    else:
                                                        contactMap.print_fasta_feature_files(arr_alt, to_dict['cm_aligned'].structSeq, 
                                                            SH2_start, SH2_stop, to_dict['cm_aligned'].return_min_residue(), 
                                                            res_start, res_end, from_dict['cm_aligned'].return_min_residue(),
                                                            fasta_header,'pTyr', PTM_file, threshold=1, append=True)


                                                    contactMap.print_fasta_feature_files(arr, from_dict['cm_aligned'].structSeq, 
                                                                 res_start, res_end, from_dict['cm_aligned'].return_min_residue(), 
                                                                 SH2_start, SH2_stop, to_dict['cm_aligned'].return_min_residue(),
                                                                fasta_header,'SH2', SH2_file, threshold=1, append=True )


                                                  
    if append_refSeq:                                
        inputfile = PTM_file+'.fasta'
        with open(inputfile, 'a') as file:
            for query_uniprotid in set(list_of_uniprotids):
                fasta_seq = SeqIO.parse(open(reference_fastafile), 'fasta')
                for fasta in fasta_seq:
                    name, sequence = fasta.id, str(fasta.seq)
                    ref_uniprot_id, ref_gene,ref_domain, ref_index, ref_ipr, ref_start, ref_end = name.split('|')
                    if query_uniprotid == ref_uniprot_id:
                        file.write('>'+name+'\n'+sequence+'\n')


def NonCanonicalFeatures(pdb_ann_file, ADJFILES_PATH, reference_fastafile, error_structures_list, append_refSeq=True, mutation = False, DOMAIN = 'SH2', 
                         filename='SH2_NC'):
    '''Generates contact features that are present across non-canonical interfaces (between two domains part of teh same protein) 
    
    Parameters
    ----------
        pdb_ann_file : str
            PDB reference file with all PDB structures annotated and filtered based on the domain of interest
        ADJFILES_PATH : str
            path to fetch adjacency files
        reference_fastafile : str
            fasta file with reference sequences of domain of interest obtained from the Uniprot reference csv file
        error_structures_list : list
            list of PDB structures that are present in the PDB reference file but not useful for contactmap analysis due to issues in the PDB structure (discontinuous chains), error generate adjacency files, unable to assign a reference sequence, etc.
        append_refSeq : boolean
            appending reference sequences to the fasta output file
        mutation : boolean
            fetches native/mutant structures. Default set to retrieve native structures
        DOMAIN : str
            the domain of interest
        filename : str
            name of fasta and feature files 
        
    Returns
    -------
        Fasta and feature files with non-canonical interface contact features'''
    
    ann = pd.read_csv(pdb_ann_file)
    list_of_uniprotids = []
    for index, row in ann.iterrows():

        PDB_ID = row['PDB_ID']
        gene = str(row['ref:gene name'])
        entity_id = row['ENTITY_ID']
        domain = str(row['ref:domains'])
        parse_domain = domain.split(';')
        species = str(row['pdbx_gene_src_scientific_name'])
        uniprot_ID = str(row['database_accession'])
        check_mutation = (pd.isnull(ann.loc[index, 'ref:variants']))
        
        print('------',PDB_ID,'------')

        if PDB_ID not in error_structures_list:
            entities = PDBHelper.PDBEntitiesClass(ann, PDB_ID)
            pdbClass = entities.pdb_dict[entity_id] #this holds information about the protein crystalized, such as domains
            chain = contactMap.chainMap(PDB_ID, entity_id)
            chain.construct(ADJFILES_PATH)

            if species != 'not found':

                if check_mutation != mutation:
                    print('mutation',row['ref:variants'])

                    if gene != 'N/A (Not Found In Reference)':

                        if domain != 'nan' and len(parse_domain) >1:

                            caligned = contactMap.translate_chainMap_to_RefSeq(chain, pdbClass)
                            if hasattr(caligned, 'refseq'):
                                value = True
                            else:
                                value = False
                                
                            list_of_uniprotids.append(uniprot_ID)

                            dom_names = []
                            SH2_dom = []
                            other_dom = []
                            for i in parse_domain:
                                name, IPR, ranges = i.split(':')
                                dom_names.append(name)
                                if DOMAIN in name:
                                    SH2_dom.append(i)
                                else:
                                    other_dom.append(i)

    #                           'SH2 = 1 and other domains >=1'
                            if len(SH2_dom) == 1 and len(other_dom) >= 1:

                                header_1, IPR, ROI_1 = SH2_dom[0].split(':')
                                ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')
                                for val in other_dom:
                                    header_2, IPR, ROI_2 = val.split(':')
                                    ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
    #                                     fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                    fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
                                    feature_header = header_2
                                    caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                     fasta_header, feature_header,filename, append=True, 
                                                                       use_ref_seq_aligned=value)

    #                           'SH2 > 1 and other domains = 0'
                            if len(SH2_dom) > 1 and len(other_dom) == 0:

                                header_1, IPR,ROI_1 = SH2_dom[0].split(':')
                                header_2, IPR,ROI_2 = SH2_dom[1].split(':')
                                ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')
                                ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
                                fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
    #                                 fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                feature_header = header_2
                                caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                 fasta_header, feature_header,filename, append=True, 
                                                                    use_ref_seq_aligned=value)

                                header_1, IPR, ROI_1 = SH2_dom[1].split(':')
                                header_2, IPR, ROI_2 = SH2_dom[0].split(':')
                                ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')
                                ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
                                fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
    #                                 fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                feature_header = header_2
                                caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                 fasta_header, feature_header,filename, append=True, 
                                                                    use_ref_seq_aligned=value)

    #                           'SH2 > 1 and other domains > 0'
                            if len(SH2_dom) > 1 and len(other_dom) > 0:

                                for val1 in SH2_dom:
                                    header_1, IPR, ROI_1 = val1.split(':')
                                    ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')

                                    for val2 in other_dom:
                                        header_2, IPR, ROI_2 = val2.split(':')
                                        ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
                                        fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
    #                                         fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                        feature_header = header_2
                                        caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                 fasta_header, feature_header,filename, append=True, 
                                                                    use_ref_seq_aligned=value)

                                header_1, IPR, ROI_1 = SH2_dom[0].split(':')
                                header_2, IPR, ROI_2 = SH2_dom[1].split(':')
                                ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')
                                ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
                                fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
    #                                 fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                feature_header = header_2
                                caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                 fasta_header, feature_header,filename, append=True, 
                                                                    use_ref_seq_aligned=value)

                                header_1, IPR, ROI_1 = SH2_dom[1].split(':')
                                header_2, IPR, ROI_2 = SH2_dom[0].split(':')
                                ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')
                                ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
                                fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
    #                                 fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                feature_header = header_2
                                caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                 fasta_header, feature_header,filename, append=True, use_ref_seq_aligned=value)

    if append_refSeq:
        inputfile = filename+'.fasta'
        with open(inputfile, 'a') as file:
            for query_uniprotid in set(list_of_uniprotids):
                fasta_seq = SeqIO.parse(open(reference_fastafile), 'fasta')
                for fasta in fasta_seq:
                    name, sequence = fasta.id, str(fasta.seq)
                    ref_uniprot_id, ref_gene,ref_domain, ref_index, ref_ipr, ref_start, ref_end = name.split('|')
                    if query_uniprotid == ref_uniprot_id:
                        file.write('>'+name+'\n'+sequence+'\n')


def make_mergedFeatureFiles(fasta_unaligned,fasta_aligned,feaFile_unaligned,feaFile_aligned,
                            feaFile_merge_aligned,feaFile_merge_unaligned, interface='NonCanonical'):
    '''
        For a given fasta and feature file with sequences extracted from structures, we can merge the features across several structures and project onto the reference sequence. The fasta files generated using 'NonCanonicalFeatures' and 'CanonicalFeatures' functions will include both the structure and reference sequences to be able to merge the features based of the reference sequence alignments. 
        Parameters
        ----------
           fasta_unaligned : str
               location of the input fasta file with unaligned sequences
            fasta_aligned : str
                location of the input fasta file with aligned sequences (can be created using any alignment software such as MAFFT, etc.)
            feaFile_unaligned : str
                location of the input feature file with features that can be projected onto the input fasta files
            feaFile_aligned : str
                location of the feature file with features where the residue positions are translated from unaligned to aligned sequence positions. This is a temporary file created while merging the features. To generate this file, one can use 'makeFeatureFile_updateSeqPos' to generate this file.
            feaFile_merge_aligned : str
                location of feature file that contains the merge features and the residue positions are with respect to the aligned sequence numbering. This is also a temporary file. 
            interface : str
                select the interface of interest - NonCanonical (default) or Canonical

        Returns
        -------
            feaFile_merge_unaligned : str
                the location of the final feature file that is generated here with merged features. The feature positions are with respect to the unaligned sequences. These features can be visualized using the input fasta file (aligned or unaligned could be used for Jalview purpose)
                '''
    
    makeFeatureFile_updateSeqPos(fasta_unaligned, fasta_aligned, feaFile_unaligned, feaFile_aligned)
    mergedFeatures(fasta_unaligned, fasta_aligned, feaFile_aligned, feaFile_merge_aligned, 
                       alignment_similarity = 85, feature_cutoff = 30, interface=interface)
    makeFeatureFile_updateSeqPos(fasta_aligned, fasta_unaligned, feaFile_merge_aligned, feaFile_merge_unaligned)
    
    if os.path.exists(feaFile_aligned):
        os.remove(feaFile_aligned)
    else:
        print("The file does not exist")
        
    if os.path.exists(feaFile_merge_aligned):
        os.remove(feaFile_merge_aligned)
    else:
        print("The file does not exist")



                        
def makeHeader(PDB_ID, entity_id, ROI_start, ROI_end, domain_of_interest, pdb_ann_file, reference_fastafile):
    '''makes a fasta header to include all the fields present in reference fasta header for a specific uniprot ID.
    
    Parameters
    ----------
        PDB_ID : str
        entity_id : int
            entity of the doamin of interest
        ROI_start : int
            starting residue of domain of interest
        ROI_end : int
            last residue of the domain of interest
        domain_of_interest : str
        pdb_ann_file : str
            PDB reference file with all PDB structures annotated and filtered based on the domain of interest
        reference_fastafile : str
            fasta file with reference sequences of domain of interest obtained from the Uniprot reference csv file
            
    Returns
    -------
        returns a reference fasta header for a specific PDB ID and domain of interest in that structure'''
    
    df = pd.read_csv(pdb_ann_file)
    domain = (df.loc[(df['PDB_ID'] == PDB_ID) & (df['ENTITY_ID'] == entity_id), ['ref:domains']] ).values.item()
    uniprot_id = (df.loc[(df['PDB_ID'] == PDB_ID) & (df['ENTITY_ID'] == entity_id), ['database_accession']] ).values.item()
    gene = (df.loc[(df['PDB_ID'] == PDB_ID) & (df['ENTITY_ID'] == entity_id), ['ref:gene name']] ).values.item()
    domainlist = domain.split(';')
    DOI = []
    for i in range(len(domainlist)):
        domain_name, IPR, ranges = domainlist[i].split(':')
        
        if domain_of_interest in domain_name:
            DOI.append(domainlist[i])
        
    for j in range(len(DOI)):
        DOI_domain, DOI_IPR, DOI_range = DOI[j].split(':')
        start, end, gap, mutations = DOI_range.split(',')
        
        if int(start) == ROI_start and int(end) == ROI_end:
            index = j+1
            header = uniprot_id + '|'+gene+'|'+DOI_domain+'|'+str(index)+'|'+DOI_IPR+'|'+start+'|'+end

    fasta_seq = SeqIO.parse(open(reference_fastafile), 'fasta')
    for fasta in fasta_seq:
        name, sequence = fasta.id, str(fasta.seq)
        ref_uniprot_id, ref_gene,ref_domain, ref_index, ref_ipr, ref_start, ref_end = name.split('|')
        if name == header:
            fasta_header = header

        else:
            if ref_uniprot_id == uniprot_id and ref_gene == gene:
                if ROI_start in range(int(ref_start)-5, int(ref_start)+5) and ROI_end in range(int(ref_end)-5, int(ref_end)+5):
#                 if ROI_start == int(ref_start) and ROI_end == int(ref_end):
                    fasta_header = name

    return fasta_header                                
                                
                                
def identityScore(aligned_sequences_list):
    '''finds a similarity score percent for a group of sequences '''
    
    score = 0
    len_sequence = len(aligned_sequences_list[0])
    for seq_len in range(0,len_sequence):
        tmp = []
        for i in range(0,len(aligned_sequences_list)):
            tmp.append(aligned_sequences_list[i][seq_len])
        if len(set(tmp)) == 1:
            score += 1
    percent = (score/len_sequence)*100
    return percent

def reference_seq(gene_of_interest, domain_of_interest, uniprot_ref_file):
    '''generates a dictionary with keys and values as fasta headers and fasta sequences extracted from data stored in the Uniprot reference file for domain of interest'''
    list_sequence = []
    list_domain = []
    for name, group in df.groupby('Gene'):
        for index, row in group.iterrows():
            if name == gene_of_interest:
                interpro_domain = row['Interpro Domains']
                sequence = row['Ref Sequence']
                uniprot_id = row['UniProt ID']
                parse_interpro_domain = interpro_domain.split(';')

                doi = []
                other = []
                for i in parse_interpro_domain:
                    domain, IPR, (start), (stop) = i.split(':')
                    if domain_of_interest in domain:
                        doi.append(i)
                    else:
                        other.append(i)
                fasta_dict = {}
                index = 1
                for j in doi:

                    domain_1, IPR, (start_1), (stop_1) = j.split(':')
                    for k in other:
                        domain_2 = k.split(':')[0]
                        newheader = uniprot_id +'|'+gene_of_interest+'|'+domain_1+'|'+domain_2+'|'+str(index)
                        domain_sequence = sequence[int(start_1)-1:int(stop_1)]
                        fasta_dict[newheader]=domain_sequence
                        index +=1
    return(fasta_dict)

def assign_ID_AA(sequence_of_domain):
    '''assigns a position value to each residue of the sequence provided as an input. The aligned seqeunces that contain '-' characters will be skipped while reporting the updated residue positions '''
    sequence_with_ID=[]
    sequence_with_ID_upd=[]
    len_of_seq = 1
    for i in sequence_of_domain:
        sequence_with_ID.append(str(i)+"-"+str(len_of_seq))
        len_of_seq +=1
        
    for j in sequence_with_ID:
        split_j = j.split('-')
        if split_j[0] != '':
            sequence_with_ID_upd.append(j)
        
    return(sequence_with_ID_upd, len(sequence_with_ID_upd))

def pair_ref_aln(sequence1, sequence2, length_of_domain):
    '''generates a dictionary with translated residue positions
    Parameters
    ----------
        sequence1 : list
            list generated from 'assign_ID_AA(sequence_of_domain)[0]' - this is for unaligned sequence numbering
        sequence2 : list
            list generated from 'assign_ID_AA(sequence_of_domain)[0]' - this is for the aligned sequence numbering
        length_of_domain : int
            the length of the domain sequence is used to check whether there are any insertions/deletions or differences in the two aligned and unaligned seqeunces. Expecting to get the same sequence whether aligned or unaligned. ''' 
    matrix_AA_ID = {}
    
    for num in range(length_of_domain):
        upd = sequence1[num].split('-')
        ref = sequence2[num].split('-')
        matrix_AA_ID[upd[1]] = ref[1]
        
    return(matrix_AA_ID)

def makeFeatureFile_updateSeqPos(fasta_file, fasta_aln_file, input_featurefile, output_featurefile):
    '''Makes a feature file with feature positions translated to the ones on the aligned sequences.
    Parameters
    ----------
        fasta_file : str
            location of the fasta file with unaligned seqeunces used as input here
        fasta_aln_file : str
            location of teh fasta file with aligned sequences used as input as well (can be aligned by any software)
        input_featurefile : str
            location of the feature file that goes with the input fasta files
    Returns
    -------
        output_featurefile : str
            location to store the output feature file with feature residue positions translated from unaligned to aligned positions. 
            For example: unaligned seq = 'AKPLYYG'; aligned seq = 'AKP--LYY-G'. If 'A' and 'L' are features, the the input feature file would have A:1 and L:4 but the output file here will show A:1 and L:6. 
            We cannot use this output feature file to project onto the fasta sequences on Jalview. But we can use this numbering for other analysis purposes'''
    
    dict_ref_header_seq = {}
    dict_aln_header_seq = {}

    ref_seq = SeqIO.parse(open(fasta_file), 'fasta')
    aln_seq = SeqIO.parse(open(fasta_aln_file), 'fasta')

    for fasta in ref_seq:
        name, sequence = fasta.id, str(fasta.seq)
        dict_ref_header_seq[name] = sequence

    for fasta in aln_seq:
        name, sequence = fasta.id, str(fasta.seq)
        dict_aln_header_seq[name] = sequence

    dict_ref_aln = {}

    for k1, v1 in dict_ref_header_seq.items():
        for k2, v2 in dict_aln_header_seq.items():
            if k1 == k2:
                dict_ref_aln[k1] = (v1, v2)
                break

    for key, value in dict_ref_aln.items():
#         print(key)
        ref = assign_ID_AA(value[0])
        aln = assign_ID_AA(value[1])
        if ref[1] == aln[1]:
            length_of_domain = ref[1]
            matrix = pair_ref_aln(ref[0], aln[0], ref[1])
            with open(output_featurefile, 'a') as file:
                for line in open(input_featurefile,'r'):
                    line.strip() 
                    line = line.split('\t')
                    feature_header = line[0]
                    header = line[1]
                    feature_1 = int(line[3])
                    feature_2 = int(line[4])
                    if key == header:
                        for unaln_val, aln_val in matrix.items():
                            if int(unaln_val) == (feature_1):
#                                 print(key, unaln_val, aln_val)
                                file.write(feature_header+"\t"+str(header)+"\t-1\t"+str(aln_val)+"\t"+str(aln_val)+'\t'+str(feature_header)+'\n')
#     print('Created feature file for aligned fasta sequences!')
    
def mergedFeatures(fasta_unaligned, fasta_aligned, features_for_alignedFasta, output_features, 
                   alignment_similarity = 85, feature_cutoff = 30, interface = 'NonCanonical'):
    '''Collapse features across PDB structures onto the reference sequence.
    Parameters
    ----------
        fasta_unaligned : str
            location of the fasta file with unaligned seqeunces used as input here
        fasta_aln_file : str
            location of teh fasta file with aligned sequences used as input as well (can be aligned by any software)
        features_for_alignedFasta : str
            location of the input feature file that have the trasnlated residue positions with respect to teh aligned sequences. This is generated using 'makeFeatureFile_updateSeqPos'.
        alignment_similarity : int
            while grouping the sequences to merge features across multiple structures, we want to make sure that identical seqeunces are under consideration. So, we use >=85% as the identity score between the sequences to create some flexibility for taking int oaccount the small differences taht arise while structure determination expriments
        feature_cutoff : int
            a feature present in more than the set threshold will be considered and will make it to the final feature set. 
        interface : str
            chose between 'NonCanonical' or 'Canonical' This is mainly to create specific headers in each of the cases. 
    Returns
    -------
        output_features : str
            location of the feature file with merged set of features that can now be visualized using a reference seqeunce file. 
            '''
    
    fasta_seq = SeqIO.parse(open(fasta_unaligned), 'fasta')
    identifier_list = []
    for fasta in fasta_seq:
        name, sequence = fasta.id, str(fasta.seq)
        splitname = name.split('|')
        if len(splitname) > 7:
            uid, gene, dom1, index, IPR, start, end, dom2,pdb = name.split('|')
            if interface == 'NonCanonical':
                identifier = uid+'|'+gene+'|'+dom1+'|'+index+'|'+IPR+'|'+start+'|'+end+'|'+dom2
            if interface == 'Canonical':
                identifier = uid+'|'+gene+'|'+dom1+'|'+index+'|'+IPR+'|'+start+'|'+end+'|lig'
            if identifier not in identifier_list:
                identifier_list.append(identifier)

    alnseq_dict = {}
    for i in identifier_list:
        tmp_list = []
        fasta_seq = SeqIO.parse(open(fasta_aligned), 'fasta')
        for fasta in fasta_seq:
            name, sequence = fasta.id, str(fasta.seq)

            if i in name:
                tmp_list.append(sequence)

        alnseq_dict[i] = tmp_list
        percent = identityScore(tmp_list)

        tmp_features = []
        tmp_headers = []
        if int(percent) >= alignment_similarity:
            for line in open(features_for_alignedFasta,'r'):
                line.strip() 
                line = line.split('\t')
                features = int(line[3])
                header = str(line[1])
                splitname = header.split('|')
                if len(splitname) > 7:
                    uid, gene, dom1, index, IPR, start, end, dom2,pdb = header.split('|')
                    if i in header:
                        tmp_features.append(features)
                        if header not in tmp_headers:
                            tmp_headers.append(header)
                            if interface == 'NonCanonical':
                                feature_header = dom2
                            if interface == 'Canonical':
                                feature_header = 'lig'
                            header_for_reference = uid+'|'+gene+'|'+dom1+'|'+index+'|'+IPR+'|'+start+'|'+end

        tmp_write = []
        with open(output_features,'a') as file:
            for fea in tmp_features:
                c = tmp_features.count(fea)
                fea_percent = 100*(c/len(tmp_headers))
                if fea_percent > feature_cutoff:
                    if fea not in tmp_write:
                        tmp_write.append(fea)
                        file.write(feature_header+'\t'+header_for_reference+'\t-1\t'+str(fea)+'\t'+str(fea)+'\t'+feature_header+'\n')
    print('Created feature file with merged features!')
