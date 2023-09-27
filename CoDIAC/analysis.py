from CoDIAC import PDBHelper, contactMap
from CoDIAC import pTyrLigand_helpers as pTyr_helpers
import pandas as pd
import numpy as np
from Bio import SeqIO

def CanonicalFeatures(pdb_ann_file, ADJFILES_PATH, reference_fastafile, error_structures_list, PTM='PTR', mutation=False, 
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
        if PDB_ID not in error_structures_list:

            for index, row in group.iterrows():
                uniprot_ID = row['database_accession']
                
                if isinstance(row['modifications'], str):
                    transDict = PDBHelper.return_PTM_dict(row['modifications'])
                    for res in transDict:
                        if PTM in transDict[res]:

                            entities = PDBHelper.PDBEntitiesClass(main, PDB_ID)
                            for entity in entities.pdb_dict.keys():
                                domains = entities.pdb_dict[entity].domains
                                for domain_num in domains:
                                    if domain_of_interest in domains[domain_num]:
                                        SH2_entity = entity 
                                        check_mutation = (pd.isnull(main.loc[index, 'ref:variants']))

                                transDict = entities.pdb_dict[entity].transDict
                                for res in transDict:
                                    if PTM in transDict[res]:
                                        lig_entity = entity
                                        
                            if check_mutation != mutation:
                                list_of_uniprotids.append(uniprot_ID)
                                print(name, SH2_entity, lig_entity)
                                pdbClass = entities.pdb_dict[lig_entity]
                                dict_of_lig = contactMap.return_single_chain_dict(main, PDB_ID, ADJFILES_PATH, lig_entity)
                                dict_of_SH2 = contactMap.return_single_chain_dict(main, PDB_ID, ADJFILES_PATH, SH2_entity)

                                cm_aligned = dict_of_lig['cm_aligned']
#                                 SH2_gene = get_gene(pdb_ann_file, SH2_entity, PDB_ID)[0]
#                                 lig_gene = get_gene(pdb_ann_file, lig_entity, PDB_ID)[0]
#                                 uid = get_gene(pdb_ann_file, SH2_entity, PDB_ID)[1]

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
                                                    fasta_header = makeHeader(PDB_ID, SH2_entity,int(SH2_start), int(SH2_stop),domain_of_interest,pdb_ann_file, reference_fastafile)+'|'+PDB_ID


                                                    if lig_entity == SH2_entity:

                                                        cm_aligned.print_fasta_feature_files(SH2_start, SH2_stop, res_start, res_end,
                                                                 fasta_header, 'pTyr',PTM_file, append=True, 
                                                                    use_ref_seq_aligned=True)

                                                        cm_aligned.print_fasta_feature_files(res_start, res_end, SH2_start, SH2_stop,
                                                                 fasta_header, 'SH2',SH2_file, append=True, 
                                                                    use_ref_seq_aligned=True)

                                                    contactMap.print_fasta_feature_files(arr_alt, to_dict['cm_aligned'].refseq, 
                                                        SH2_start, SH2_stop, to_dict['cm_aligned'].return_min_residue(), 
                                                        res_start, res_end, from_dict['cm_aligned'].return_min_residue(),
                                                        fasta_header,'pTyr', PTM_file, threshold=1, append=True)


                                                    contactMap.print_fasta_feature_files(arr, from_dict['cm_aligned'].structSeq, 
                                                                 res_start, res_end, from_dict['cm_aligned'].return_min_residue(), 
                                                                 SH2_start, SH2_stop, to_dict['cm_aligned'].return_min_residue(),
                                                                fasta_header,'SH2', SH2_file, threshold=1, append=True )
                                    
    inputfile = PTM_file+'.fasta'
    with open(inputfile, 'a') as file:
        for query_uniprotid in set(list_of_uniprotids):
            fasta_seq = SeqIO.parse(open(reference_fastafile), 'fasta')
            for fasta in fasta_seq:
                name, sequence = fasta.id, str(fasta.seq)
                ref_uniprot_id, ref_gene,ref_domain, ref_index, ref_ipr, ref_start, ref_end = name.split('|')
                if query_uniprotid == ref_uniprot_id:
                    file.write('>'+name+'\n'+sequence+'\n')
                    

def NonCanonicalFeatures(pdb_ann_file, ADJFILES_PATH, reference_fastafile, error_structures_list, mutation = False, DOMAIN = 'SH2', 
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
                                                                       use_ref_seq_aligned=True)

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
                                                                    use_ref_seq_aligned=True)

                                header_1, IPR, ROI_1 = SH2_dom[1].split(':')
                                header_2, IPR, ROI_2 = SH2_dom[0].split(':')
                                ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')
                                ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
                                fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
    #                                 fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                feature_header = header_2
                                caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                 fasta_header, feature_header,filename, append=True, 
                                                                    use_ref_seq_aligned=True)

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
                                                                    use_ref_seq_aligned=True)

                                header_1, IPR, ROI_1 = SH2_dom[0].split(':')
                                header_2, IPR, ROI_2 = SH2_dom[1].split(':')
                                ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')
                                ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
                                fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
    #                                 fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                feature_header = header_2
                                caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                 fasta_header, feature_header,filename, append=True, 
                                                                    use_ref_seq_aligned=True)

                                header_1, IPR, ROI_1 = SH2_dom[1].split(':')
                                header_2, IPR, ROI_2 = SH2_dom[0].split(':')
                                ROI_10, ROI_11, gap_1, mut_1 = ROI_1.split(',')
                                ROI_20, ROI_21, gap_2, mut_2 = ROI_2.split(',')
                                fasta_header = makeHeader(PDB_ID, entity_id,int(ROI_10), int(ROI_11),DOMAIN,pdb_ann_file, reference_fastafile)+'|'+header_2+'|'+PDB_ID
    #                                 fasta_header = uniprot_ID+'|'+gene+'|'+header_1+'|'+header_2+'|'+PDB_ID
                                feature_header = header_2
                                caligned.print_fasta_feature_files(int(ROI_10), int(ROI_11), int(ROI_20), int(ROI_21),
                                                                 fasta_header, feature_header,filename, append=True, use_ref_seq_aligned=True)


    inputfile = filename+'.fasta'
    with open(inputfile, 'a') as file:
        for query_uniprotid in set(list_of_uniprotids):
            fasta_seq = SeqIO.parse(open(reference_fastafile), 'fasta')
            for fasta in fasta_seq:
                name, sequence = fasta.id, str(fasta.seq)
                ref_uniprot_id, ref_gene,ref_domain, ref_index, ref_ipr, ref_start, ref_end = name.split('|')
                if query_uniprotid == ref_uniprot_id:
                    file.write('>'+name+'\n'+sequence+'\n')

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
            print(domain_name)
            DOI.append(domainlist[i])
            print(domainlist[i], i)
        
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
            fasta_header = header+'|'+PDB_ID

        else:
            if ref_uniprot_id == uniprot_id and ref_gene == gene:
                if ROI_start in range(int(ref_start)-5, int(ref_start)+5) and ROI_end in range(int(ref_end)-5, int(ref_end)+5):
#                 if ROI_start == int(ref_start) and ROI_end == int(ref_end):
                    fasta_header = name+'|'+PDB_ID

    return fasta_header                                
                                
                                
def similarityScore(aligned_sequences_list):
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
    
    list_sequence = []
    list_domain = []
    for name, group in df.groupby('Gene'):
        for index, row in group.iterrows():
            if name == gene_of_interest:
                interpro_domain = row['Interpro Domains']
                sequence = row['Ref Sequence']
                uniprot_id = row['UniProt ID']
                parse_interpro_domain = interpro_domain.split(';')

                sh2 = []
                other = []
                for i in parse_interpro_domain:
                    domain, IPR, (start), (stop) = i.split(':')
                    if domain_of_interest in domain:
                        sh2.append(i)
                    else:
                        other.append(i)
                fasta_dict = {}
                index = 1
                for j in sh2:

                    domain_1, IPR, (start_1), (stop_1) = j.split(':')
                    for k in other:
                        domain_2 = k.split(':')[0]
                        newheader = uniprot_id +'|'+gene_of_interest+'|'+domain_1+'|'+domain_2+'|'+str(index)
                        domain_sequence = sequence[int(start_1)-1:int(stop_1)]
                        fasta_dict[newheader]=domain_sequence
                        index +=1
    return(fasta_dict)

def assign_ID_AA(sequence_of_domain):
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
    
    matrix_AA_ID = {}
    
    for num in range(length_of_domain):
        upd = sequence1[num].split('-')
        ref = sequence2[num].split('-')
        matrix_AA_ID[upd[1]] = ref[1]
        
    return(matrix_AA_ID)

def makeFeatureFile_alignedSeq(fasta_file, fasta_aln_file, input_featurefile, output_featurefile):
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
    print('Created feature file for aligned fasta sequences!')
    
def mergedFeatures(fasta_unaligned, fasta_aligned, features_for_alignedFasta, output_features, 
                   alignment_similarity = 85, feature_cutoff = 30):
    
    fasta_seq = SeqIO.parse(open(fasta_unaligned), 'fasta')
    identifier_list = []
    for fasta in fasta_seq:
        name, sequence = fasta.id, str(fasta.seq)
        uid, gene, dom1, dom2, pdb = name.split('|')
        identifier = gene+'|'+dom1+'|'+dom2
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
        percent = similarityScore(tmp_list)

        tmp_features = []
        tmp_headers = []
        if int(percent) >= alignment_similarity:
            for line in open(features_for_alignedFasta,'r'):
                line.strip() 
                line = line.split('\t')
                features = int(line[3])
                header = str(line[1])
                uniprot, gene, domain_1, domain_2, PDB = header.split('|')
                search_header = gene+'|'+domain_1+'|'+domain_2
                if i in search_header:
                    tmp_features.append(features)
                    if header not in tmp_headers:
                        tmp_headers.append(header)
                        feature_header = domain_2
                        header_with_uniprot = uniprot+'|'+gene+'|'+domain_1+'|'+domain_2

    #     print(i, tmp_headers, feature_header)
        tmp_write = []
        with open(output_features,'a') as file:
            for fea in tmp_features:
                c = tmp_features.count(fea)
                fea_percent = 100*(c/len(tmp_headers))
                if fea_percent > feature_cutoff:
                    if fea not in tmp_write:
                        tmp_write.append(fea)
                        file.write(feature_header+'\t'+header_with_uniprot+'\t-1\t'+str(fea)+'\t'+str(fea)+'\t'+feature_header+'\n')
    print('Created feature file with merged features for aligned seq positions!')