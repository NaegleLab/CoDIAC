import pandas as pd
import requests

def appendRefFile(input_RefFile, outputfile):
    """
    Takes a reference file generated by CODAC.py and adds interpro domain metadata as a new column
    Parameters
    ----------
    input_RefFile: string
        name of the input reference file generated from the  makeRefFile function in CODAC.py
    outputfile: string
        name of the file to be outputted by this function
    Returns
    -------
        .csv with name specified in outputfile containing interpro annotated domain metadata
    """
    df = pd.read_csv(input_RefFile)
    uniprotList = df['UniProt ID'].to_list()
    
    domainMD = generateDomainMetadata(uniprotList)
    processed_MD = process_proteins(domainMD)
    metadata_string_list, domain_arch_list = generate_domain_metadata_string_list(processed_MD, uniprotList)
    
    df['Interpro Domains'] = metadata_string_list
    df['Interpro Domain Architecture'] = domain_arch_list
    df.to_csv(outputfile, index=False)
    print('Interpro metadata succesfully incorporated')

def generate_domain_metadata_string_list(processed_metadata, uniprot_list):
    """
    Condenses protein metadata into strings and outputs them as a list corresponding with indices of accessions in uniprot_list
    Parameters
    ----------
    processed_metadata: dict
        dictionary of protein accessions containing domain metadata
    uniprot_list: list
        list of uniprot IDs
    Returns
    -------
    metadata_string_list: list of lists of strings
        list containing lists of domain strings that condense domain metadata information
    """
    metadata_string_list = []
    domain_arch_list = []
    for i in range(len(uniprot_list)):
        string_list = []
        accession = uniprot_list[i]
        domain_metadata_list = processed_metadata[accession]
        for domain_dict in domain_metadata_list:
            short_name = domain_dict['short']
            id = domain_dict['accession']
            start = domain_dict['start']
            end = domain_dict['end']
            metadata_string = short_name+':'+id+':'+str(start)+':'+str(end)
            string_list.append(metadata_string)
        domain_arch = return_domain_architecture(string_list)
        domain_arch_list.append(domain_arch)
        metadata_string_list.append(';'.join(string_list))
    return metadata_string_list, domain_arch_list

def generateDomainMetadata(uniprot_accessions):
    interpro_url = "https://www.ebi.ac.uk/interpro/api"
    domain_database = ["smart", "pfam", "profile", "prosite", "prints", "cdd", "tigrfams", "sfld", "panther"]

    metadata = {}

    for protein_accession in uniprot_accessions:
        
        new_accession_flag = False

        url = interpro_url + "/entry/interpro/protein/uniprot/" + protein_accession
        try:
            resp = requests.get(url).json()
        except:
            metadata[protein_accession] = []  
            continue         
        
        current_accession = [protein_accession]

        entry_list = []

        for entry in resp['results']:

            if entry['metadata']['type'] != 'domain':
                continue
            
            info_dict = {}

            interpro_dict = collect_data(entry, current_accession, domain_database)

            info_dict['interpro'] = interpro_dict

            entry_list.append(info_dict)
        
        # Accessing domain information from other databases 
        for database in domain_database:
            db_url = interpro_url + "/entry/" + database + "/protein/uniprot/" + protein_accession
            
            try:
                db_domain_resp = requests.get(db_url).json()
            except:
                continue

            for entry in db_domain_resp['results']:

                if entry['metadata']['type'] != 'domain':
                    continue

                db_dict = collect_data(entry, current_accession)
                integrated = entry['metadata']['integrated']
                index = 0
                for i in range(len(entry_list)):
                    if 'interpro' in entry_list[i]:
                        if entry_list[i]['interpro']['accession'] == integrated:
                            index = i
                            break
                    if i == len(entry_list)-1:
                        index = i+1

                if index==len(entry_list):
                    # This part is the necessary part
                    info_dict = {}
                    info_dict[database] = db_dict
                    entry_list.append(info_dict)
                else:
                    # add to the already existing dictionary that contains the interpro domain info
                    entry_list[index][database] = db_dict   
        metadata[protein_accession] = entry_list 
    return metadata

# domain_database added to figure out if there are any databases we might be missing
def collect_data(entry, protein_accession, domain_database=None):
    
    # can be removed later
    if domain_database is not None:
        keys = entry['metadata']['member_databases'].keys()

        for item in keys:
            if item not in domain_database:
                print(item)

    # back to normal code
    dictionary = {}
    
    dictionary['name'] = entry['metadata']['name']
    dictionary['accession'] = entry['metadata']['accession']
    # try-except can be removed post testing
    try:
        num_boundaries = len(entry['proteins'][0]['entry_protein_locations'])
        dictionary['num_boundaries'] = num_boundaries
        # list that contains all the discontinous boundaries
        boundaries_list = []
        for i in range(num_boundaries):
            boundaries_dict = {}
            bounds = entry['proteins'][0]['entry_protein_locations'][i]['fragments'][0]
            boundaries_dict['start'] = bounds['start']
            boundaries_dict['end'] = bounds['end']
            boundaries_list.append(boundaries_dict)

            if bounds['dc-status'] != "CONTINUOUS":
                print("Domain not continuous")
                print(protein_accession + " : " + dictionary["accession"])

        dictionary['boundaries'] = boundaries_list
    except:
        print('boundaries dictionary raises error')
    
    return dictionary

# complete
def process_proteins(PROTEINS):
    """
    Takes in a dictionary of proteins accessions containing domain annotations for each protein and outputs a dictionary of proteins with domains for that protein organized and collapsed such that no domain intersects with another domain
    Parameters
    ----------
    PROTEINS: dict
        dictionary of protein accessions containing a list of domains
    Returns
    -------
    processed_proteins_dict: dict
        dictionary of protein accessions containing a list of collapsed and arranged domains
    """
    interpro_dict = collapse_InterPro_Domains(PROTEINS) #inplace deletion of domainMD, may want to move this to process_proteins?

    processed_proteins_dict = {}
    for protein_accession, domains_list in PROTEINS.items():
        # domains_list is a list of unprocessed domains which can be processed with make_interpro_annotations
        processed_domains_list = make_interpro_annotations(domains_list)
        
        # The list of processed domains can be associated with its accession
        processed_proteins_dict[protein_accession] = processed_domains_list
    
    return processed_proteins_dict

# complete
def make_interpro_annotations(domains_list):
    """
    Creates a list of domains arranged by sequence
    Parameters
    ----------
    domains_list : list
        A list of domain dictionaries
    Returns
    -------
    overall_dict : list of dictionaries
        
    """
    interpro_domains_list = []

    for i in range(len(domains_list)):
        if 'interpro' in domains_list[i].keys():
            interpro_domains_list.append(domains_list[i])

    # domains_list is processed such that the domains that intersect can be collapsed based on hierarchy
    collapsed_domains_list, hierarchy_dict = collapse_domains(interpro_domains_list)
    #print(collapsed_domains_list)

    formatted_domains_list = handle_bounds(collapsed_domains_list)
    
    formatted_domains_list.sort(key=sort_function)

    return formatted_domains_list

# documentation incomplete
def handle_bounds(domain_list):
    """
    Handles instances where there may be discontinuous boundaries for a single domain
    Returns
    -------
    processed_domain_list : list
        list of dictionaries containing one or more regions annotated as the same domain
    """ 
    formatted_domain_list = []

    for i in range(len(domain_list)):
        domain_dict = domain_list[i]['interpro']

        num_boundaries = domain_dict['num_boundaries']
        name = domain_dict['name']
        short = domain_dict['short']
        accession = domain_dict['accession']
        bounds = domain_dict['boundaries']

        if num_boundaries > 1:
            for i in range(num_boundaries):
                new_domain_dict = {}
                new_domain_dict['short'] = short
                new_domain_dict['accession'] = accession
                new_domain_dict['start'] = bounds[i]['start']
                new_domain_dict['end'] = bounds[i]['end']
                new_domain_dict['name'] = name
                formatted_domain_list.append(new_domain_dict)
        else:
            new_domain_dict = {}
            new_domain_dict['short'] = short
            new_domain_dict['accession'] = accession
            new_domain_dict['start'] = bounds[0]['start']
            new_domain_dict['end'] = bounds[0]['end']
            new_domain_dict['name'] = name
            formatted_domain_list.append(new_domain_dict)
        
    return formatted_domain_list

# documentation incomplete
def sort_function(formatted_domain_dict):
    return formatted_domain_dict['start']


# complete
def compare_hierarchy(hierarchy_dict):
    """
    Checks the domain hierarchies of all domains present within a protein and collapses domains into the domain at the top of the hierarchy
    Parameters
    ----------
    hierarchy_dict : dict
        Dictionary containing hierarchy information where the key is the interpro accession of the domain and the value is the dictionary with data regarding the domain's hierarchy
    Returns
    -------
    collapsed_domain_mapper : dict
        A dictionary that maps the collapsed domains to the domain they have been collapsed under
    """
    skipped_accessions = []
    collapsed_domain_mapper = {}
    interpro_accession_list = list(hierarchy_dict.keys())
    for i in range(len(interpro_accession_list)):
        accession = interpro_accession_list[i]
        # check if the accession has already been collapsed and skip iteration if collapsed
        if accession in skipped_accessions:
            continue
        
        compare_to_list = interpro_accession_list[i+1:]
        collapsable_domains = []
        collapsable_domain_found = False
        for j in range(len(compare_to_list)):
            compare_to_accession = compare_to_list[j]
            # check if the accession has already been collapsed and skip iteration if collapsed
            if compare_to_accession in skipped_accessions:
                continue
            if hierarchy_dict[accession] == hierarchy_dict[compare_to_accession]:
                collapsable_domain_found = True
                collapsable_domains.append(compare_to_accession)
                skipped_accessions.append(compare_to_accession)

        if collapsable_domain_found:
            collapsable_domains.append(accession)
            # implement breadth first search to find the topmost domain within the hierarchy
            topmost_accession = BFS(collapsable_domains, hierarchy_dict[accession])
            collapsable_domains.remove(topmost_accession)
            collapsed_domain_mapper[topmost_accession] = collapsable_domains
    
    return collapsed_domain_mapper


# complete            
def BFS(accession_list, hierarchy_tree):
    """
    A function that performs breadth-first search on a provided hierarchy tree
    Parameters
    ----------
    accession_list : list
        list of accessions that the algorith is searching for within the tree
    hierarchy_tree : dict
        a dictionary that contains the hierarchy data
    Returns
    -------
    topmost_accession : str
        the accession found at the topmost level of the tree that also exists in accession_list
    """
    queue = [hierarchy_tree]
    while len(queue) != 0:
        tree = queue.pop(0)
        if tree['accession'] in accession_list:
            topmost_accession = tree['accession']
            return topmost_accession
        else:
            children = tree['children']
            for i in range(len(children)):
                queue.append(children[i])


# complete
def collapse_domains(domain_list):
    """
    Takes the domains within an interpro protein that have similar/intersecting boundaries
    and collapses them based on domain hierarchy
    Parameters
    ----------
    domain_list : list
        A list of domain dictionaries
    Returns
    -------
    collapsed_domain_list : list
        A list of dictionaries where domains at the top of the hierarchy table are kept and any underlying domains have been removed
    """    
    hierarchy_dict = {}
    interpro_accessions_list = []
    interpro_url = 'https://www.ebi.ac.uk/interpro/api'

    for i in range(len(domain_list)):
        interpro_accessions_list.append(domain_list[i]['interpro']['accession'])
    
    for i in range(len(interpro_accessions_list)):
        interpro_accession = interpro_accessions_list[i]
        url = interpro_url + "/entry/interpro/" + interpro_accession
        resp = requests.get(url).json()

        domain_list[i]['interpro']['short'] = resp['metadata']['name']['short']
        hierarchy_dict[interpro_accession] = resp['metadata']['hierarchy']
    
    collapsed_domain_mapper = compare_hierarchy(hierarchy_dict)
    
    to_remove = []
    for key, mapped_list in collapsed_domain_mapper.items():
        for i in range(len(domain_list)):
            if domain_list[i]['interpro']['accession'] in mapped_list:
                to_remove.append(i)
    
    to_remove.sort(reverse=True)
    for i in range(len(to_remove)):
        domain_list.pop(to_remove[i])

    return domain_list, hierarchy_dict

def return_domain_architecture(domain_list):
    """
    Given a domain_list, list of domain information short_name:id:start:end return a domain 
    architecture, which is the | separated list of domain names, in the order they appear in protein

    """
    #Domain Architecture
    domdict = {}
    for domain_info in domain_list:
        name, id, start, end = domain_info.split(':')
        start = int(start)
        domdict[start] = end, name

    sorted_dict = dict(sorted(domdict.items(),reverse=False))
    domain_arch = []
    for key, value in sorted_dict.items():
        #sort_start = key
        #sort_end = value[0]
        domain = value[1]
        domain_arch.append(domain)

    
    final_domarch = '|'.join(domain_arch)   
    return final_domarch

def calculate_domain_overlap(interpro, overlap_threshold, length_larger_than):
    """
    Given a dict of intpro domains, where keys are the position they are found in the list for domain metadata, 
    calculate the percent overlap between all domains -- this skips multiple entries of the same domain from the calculation
    Returns a list of domain entries to delete, because they are contained at > threshold in another domain.
    """

    #given a number of possible domains, find the domains that overlap each other by at least 50% of the smallest domain
    keys = interpro.keys()
    to_remove = []
    for i in range(0, len(keys)):
        boundaries_i = interpro[i]['boundaries']
        if len(boundaries_i) > 1:
            continue

        for j in range(i+1, len(keys)):
            boundaries_j = interpro[j]['boundaries']
            if len(boundaries_j)> 1:
                    #print("Skipping double domains")
                    continue
            else:
                set_i = set(range(boundaries_i[0]['start'], boundaries_i[0]['end']))
                set_j  = set(range(boundaries_j[0]['start'], boundaries_j[0]['end']))
                if len(set_i) < len(set_j):
                    if length_larger_than*len(set_i) < len(set_j):
                        intersection = len(set_i.intersection(set_j))/len(set_i)
                        if intersection > overlap_threshold:
                            #print("Should replace domain %d with %d"%(i, j))
                            to_remove.append(i)
                else:
                    if length_larger_than*len(set_j) < len(set_i):
                        intersection = len(set_j.intersection(set_i))/len(set_j)
                        if intersection > overlap_threshold:
                            #print("Should replace domain %d with %d"%(j, i))
                            to_remove.append(j)
                #temp_dict[j] = intersection
                #overlap[i] = temp_dict
    return list(set(to_remove))

def collapse_InterPro_Domains(domain_metadata, overlap_threshold=0.8, length_larger_than = 1.4):
    """
    Given a domain metadata dictionary (keys are protein ID and points to list of domain dicts)
    Walk through and first find any larger protein domains that cover smaller interpro domains. Keep 
    the largets of the domains only if they overlap by overlap_threshold or more and the larger protein 
    is length_larger_than x longer (for example 1.4 is 40% longer)
    """

    #first assemble all the interpro domains across the list. Keep track of their entry number, and the information
    interpro_dict = {}
    database_key = 'interpro'
    for protein in domain_metadata:
        domain_list = domain_metadata[protein]
        for i in range(0, len(domain_list)):
            domain_dict_temp = domain_list[i]
            if database_key in domain_list[i]:
                interpro_dict[i] = domain_list[i][database_key]
                #interpro_dict_list[-1]['name'] = 'test' #test - this can modify dict, it's pointer to dictionary
        indexes_to_remove = calculate_domain_overlap(interpro_dict, overlap_threshold, length_larger_than)
        for index in indexes_to_remove:
            del domain_list[index][database_key]
    return interpro_dict
                
            