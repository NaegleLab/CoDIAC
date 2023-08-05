import requests
import logging
import os
import pandas as pd
import copy


def fetch_uniprotids(interpro_ID, REVIEWED=True, species='Homo sapiens'):
    """
    Given an InterPro_ID, fetch all the records (or only reviewed records)
    for all species or for a specific taxonomy. 

    Parameters
    ----------
    interpro_ID: str
        InterPro ID to search for
    REVIEWED: bool
        If TRUE, only reviewed records will be returned 
    species: string
        Using scientific name under the Uniprot taxonomy to define species. 
        See here for taxonomy names: https://www.uniprot.org/taxonomy?query=*

    Returns
    -------
    uniprot_ID_list: list
        list of all uniprot IDs found in search. If a species was set, all uniprot IDs for a species
        will be returned in this list, otherwise, all species from search will be returned.
    species_dict: dict
        Dictionary, with top keys equal to the species scientific name and 
        points to an inner dict that keeps track of the database source 'reviewed' or 'unreviewed'
        and has lists of the uniprot IDs found for that species under that database source.

    Examples:
    fetch_uniprotids('IPR000980', REVIEWED=True, species='Homo sapiens')
    will fetch uniprot ids for human only that are reviewed

    fetch_uniprotids('IPR000980', REVIEWED=False, species='all')
    Will fetch all uniprot ids with any record of that interpro id
    """

    count_data = 0
    #species = species.replace(" ", "+")
    interpro_url = "https://www.ebi.ac.uk/interpro/api"
    if REVIEWED:
        url = interpro_url +  f"/protein/reviewed/entry/interpro/{interpro_ID}/" 
    else:
        url = interpro_url + f"/protein/UniProt/entry/interpro/{interpro_ID}/" 
    
    if species!='all' or species!='All' or species!='ALL':
        species_temp = species.replace(" ", "+")
        url = url + f"?search={species_temp}&page_size=200"
    else:
        url = url+ f"?&page_size=200"
        if not REVIEWED:
            print("WARNING: About to search for all records for all species, this will take a while...")
    logging.basicConfig(filename='error.log',
                    format='%(asctime)s %(message)s',
                    encoding='utf-8',
                    level=logging.WARNING)
    
    fetch_all_results = []

    try:
        response = requests.get(url)
        data = response.json()
        count_data = data['count']

        fetch_all_results = fetch_all_results + data['results']

        while data['next'] is not None:

            print("Found next page and downloading", data['next'])
            response = requests.get(data['next'])
            data = response.json()
            fetch_all_results = fetch_all_results + data['results']

    except requests.exceptions.RequestException as err:
        print('Found request exceptions...')
        print('Most likely error is species formatting, check %s'%(species))
        logging.warning(err)

    UNIPROT_ID_LIST = []
    uniprot_dict = {}
    species_dict = {}

    for entry in fetch_all_results:
        Uniprot_Accession = (entry['metadata']['accession'])
        source_database = (entry['metadata']['source_database'])
        taxId = (entry['metadata']['source_organism']['taxId'])
        Scientific_name = entry['metadata']['source_organism']['scientificName']
        fullname = (entry['metadata']['source_organism']['fullName'] )
        if species!='all' or species!='All' or species!='ALL':
            if Scientific_name==species:
                UNIPROT_ID_LIST.append(Uniprot_Accession)
        else:
            UNIPROT_ID_LIST.append(Uniprot_Accession)
        uniprot_dict[Uniprot_Accession] = Scientific_name
        if Scientific_name not in species_dict:
            species_dict[Scientific_name] = {}
        if source_database not in species_dict[Scientific_name]:
            species_dict[Scientific_name][source_database] = []
        species_dict[Scientific_name][source_database].append(Uniprot_Accession)
    print("Fetched %d Uniprot IDs linked to %s, where count expected to be %d"%(len(UNIPROT_ID_LIST),interpro_ID, count_data))
    return(UNIPROT_ID_LIST, species_dict)


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
    print("Fetching domains..")
    domainMD = generateDomainMetadata(uniprotList)
    print("Processing domains...")
    processed_MD = process_proteins(domainMD)
    metadata_string_list, domain_arch_list = generate_domain_metadata_string_list(processed_MD, uniprotList)
    
    df['Interpro Domains'] = metadata_string_list
    df['Interpro Domain Architecture'] = domain_arch_list
    df.to_csv(outputfile, index=False)
    print('Interpro metadata succesfully incorporated')
    return df

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
def process_proteins(PROTEINS, VERBOSE=True, log_file = 'process_interpro_log.txt'):
    """
    Takes in a dictionary of proteins accessions containing domain annotations for each protein and outputs a dictionary of proteins with domains for that protein organized and collapsed such that no domain intersects with another domain
    Parameters
    ----------
    PROTEINS: dict
        dictionary of protein accessions containing a list of domains
    VERBOSE: boolean
        If VERBOSE, print logging statements to log_file output about decisions being made at the coherence of domains and boundaries
    log_file: str
        Location of verbose output. If no file provided, this will default to print logging statements as append to a default file
    Returns
    -------
    processed_proteins_dict: dict
        dictionary of protein accessions containing a list of collapsed and arranged domains
    """

    domainMD_new = {}
    if VERBOSE:
        with open(log_file, "+a") as log:
            log.write("Collapsing domains based on mulit-domains within a larger domain\n")
            log.write("Parameters: \n")
    for protein in PROTEINS:
        #print("Working on %s"%(protein))
        if VERBOSE:
            with open(log_file, "+a") as log:
                log.write(protein+"\n")
        domain_list = PROTEINS[protein]
        domain_list_reduced = collapse_InterPro_Domains(domain_list, VERBOSE, log_file) 
        domainMD_new[protein] = domain_list_reduced
    #interpro_dict = collapse_InterPro_Domains(PROTEINS) #inplace deletion of domainMD, may want to move this to process_proteins?

    processed_proteins_dict = {}
    for protein_accession, domains_list in domainMD_new.items():
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

def return_interpro_domain_list(domain_list):
    """
    From a list of domain dicts that comes from get meta data, get just the subset of domains that are interpro
    """
    interpro_domain_list = []
    for i in range(len(domain_list)):
        if 'interpro' in domain_list[i].keys():
            interpro_domain_list.append(domain_list[i])
    return interpro_domain_list

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

        num_boundaries = len(domain_dict['boundaries'])#domain_dict['num_boundaries']
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
    Checks the domain hierarchies of all domains present within a protein and collapses domains into the domain at the 
    top of the hierarchy
    Parameters
    ----------
    hierarchy_dict : dict
        Dictionary containing hierarchy information where the key is the interpro accession of the domain and 
        the value is the dictionary with data regarding the domain's hierarchy
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


def collapse_InterPro_Domains(domain_list, VERBOSE, log_file, overlap_threshold=0.7):
    """
    Given a domain metadata dictionary (keys are protein ID and points to list of domain dicts), remove
    all domains that overlap a domain by more than the fraction defined by overlap_threshold. Behavior: 
    Keeps the first instance of the domain (as returned by Interpro) and removes later instances, which are 
    less generalized/more specific to proteins.

    Parameters
    ----------
    domain_list: list
        List of domain dictionaries from InterPro for a single protein 
    VERBOSE: boolean
        If VERBOSE, print logging statements to log_file output about decisions being made at the coherence of domains and boundaries
    log_file: str
        Location of verbose output. If no file provided, this will default to print logging statements as append to a default file
    overlap_threshold: float
        Percent of the smaller domain that must overlap with a larger domain before considering collapsing these. 
  
    Returns
    -------
    domain_list_new: list
        List of domain dictionaries with collapsed entries removed. Collapsed when more than one domain is found within a larger
        domain that has > overlap_threshold and is smaller by larger_than_length.
    """

    #first assemble all the interpro domains across the list. Keep track of their entry number, and the information
    #let's make a new version of this, as a dict, where the index values cannot be changed as a result of 
    #deleting values as we go. We'll cast back to a list for export
    domain_list_temp = copy.deepcopy(domain_list)
    domain_dict = {}
    for i in range(0, len(domain_list_temp)):
        domain_dict[i] = domain_list_temp[i].copy()

    interpro_dict = {}
    database_key = 'interpro'
    #for protein in domain_metadata:
        #print("Working on %s"%(protein))
        #domain_list = domain_metadata[protein]
    for key in domain_dict.keys():
        if database_key in domain_dict[key]:
            interpro_dict[key] = domain_dict[key][database_key]

    keys_to_remove = return_keys_to_remove_by_overlap(interpro_dict, overlap_threshold)

    domainVal_report = {}
    if keys_to_remove is not None:
        for key in keys_to_remove:#.sort(reverse=True):
            if key in domain_dict:
                boundary_starts = keys_to_remove[key]
                if len(boundary_starts)==len(domain_dict[key]['interpro']['boundaries']):
                    del domain_dict[key]#[database_key] #can delete the entire entry
                    domainVal_report[key] = 'all instances'
                elif len(boundary_starts)>len(domain_dict[key]['interpro']['boundaries']):
                    print("ERROR: boundaries to remove is larger than available")
                    print(boundary_starts)
                else:
                    #sort the boundaries in reverse order
                    domainVal_report[key] = 'some instances:'
                    for start in boundary_starts:
                        #find the correct value in the range to remove
                        for i in range(0, len(domain_dict[key]['interpro']['boundaries'])):
                            if start == domain_dict[key]['interpro']['boundaries'][i]['start']:
                                domainVal_report[key] += " %d"%(start)
                                #print("DEBUG: removing %d position with start %d, matching %d"%(i, start, domain_dict[key]['interpro']['boundaries'][i]['start']))
                                domain_dict[key]['interpro']['boundaries'].pop(i)
                                break
                    domain_dict[key]['interpro']['num_boundaries'] = len(domain_dict[key]['interpro']['boundaries'])
            else:
                print("ERROR: Index %d no longer in domain list"%(key))
        #re list the domains
    domain_list_new = []
    for key in domain_dict:
        domain_list_new.append(domain_dict[key])

    #write the report if VERBOSE
    if VERBOSE:
        with open(log_file, "+a") as log:
            log.write("Removal based on domain overlap\n")
            for domain_val in keys_to_remove:
                log.write("\t REMOVED: %s, %s\n"%(return_info_string(interpro_dict, domain_val), domainVal_report[domain_val]))
    return domain_list_new


def return_keys_to_remove_by_overlap(interpro_dict, overlap_threshold):
        """
        
        Returns
        -------
        keys_to_remove: dict
            Dictionary with key equal to the dictionary entry of an interpro domain to be removed
            and values equal to a list of unique start positions of boundaries to be removed
            from the domain set
        """

        keys = interpro_dict.keys()
        # keep track of the start position of the domain number and the boundary number that need to be removed
        to_remove = {}
        for i in range(0, len(keys)):
                boundaries_i = interpro_dict[i]['boundaries']
                for j in range(i+1, len(keys)): 
                        #j is always a later set than i, so if it overlaps, remove it
                        boundaries_j = interpro_dict[j]['boundaries']
                        for boundary_j in boundaries_j:
                                for boundary_i in boundaries_i:
                                        #print('DEBUG comparing x to y')
                                        #print(boundary_i)
                                        #print(boundary_j)
                                        intersection = return_overlap(boundary_i, boundary_j)
                                        if intersection >= overlap_threshold:
                                                if j not in to_remove:
                                                        to_remove[j] = []
                                                to_remove[j].append(boundary_j['start'])
        #clean up and make boundaries to remove a set
        for key in to_remove:
                to_remove[key] = list(set(to_remove[key]))
        return to_remove               

                
def return_info_string(interpro_dict, dict_item_number):
    """
    Given an intepro_dict and the item number, return a string of info about the domain, specifically 
    the short name, interpro number, and start, stop 
    """
    return ':'.join([interpro_dict[dict_item_number]['name'], interpro_dict[dict_item_number]['accession'], str(interpro_dict[dict_item_number]['boundaries'])])


def return_overlap(boundary_set_1, boundary_set_2):
    """ 
    Given a boundary dictionary with 'start' and 'end' return the overlap fraction of the 
    set (measuring the intersection of the two boundaries and normalizing by the smaller domain)

    """
    set_1 = set(range(boundary_set_1['start'], boundary_set_1['end']))
    set_2 = set(range(boundary_set_2['start'], boundary_set_2['end']))
    intersection = len(set_1.intersection(set_2))/min(len(set_1), len(set_2))
    return intersection






def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument("InterPro_ID",help="InterPro ID you wish to enter", type=str)
    parser.add_argument("Output_Dir", help="Directory to write the file with Uniprot IDs", type=str)

    args=parser.parse_args()

    df = pd.DataFrame()
    df['Uniprot ID'] = fetch_uniprotids(args.InterPro_ID)
    PATH = args.Output_Dir+'/UniProt_IDs.csv'
    df.to_csv(PATH, index=False)

if __name__ =='__main__':
    Main()