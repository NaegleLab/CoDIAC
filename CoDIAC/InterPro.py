from collections import defaultdict
import requests
import logging
import pandas as pd

import timeit


def fetch_uniprotids(interpro_ID, REVIEWED=True, species='Homo sapiens'):
    """
    Given an InterPro_ID, fetch all the records (or only reviewed records)
    for all species or for a specific taxonomy. 

    Examples:
    Examples, use this module like this, the first being a more restrictive human with reviewed, versus all 
    records associated within an Interpro ID:
    .. code-block:: python
        fetch_uniprotids('IPR000980', REVIEWED=True, species='Homo sapiens') # human proteins with reviewed records
        fetch_uniprotids('IPR000980', REVIEWED=False, species='all') #all species records, reviewed and unreviewed
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
    """
    count_data = 0 # count of records expected to be found
    #species = species.replace(" ", "+")
    interpro_url = "https://www.ebi.ac.uk/interpro/api"
    if REVIEWED: # if reviewed, we need to change the URL
        url = ''.join([interpro_url, "/protein/reviewed/entry/interpro/", interpro_ID, "/"])
    else:
        url = ''.join([interpro_url, "/protein/UniProt/entry/interpro/", interpro_ID, "/"])
    # if species is defined, we need to add this to the URL
    if species.lower() !='all':
        species_temp = species.replace(" ", "+")
        url = ''.join([url, "?search=", species_temp, "&page_size=200"])
    else: # if all species, we need to add a page size to the URL
        url = ''.join([url, "?&page_size=200"])
        if not REVIEWED:
            print("WARNING: About to search for all records for all species, this will take a while...")
    logging.basicConfig(filename='error.log',
                    format='%(asctime)s %(message)s',
                    encoding='utf-8',
                    level=logging.WARNING) # set up logging for requests exceptions
    fetch_all_results = []  # list to store all results

    try:
        with requests.Session() as session:
            response = session.get(url)  # make the request
            data = response.json()  # convert to json
            count_data = data['count']

            fetch_all_results.extend(data['results'])  # use extend instead of +

        # if there are more pages, we need to fetch these as well
        while data['next'] is not None:
            print("Found next page and downloading", data['next'])
            response = session.get(data['next'])
            data = response.json()
            fetch_all_results.extend(data['results'])  # use extend instead of +

    except requests.exceptions.RequestException as err:
        print('Found request exceptions...')
        print('Most likely error is species formatting, check %s' % species)
        logging.warning(err)

    UNIPROT_ID_LIST = []  # list to store all uniprot IDs
    uniprot_dict = {}  # dictionary to store uniprot IDs and their associated species
    species_dict = defaultdict(lambda: defaultdict(list))  # dictionary to store species and their associated uniprot IDs

    all_species = {'all', 'All', 'ALL'}

    for entry in fetch_all_results:  # loop through all results
        Uniprot_Accession = entry['metadata']['accession']
        source_database = entry['metadata']['source_database']
        Scientific_name = entry['metadata']['source_organism']['scientificName']

        if species not in all_species and species not in Scientific_name:
            continue

        UNIPROT_ID_LIST.append(Uniprot_Accession)
        uniprot_dict[Uniprot_Accession] = Scientific_name
        species_dict[Scientific_name][source_database].append(Uniprot_Accession)

    print(f"Fetched {len(UNIPROT_ID_LIST)} Uniprot IDs linked to {interpro_ID}, where count expected to be {count_data}")
    return(UNIPROT_ID_LIST, species_dict)


def collect_data(entry):
    """
    Given an entry from the InterPro API, collect the data for a protein accession
    and return a dictionary with the keys 'name', 'accession', 'num_boundaries', 'boundaries'
    where 'boundaries' is a list of dictionaries with keys 'start' and 'end'.
    If 'short_name' is in the extra fields, this will be added to the dictionary as 'short'
    
    Parameters
    ----------
    entry: dict
        dictionary from the InterPro API
    protein_accession: str
        Uniprot accession ID for a protein
    domain_database: str
        Domain database to search for, default is None
    Returns
    -------
    dictionary: dict
        dictionary with the keys 'name', 'accession', 'num_boundaries', 'boundaries'
        where 'boundaries' is a list of dictionaries with keys 'start' and 'end'

    """
    entry_protein_locations = entry['proteins'][0]['entry_protein_locations']
    if entry_protein_locations is None:
        entry_protein_locations = []

    num_boundaries = len(entry_protein_locations)
    if num_boundaries == 0:
        return None

    dictionary = { 
        'name': entry['metadata']['name'],
        'accession': entry['metadata']['accession'],
        'num_boundaries': num_boundaries,
        'boundaries': [
            {
                'start': bounds['start'],
                'end': bounds['end']
            } 
            for i in range(num_boundaries)
            for bounds in [entry_protein_locations[i]['fragments'][0]]
            if bounds['dc-status'] == "CONTINUOUS"
        ]
    }
    if 'extra_fields' in entry and 'short_name' in entry['extra_fields']:
        dictionary['short'] = entry['extra_fields']['short_name']
    return dictionary

def fetch_InterPro_json(protein_accessions):
    """
    Instantiates an api fetch to InterPro database for domains. Returns a dictionary of those repsonses with keys equal 
    to the protein accession run.
    
    """
    interpro_url = "https://www.ebi.ac.uk/interpro/api"
    extra_fields = ['short_name']
    response_dict = {}
    # code you want to evaluate
    with requests.Session() as session:
        for protein_accession in protein_accessions:
            url = interpro_url + "/entry/interpro/protein/uniprot/" + protein_accession + "?extra_fields=" + ','.join(extra_fields)
            try:
                response_dict[protein_accession] = session.get(url).json()

            except Exception as e:
                if session.get(url).status_code == 204:
                    response_dict[protein_accession] = {}
                    print(f"An empty response was received for {protein_accession} resulting in empty domain architecture.")
                else:
                    print(f"Error processing {protein_accession}: {e}")  # Debugging line
    return response_dict


def get_domains(protein_accessions):
    """
    Given a uniprot accession (protein_accession), return a list of domain dictionaries
    each domain dictionary has keys 'name', 'start', 'end', 'accession' (InterPro ID), 'num_boundaries' (number of this type found)
    These domains are in the order as returned by InterPro, where InterPro returns the parent nodes first. Once 
    we find domains that begin to overlap in the API response, we stop adding those to the final set of domains. 

    Parameters
    ----------
    protein_accession: str
        Uniprot accession ID for a protein
    Returns
    -------
    domain_dict: dict of list of dicts
        outer key values are the individual protein accessions
        these point to a list of dictionaries, each dictionary is a domain entry with keys 'name', 'start', 'end', 'accession', 'num_boundaries'
        this list is ordered by start positions of domains
    domain_string_dict: dict of lists of strings
        outer key values are the individual protein accessions
        these point to a list of stirngs, each string is information for the domain in this manner
        short_name:interpro_id:start:end
    arch_dict: dict of strings
        outer key values are the individual protein accessions
        these point to a string that is the domain architecture, | separated list of domain names
    """
    resp_dict = fetch_InterPro_json(protein_accessions) #pack and unpack as a list for a single domain fetch
    domain_dict = {}
    domain_string_dict = {}
    arch_dict = {}
    for protein_accession in resp_dict:
        domain_dict[protein_accession], domain_string_dict[protein_accession], arch_dict[protein_accession] = get_domains_from_response(resp_dict[protein_accession])
    return domain_dict, domain_string_dict, arch_dict

def get_domains_from_response(resp):
    """
    Given a response from the InterPro API for a single protein search, return a list of domain dictionaries
    each domain dictionary has keys 'name', 'start', 'end', 'accession' (InterPro ID), 'num_boundaries' (number of this type found)
    These domains are in the order as returned by InterPro, where InterPro returns the parent nodes first. Once 
    we find domains that begin to overlap in the API response, we stop adding those to the final set of domains. 
    This returns the ordere list of domains and a list of domain information strings, based on start site.

    Parameters
    ----------
    resp: dict
        response from the InterPro API (json)
    Returns
    -------
    sorted_domain_list: list
        list of dictionaries, each dictionary is a domain entry with keys 'name', 'start', 'end', 'accession', 'num_boundaries'
    domain_string_list: list
        list of domain information short_name:id:start:end
    domain_arch: string
        domain architecture as a string, | separated list of domain names
    """
    # An empty response passed from the Interpro API will bypass all this
    if resp:
        entry_results = resp['results']
        d_dict = {} # Dictionary to store domain information for each entry
        d_resolved = []
        for i, entry in enumerate(entry_results):
        #for i, entry in enumerate(entry_list):
            if entry['metadata']['type'] == 'domain': #get domain level only features
                d_dict[i] = collect_data(entry)
        
        values = list(d_dict.keys())
        if values:
            d_resolved+=return_expanded_domains(d_dict[values[0]]) # a list now: kick off the resolved domains, now start walking through and decide if taking a new domain or not.
        
        for domain_num in values[1:]:
        
            d_resolved = resolve_domain(d_resolved, d_dict[domain_num])

        #having resolved, now let's sort the list and get the domain string information
        sorted_domain_list, domain_string_list, domain_arch = sort_domain_list(d_resolved)
    else:
        sorted_domain_list = []
        domain_string_list = []
        domain_arch = ''
    return sorted_domain_list, domain_string_list, domain_arch
  
def return_expanded_domains(domain_entry):
    """
    Given a domain entry, such as from collect_data, return a list of expanded domains where there is only 
    one boundary per set. This will reset the dictionary, such that instead of 'boundaries' with a list of ['start': x, 'end': y]
    it will be a single boundary with 'start': x, 'end': y as keys in the dictionary.
    """
    boundary_list = domain_entry['boundaries']
    domain_new = domain_entry.copy()
    domain_new['start'] = boundary_list[0]['start']
    domain_new['end'] = boundary_list[0]['end']
    domain_list = []
    domain_list.append(domain_new)
    #make a new dictionary with the start and end values of subsequent domains (Then go through and pop boundaries off all)
    if len(boundary_list) > 1:
        for i in range(1, len(boundary_list)):
            domain_temp = domain_entry.copy()
            domain_temp['start'] = boundary_list[i]['start']
            domain_temp['end'] = boundary_list[i]['end']
            domain_list.append(domain_temp)
    for domain in domain_list: #remove the boundaries term now
        domain.pop('boundaries')
    return domain_list

def resolve_domain(d_resolved, dict_entry, threshold=0.5):
    """
    Given a list of resolved domains and a new domain entry, resolve the new domain entry with the existing domains
    Keep the new domain entry if it does not overlap by more than threshold% with any existing domain. Default threshold is
    50% (or 0.5)
    Parameters
    ----------
    d_resolved: list
        list of dictionaries, each dictionary is a domain entry with keys 'name', 'start', 'end', 'accession', 'num_boundaries'
    dict_entry: dict
        dictionary of domain information with that comes from collect_data on an entry. This function expands multiple domain entries, meaning that these are prioritized as they are encountered first
    threshold: float
        threshold for rejecting domains by overlap, default is 0.5, should be between 0 and 1
    Returns
    -------
    d_resolved: list
        list of dictionaries, each dictionary is a domain entry with keys 'name', 'start', 'end', 'accession', 'num_boundaries'
    """
    # d_resolved is a list of dictionaries, each dictionary is a domain entry
    #setup the existing boundaries that are in d_resolved
    if threshold < 0 or threshold > 1:
        threshold = 0.5 #set to default
        print("WARN: Threshold must be between 0 and 1 for rejecting domains by overlap, setting to default of 0.5")

    boundary_array = []
    for domain in d_resolved:
        boundary_array.append(set(range(domain['start'], domain['end'])))
    
    #expand the dict_entry as well to a list
    new_domains = return_expanded_domains(dict_entry)
    #print("DEBUG: have these new domains")
    #print(new_domains)

    #first expand if multiple boundaries exist in the dict_entry
    for domain in new_domains:
        range_new = set(range(domain['start'], domain['end']))
        found_intersecting = False
        for range_existing in boundary_array:
            #check if the set overlap between the new range and the existing range is greater than 
            # 50% of the new range. If so, do not add the new range.
            if len(range_new.intersection(range_existing))/len(min(range_new, range_existing)) > threshold:
                found_intersecting = True
                break
        if not found_intersecting:
            d_resolved.append(domain)
    return d_resolved

def sort_domain_list(domain_list):
    """
    Given a list of resolved domains, return a string of domain information short_name:id:start:end and a sorted list of
    the domains according to the start site. 

    Parameters
    ----------
    domain_list: list
        list of domain dictionaries, where the dictonaries have keys 'name', 'accession', 'short', 'start', 'end'
    Returns
    -------
    sorted_domain_list: list
        list of domain dictionaries, now sorted by the start positions.
    domain_string_list: list
        list of domain information short_name:id:start:end
    domain_arch: string
        domain architecture as a string, | separated list of domain names

    """
    domdict = {}
    for domain in domain_list:
        start = int(domain['start'])
        if start in domdict:
            print("ERROR: More than one domain have the same start position!")
        domdict[start] = domain

    sorted_dict = dict(sorted(domdict.items(),reverse=False))
    sorted_domain_list = []
    domain_string_list = []
    domain_arch_names = []
    for key, value in sorted_dict.items():
        sorted_domain_list.append(value)
        domain_string_list.append(value['short']+':'+value['accession']+':'+str(key)+':'+str(value['end']))
        domain_arch_names.append(value['short'])
    return sorted_domain_list, domain_string_list, '|'.join(domain_arch_names)


def generateDomainMetadata_wfilter(uniprot_accessions):
    """
    NO LONGER USED: Leaving this here, this is code that focuses on hierarchy and children to determine the list 
    of proteins. May 22, 2024 KMN. Written by LC. 

    Given a list of uniprot accessions, return a dictionary of protein accessions containing domain metadata.
    This function will return a dictionary with keys as the protein accession and values as a list of dictionaries
    where each dictionary contains domain metadata. This metadata is collected from the InterPro API and includes
    the keys 'interpro', 'num_children'. The 'interpro' key contains a dictionary with keys 'name', 'accession', 'num_boundaries', 'boundaries'
    where 'boundaries' is a list of dictionaries with keys 'start' and 'end'. If 'short_name' is in the extra fields, this will be added to the dictionary as 'short'
    This version of code uses hierarchy and children nodes to select InterPro domains.

    Parameters
    ----------
    uniprot_accessions: list
        list of uniprot accessions
    Returns
    -------
    metadata: dict
        dictionary of protein accessions containing domain metadata
    """


    interpro_url = "https://www.ebi.ac.uk/interpro/api"
    extra_fields = ['hierarchy', 'short_name']
    metadata = {}

    with requests.Session() as session:
        for protein_accession in uniprot_accessions:
            #print(f"Processing {protein_accession}")  # Debugging line

            url = interpro_url + "/entry/interpro/protein/uniprot/" + protein_accession + "?extra_fields=" + ','.join(extra_fields)

            try:
                resp = session.get(url).json()
            except Exception as e:
                print(f"Error processing {protein_accession}: {e}")  # Debugging line
                metadata[protein_accession] = []  
                continue             

            current_accession = [protein_accession]

            # pulling out top hierarchy domains from extra fields
            entry_results = resp['results']
            top_hierarchy = [] # Getting top hierarchy domains for each entry in the results
            num_children_list = [] # List to store number of children for each domain
            for i, entry in enumerate(entry_results):
                if 'extra_fields' in entry and 'hierarchy' in entry['extra_fields']:
                    if 'children' in entry['extra_fields']['hierarchy']:
                        num_children = len(entry['extra_fields']['hierarchy']['children'])
                        num_children_list.append(num_children)
                    top_hierarchy.append(entry['extra_fields']['hierarchy']['accession'])
                else: 
                    print(f"Error processing {protein_accession}: No hierarchy found")
            #print(top_hierarchy)
            entry_list = [
                {'interpro': data, 'num_children': num_children_list[i]}
                for i, entry in enumerate(resp['results'])
                if entry['metadata']['type'] == 'domain' and entry['metadata']['accession'] in top_hierarchy
                for data in [collect_data(entry)]
                if data is not None
            ]
        
            metadata[protein_accession] = entry_list 
           # print(entry_list)

    return metadata

# Filtering returned metadata 

def filter_domains(metadata, threshold=0.15):
    """NO LONGER USED: Leaving this here, this is code that focuses on hierarchy and children to determine the list 
    of proteins. May 22, 2024 KMN. Written by LC. """
    for protein_accession, domains in metadata.items():
        # Sort by end position and then by size in descending order
        domains.sort(key=lambda x: (-x['interpro']['boundaries'][0]['end'], -(x['interpro']['boundaries'][0]['end'] - x['interpro']['boundaries'][0]['start']), -x['num_children']))
        filtered_domains = []
        for domain in domains:
            overlap = False
            for existing_domain in filtered_domains:
                start_max = max(domain['interpro']['boundaries'][0]['start'], existing_domain['interpro']['boundaries'][0]['start'])
                end_min = min(domain['interpro']['boundaries'][0]['end'], existing_domain['interpro']['boundaries'][0]['end'])
                overlap_length = max(0, end_min - start_max)
                domain_length = domain['interpro']['boundaries'][0]['end'] - domain['interpro']['boundaries'][0]['start']
                if overlap_length / domain_length > threshold:
                    # If the current domain has more children, replace the existing domain
                    if domain['num_children'] > existing_domain['num_children']:
                        filtered_domains.remove(existing_domain)
                        filtered_domains.append(domain)
                    overlap = True
                    break
            if not overlap:
                filtered_domains.append(domain)
        filtered_domains.sort(key=lambda x: x['interpro']['boundaries'][0]['start'])
        metadata[protein_accession] = filtered_domains
    return metadata


# NO LONGER USED, this string generation goes with the code that focuses on hierarch. 
def generate_domain_metadata_string_list(metadata, uniprot_list):
    """
    NO LONGER USED: Leaving this here, this is code that focuses on hierarchy and children to determine the list 
    of proteins. May 22, 2024 KMN. Written by LC.

    Condenses protein metadata into strings and outputs them as a list corresponding with indices of accessions in uniprot_list
    Parameters
    ----------
    processed_metadata: dict
        dictionary of protein accessions containing domain metadata. Following hierarchy and overlap filtering
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
        domain_metadata_list = metadata[accession]
        # Sort domain metadata by start position
        domain_metadata_list.sort(key=lambda x: x['interpro']['boundaries'][0]['start'])
        for domain_dict in domain_metadata_list:
            short_name = domain_dict['interpro']['short']
            id = domain_dict['interpro']['accession']
            # Iterate over boundaries
            # here, also sort the domains in the list according to start positions.
            for boundary in domain_dict['interpro']['boundaries']:
                start = boundary['start']
                end = boundary['end']
                metadata_string = short_name+':'+id+':'+str(start)+':'+str(end)
                string_list.append(metadata_string)
        sorted_domain_list = sort_domain_list(string_list)
        domain_arch = return_domain_architecture(string_list)
        domain_arch_list.append(domain_arch)
        metadata_string_list.append(';'.join(sorted_domain_list))
    return metadata_string_list, domain_arch_list





# def sort_domain_list(domain_string_list):
#     """
#     Given a domain_list, list of domain information short_name:id:start:end, 
#     sort the list according to start positions of the domains and return a new list

#     Parameters
#     ----------
#     domain_string_list: list
#         list of domain information short_name:id:start:end as a string, array of strings
#     Returns
#     -------
#     sorted_domain_string_list: list
#         list of domain information sorted according to start positions as short_name:id:start:end as a string, array of strings

#     """
#     domdict = {}
#     for domain_info in domain_string_list:
#         name, uniprot_id, start, end = domain_info.split(':')
#         start = int(start)
#         domdict[start] = name, uniprot_id, end

#     sorted_dict = dict(sorted(domdict.items(),reverse=False))
#     sorted_domain_list = []
#     for key, values in sorted_dict.items():
#         domain = values[0]+':'+values[1]+':'+str(key)+':'+str(values[2])
#         sorted_domain_list.append(domain)
#     return sorted_domain_list

def appendRefFile(input_RefFile, outputfile):
    '''
    Takes a reference file generated made by CoDIAC.UniProt.makeRefFile and adds 
    interpro domain metadata as a new column (i.e. this appends domain information defined by InterPro
    to the Uniprot reference)

    Parameters
    ----------
    input_RefFile: string
        name of the input reference file generated from the  makeRefFile function in CODAC.py
    outputfile: string
        name of the file to be outputted by this function
        
    Returns
    -------
    df: Pandas Dataframe
        In addition to printing the dataframe to a CSV file (as defined by outputfile)
        this returns the dataframe that is prented
    '''
    df = pd.read_csv(input_RefFile)
    uniprotList = df['UniProt ID'].to_list()
    print("Fetching domains..")
    domain_dict, domain_string_dict, domain_arch_dict = get_domains(uniprotList)
    print("Appending domains to file..")    
    for i in range(len(uniprotList)):
        df.at[i, 'Interpro Domains'] = ';'.join(domain_string_dict[uniprotList[i]])
        df.at[i, 'Interpro Domain Architecture'] = domain_arch_dict[uniprotList[i]]
    #df['Interpro Domains'] = metadata_string_list
    #df['Interpro Domain Architecture'] = domain_arch_list
    df.to_csv(outputfile, index=False)
    print('Interpro metadata succesfully incorporated')
    return df

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