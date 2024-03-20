from collections import defaultdict
import requests
import logging
import os
import pandas as pd
import copy

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

        if species not in all_species and Scientific_name != species:
            continue

        UNIPROT_ID_LIST.append(Uniprot_Accession)
        uniprot_dict[Uniprot_Accession] = Scientific_name
        species_dict[Scientific_name][source_database].append(Uniprot_Accession)

    print(f"Fetched {len(UNIPROT_ID_LIST)} Uniprot IDs linked to {interpro_ID}, where count expected to be {count_data}")
    return(UNIPROT_ID_LIST, species_dict)


def collect_data(entry, protein_accession, domain_database=None):
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


def generateDomainMetadata_wfilter(uniprot_accessions):
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
                for data in [collect_data(entry, current_accession)]
                if data is not None
            ]
        
            metadata[protein_accession] = entry_list 
           # print(entry_list)

    return metadata

# Filtering returned metadata 

def filter_domains(metadata, threshold=0.15):
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



def generate_domain_metadata_string_list(metadata, uniprot_list):
    """
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
            for boundary in domain_dict['interpro']['boundaries']:
                start = boundary['start']
                end = boundary['end']
                metadata_string = short_name+':'+id+':'+str(start)+':'+str(end)
                string_list.append(metadata_string)
        domain_arch = return_domain_architecture(string_list)
        domain_arch_list.append(domain_arch)
        metadata_string_list.append(';'.join(string_list))
    return metadata_string_list, domain_arch_list


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
    domainMD = generateDomainMetadata_wfilter(uniprotList)
    print("Processing domains...")
    #processed_MD = process_proteins(domainMD)
    filtered_MD = filter_domains(domainMD)
    metadata_string_list, domain_arch_list = generate_domain_metadata_string_list(filtered_MD, uniprotList)
    
    df['Interpro Domains'] = metadata_string_list
    df['Interpro Domain Architecture'] = domain_arch_list
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