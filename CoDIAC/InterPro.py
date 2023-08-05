import requests
import logging
import os
import requests
import logging


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