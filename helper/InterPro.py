import requests
import logging

'''Given an Interpro Accession ID, this class returns a list of Uniprot Accession IDs associated to this Interpro Acc ID'''
    
def fetch_uniprotid(InterPro_ID, specie='Homo sapiens'):
    
    fetch_all_results = []
    
    base_url = f"https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/InterPro/{InterPro_ID}/?page_size=200"

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

        if Scientific_name == specie:

            UNIPROT_ID_LIST.append(Uniprot_ACC)
    print('Found {count} Uniprot IDs linked to {IPR_ID}'.format(count=len(UNIPROT_ID_LIST), IPR_ID=InterPro_ID))        
    return(UNIPROT_ID_LIST)