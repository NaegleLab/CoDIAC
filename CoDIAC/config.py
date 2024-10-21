import json
import requests
from pathlib import Path

# Or use this test set of ids that have small files (To use, delete the '#' in the next line)
item_ids = [27264570]

#Set the base URL
BASE_URL = 'https://api.figshare.com/v2'

def install_resource_files(package_directory):

    """ Downloads large data files from Figshare into respective directories that contains data generated from Phosphosite plus and proteomeScout databases. """ 
    
    file_info = [] #a blank list to hold all the file metadata
    
    for i in item_ids:
        r = requests.get(BASE_URL + '/articles/' + str(i) + '/files')
        file_metadata = json.loads(r.text)
        for j in file_metadata: #add the item id to each file record- this is used later to name a folder to save the file to
            j['item_id'] = i
            file_info.append(j) #Add the file metadata to the list
    
    #Download each file to a subfolder named for the article id and save with the file name
    for k in file_info:
        response = requests.get(BASE_URL + '/file/download/' + str(k['id']))
        if k['name'] == 'data.tsv':
            # package_directory = os.path.dirname(os.path.abspath(__file__))
            DIRECTORY = package_directory + '/data/proteomescout_everything_20190701'
            Path(DIRECTORY).mkdir(exist_ok=True) 
            open(DIRECTORY + '/' + k['name'], 'wb').write(response.content)
        if k['name'] == 'phosphositeplus_data.csv':
            # package_directory = os.path.dirname(os.path.abspath(__file__))
            DIRECTORY = package_directory + '/data' 
            Path(DIRECTORY).mkdir(exist_ok=True) 
            open(DIRECTORY + '/' + k['name'], 'wb').write(response.content)
        