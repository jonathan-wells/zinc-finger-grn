#!/usr/bin/env python3

import requests

species_dict = {
    'Limnephilus_lunatus': '1218281', # Cinnamon sedge
    'Ctenocephalides_felis': '7515', # Cat flea
    'Rana_temporaria': '8407', # Common frog
    'Scyliorhinus_canicula': '7830', # Small-spotted catshark
    'Archivesica_marissinica': '2291877' # Cold seep clam 
    }

url = 'https://dfam.org/api/families'
for species, taxid in species_dict.items():
    params = {
        'clade': taxid,
        'format': 'fasta',
        'clade_relatives': 'both',
        'include_raw': 'true'
    }
    response = requests.get(url, params=params)
    results = response.json()['body']
    with open(f'{species}_dfam_library.fa', 'w') as outfile:
        for record in results:
            outfile.write(record)
