import requests
import re
import pandas as pd

# region DF config
FILE = "83344-9 proteinGroups-perseus.xlsx"
SHEET = "83344-9 proteinGroups-perseus"
COLUMNS = "A:D"
HEADERS = 1
NUM_OF_ROWS = 619


# endregion
def get_protein_data(protein_id: str) -> str:
    url = f"https://www.uniprot.org/uniprot/{protein_id}.xml"
    response = requests.get(url)
    response.raise_for_status()
    return response.content.decode("utf-8")


def extract_subcellular_locations(protein_data_from_xml: str) -> list:
    subcellular_location_tag_regex = r'(?<=<comment type="subcellular location">)(.*?)(?=</comment>)'
    location_tag_regex = r'<location.*?>(.*?)<\/location>'
    locations = list()

    matches = re.findall(subcellular_location_tag_regex, protein_data_from_xml, re.DOTALL)
    for match in matches:
        subcellular_locations = re.findall(location_tag_regex, match, re.DOTALL)
        for subcellular_location in subcellular_locations:
            # print(subcellular_location)
            locations.append(subcellular_location)
    return locations


def get_subcellular_locations(protein_id: str) -> list:
    protein_data = get_protein_data(protein_id)
    return extract_subcellular_locations(protein_data)


# locations = get_subcellular_locations('P19338')
# print(locations)

df = pd.read_excel(FILE, sheet_name=SHEET, usecols=COLUMNS, header=HEADERS, nrows=NUM_OF_ROWS)

df['Protein IDs'] = df['Protein IDs'].apply(lambda x: x.split(';') if ';' in x else [x])
for idx, row in df.iterrows():
    for protein_id in row['Protein IDs']:
        subcellular_locations = get_subcellular_locations(protein_id)
        print(subcellular_locations)
        pass
