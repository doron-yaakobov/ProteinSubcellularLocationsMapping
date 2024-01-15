import requests
import re


def get_protein_data(protein_id: str) -> str:
    url = f"https://www.uniprot.org/uniprot/{protein_id}.xml"
    response = requests.get(url)
    response.raise_for_status()
    return response.content.decode("utf-8")


def extract_subcellular_locations(protein_data_from_xml: str) -> list:
    subcellular_location_tag_regex = r'(?<=<comment type="subcellular location">)(.*?)(?=</comment>)'
    location_tag_regex = r'(?<=<location>)(.*?)(?=</location>)'
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
