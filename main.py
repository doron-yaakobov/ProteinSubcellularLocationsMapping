import requests
import re
import pandas as pd
import json
from collections import Counter
import matplotlib.pyplot as plt

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


def build_gene_locations_json(df):
    gene_locations = {}
    for idx, row in df.iterrows():
        gene = row['Gene names']
        for protein_id in row['Protein IDs']:
            subcellular_locations = get_subcellular_locations(protein_id)
            if gene not in gene_locations:
                gene_locations[gene] = subcellular_locations
            else:
                gene_locations[gene] = list(set(gene_locations[gene] + subcellular_locations))

        if not gene_locations[gene]:  # empty
            print(f"ERROR! couldn't define subcellular locations for {gene}:\t{row['Protein IDs']}")

        with open('gene_locations.json', 'w') as f:
            json.dump(gene_locations, f)

    print(f"Json Dumped")


def calculate_locations_counts(gene_locations: dict):
    all_values = [value for values in gene_locations.values() for value in values]
    return Counter(all_values)


def calculate_locations_counts_uniques(gene_locations: dict):
    result = Counter()

    for key, values in gene_locations.items():
        if len(values) == 1:
            result[values[0]] += 1

    return result

if __name__ == "__main__":
    df = pd.read_excel(FILE, sheet_name=SHEET, usecols=COLUMNS, header=HEADERS, nrows=NUM_OF_ROWS)
    # df = df[431:] # splitting the work
    df['Protein IDs'] = df['Protein IDs'].apply(lambda x: x.split(';') if ';' in x else [x])

    # build_gene_locations_json(df)

    # region merge jsons:
    ''' use only if you splited the jsons, else - disable one of the loads.'''

    with open(r'gene_locations.json', 'r') as f:
        data = json.load(f)
    with open(r'gene_locations_continued.json', 'r') as f:
        data_continued = json.load(f)

    subcellular_locations = {**data, **data_continued}
    # endregion

    # region find empty locations:
    keys_with_empty_values = [key for key, value in subcellular_locations.items() if not value]
    with open('genes_without_locations.txt', 'w') as file:
        for item in keys_with_empty_values:
            file.write(f"{item}\n")

    print(f"found {len(keys_with_empty_values)} genes with empty locations. saved.")
    # endregion

    subcellular_locations_count_with_duplications = calculate_locations_counts(subcellular_locations)

    # region plot with duplications:
    with_duplicates_df = \
        pd.DataFrame({'subcellular_location': list(subcellular_locations_count_with_duplications.keys()),
                      'amount': list(subcellular_locations_count_with_duplications.values())}) \
            .sort_values(by='amount', ascending=False)
    with_duplicates_df.to_excel("genes_locations_with_duplicates.xlsx", index=False)
    print("Saved genes_locations_with_duplicates.xlsx \nplotting head.")

    with_duplicates_df_head = with_duplicates_df.head(5)
    with_duplicates_df_head.plot.pie(y='amount', labels=with_duplicates_df_head['subcellular_location'],
                                     autopct='%1.1f%%', legend=False)
    plt.ylabel("")
    plt.xlabel("")
    plt.title("Genes Source Locations With Duplicates")
    plt.savefig("Genes Source Locations With Duplicates.png")
    print("Figure 'Genes Source Locations With Duplicates.png' saved.")
    # endregion


    subcellular_locations_count_unique = calculate_locations_counts_uniques(subcellular_locations)

    # region plot Uniques:
    uniques_df = \
        pd.DataFrame({'subcellular_location': list(subcellular_locations_count_unique.keys()),
                      'amount': list(subcellular_locations_count_unique.values())}) \
            .sort_values(by='amount', ascending=False)
    uniques_df.to_excel("genes_locations_uniques.xlsx", index=False)
    print("Saved genes_locations_uniques.xlsx \nplotting head.")

    uniques_df_head = uniques_df.head(5)
    uniques_df_head.plot.pie(y='amount', labels=uniques_df_head['subcellular_location'],
                                     autopct='%1.1f%%', legend=False)
    plt.ylabel("")
    plt.xlabel("")
    plt.title("Genes Source Locations Without Duplicates")
    plt.savefig("Genes Source Locations Without Duplicates.png")
    print("Figure 'Genes Source Locations Without Duplicates.png' saved.")
    # endregion

    # plt.pie(sizes, labels=labels)
    # plt.title('Pie Chart')
    plt.show()

    # endregion

    pass
