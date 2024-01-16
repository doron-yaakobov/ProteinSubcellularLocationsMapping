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


def get_subcellular_locations(protein_id: str) -> list:
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

    protein_data = get_protein_data(protein_id)
    return extract_subcellular_locations(protein_data)


def build_gene_locations_json(df) -> None:
    gene_locations = {}
    for idx, row in df.iterrows():
        gene = row['Gene names']
        for protein_id in row['Protein IDs']:
            subcellular_locations = get_subcellular_locations(protein_id)
            if gene not in gene_locations:
                gene_locations[gene] = subcellular_locations
            else:
                gene_locations[gene] = list(set(gene_locations[gene] + subcellular_locations))
        # Alert on missing gene locations
        if not gene_locations[gene]:
            print(f"ERROR! couldn't define subcellular locations for {gene}:\t{row['Protein IDs']}")
        # save results
        with open('gene_locations.json', 'w') as f:
            json.dump(gene_locations, f)
    print(f"build_gene_locations_json() - Json Dumped!")


# def calculate_locations_counts(gene_locations: dict):
#     all_values = [value for values in gene_locations.values() for value in values]
#     return Counter(all_values)
#
#
# def calculate_locations_counts_uniques(gene_locations: dict):
#     result = Counter()
#
#     for key, values in gene_locations.items():
#         if len(values) == 1:
#             result[values[0]] += 1
#
#     return result


def count_locations(gene_locations: dict, unique_vals_only: bool = False):
    result = Counter()
    if unique_vals_only:
        for key, values in gene_locations.items():
            if len(values) == 1:
                result[values[0]] += 1
    else:
        all_values = [value for values in gene_locations.values() for value in values]
        result = Counter(all_values)
    return result


def analyze_subcellular_locations(counted_locations: dict, output_xlsx: str, output_png: str) -> None:
    df = pd.DataFrame({
        'subcellular_location': list(counted_locations.keys()),
        'amount': list(counted_locations.values())
    }).sort_values(by='amount', ascending=False)
    df.to_excel(output_xlsx, index=False)
    print(f"Saved {output_xlsx}\n")

    # generate plot:
    df = df.head(5)
    df.plot.pie(y='amount', labels=df['subcellular_location'], autopct='%1.1f%%', legend=False)
    plt.ylabel("")
    plt.xlabel("")
    plt.title(f"{output_png}")
    plt.savefig(f"{output_png}.png")
    print(f"Figure '{output_png}.png' saved.")


if __name__ == "__main__":
    # region get df
    df = pd.read_excel(FILE, sheet_name=SHEET, usecols=COLUMNS, header=HEADERS, nrows=NUM_OF_ROWS)
    df['Protein IDs'] = df['Protein IDs'].apply(lambda x: x.split(';') if ';' in x else [x])
    # endregion

    # build_gene_locations_json(df)

    # region read gene locations.json
    with open(r'gene_locations.json', 'r') as f:
        data = json.load(f)
    subcellular_locations = {**data}

    # region merge jsons:
    ''' use only if you splited the jsons, else - disable.'''
    with open(r'gene_locations_continued.json', 'r') as f:
        data_continued = json.load(f)

    subcellular_locations = {**data, **data_continued}
    # endregion

    df = pd.DataFrame.from_dict(subcellular_locations, orient='index')
    df = df.transpose()
    df.to_excel('Sub-Cellular Locations.xlsx', index=False)
    # endregion

    # region Identify empty locations:
    genes_without_locations = [key for key, value in subcellular_locations.items() if not value]
    with open('genes_without_locations.txt', 'w') as file:
        for item in genes_without_locations:
            file.write(f"{item}\n")

    print(f"found {len(genes_without_locations)} genes with empty locations. saved.")
    # endregion

    subcellular_locations_count_with_duplications = count_locations(subcellular_locations)

    analyze_subcellular_locations(subcellular_locations_count_with_duplications,
                                  "genes_locations_with_duplicates.xlsx",
                                  "Genes Source Locations With Duplicates")

    subcellular_locations_count_unique = count_locations(subcellular_locations, unique_vals_only=True)
    analyze_subcellular_locations(subcellular_locations_count_unique,
                                  "genes_locations_uniques.xlsx",
                                  "Genes Source Locations Without Duplicates")

    plt.show()

    # endregion
