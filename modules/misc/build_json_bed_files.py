
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 19:07:52 2023

@author: Javier Perez Florido, Edurne Urrutia
"""
import csv
import json
from natsort import natsorted
import requests

def read_csv(in_csv, category):
    """
    Read a CSV file and store information in a dictionary

    Args:
        in_csv (str): Path to CSV file

    Returns:
        dict, list: A dictionary with information from CSV file and a list of gene symbol
    """
    # Create dictionary and gene list to store info
    if category == 'pr':
        cat_str = 'personal'
    elif category == 'rr':
        cat_str = 'reproductive'
    genes_dct = {
        "category": f"Secondary findings of {cat_str} risk",
        "genes": []
    }

    genes_lst = []

    # Read CSV file and store it in the dictionary
    with open(in_csv, 'r', encoding='latin1') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            gene_symbol = row['Gene']
            gene_info = {
                'gene_symbol': row['Gene'],
                'phenotype': row['Phenotype'],
                'ACMG_version': row.get('ACMG SF List Version', '') ,
                'OMIM_disorder': row['OMIM Disorder'],
                'inheritance': row['Inheritance'],
                'variants_to_report': row.get('Variants to Report', ''),
                'specific_variant_GRCh38': row.get('Specific variant GRCh38', ''),
                'specific_variant_GRCh37': row.get('Specific variant GRCh37', ''),
                'specific_consequence': row.get('Specific consequence', ''),
            }
            genes_dct["genes"].append(gene_info)
            genes_lst.append(gene_symbol)

    return(genes_dct, genes_lst)


def get_gene_location_ensembl(gene_symbol, assembly):
    """Fetch chromosome location of a gene from Ensembl REST API."""

    # Define the Ensembl server based on genome version
    if assembly == "37":
        server = "https://grch37.rest.ensembl.org"  # GRCh37 (hg19) Ensembl server
    elif assembly == "38":
        server = "https://rest.ensembl.org"  # GRCh38 (hg38) Ensembl server

    url = f"{server}/lookup/symbol/human/{gene_symbol}?content-type=application/json"

    response = requests.get(url)

    if response.status_code != 200:
        return f"Error fetching data for {gene_symbol} ({assembly}): {response.status_code}"

    data = response.json()

    result = {}

    result['Gene_symbol'] = gene_symbol
    result['Chromosome'] = data['seq_region_name']
    result['Start'] = data['start']
    result['End'] = data['end']

    return result


def write_bed_file(assembly, genes_lst, category, categories_path):
    """
    Write gene information to a BED format file

    Args:
        assembly (str): reference genome version ("37" or "38").
        genes_lst (list): Gene symbol list
        category (str): category
        categories_path (str): Path to category directory

    Returns:
        None
    """
    gene_coords = []
    #### Some genes change name between assemblies
    for gene in genes_lst:
        if gene == 'MMUT' and assembly == '37':
            gene = 'MUT'
        elif gene == 'ELP1' and assembly == '37':
            gene = 'IKBKAP'
        elif gene == 'G6PC' and assembly == '38':
            gene = 'G6PC1'
        elif gene == 'GBA' and assembly == '38':
            gene = 'GBA1'

        gene_pos = get_gene_location_ensembl(gene, assembly)

        gene_coords.append((gene_pos['Chromosome'], int(gene_pos['Start']), int(gene_pos['End']), gene))
    sorted_coords = natsorted(gene_coords)

    filename = f"{categories_path}{category.upper()}/{category}_risk_genes_GRCh{assembly}.bed"
    with open(filename, "w") as bed_file:
        for chrom, start, end, gene in sorted_coords:
            bed_file.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
        print(f"BED file '{filename}' generated successfully.")


def generate_json_from_fg_csv(csv_file, assembly, categories_path):
    """
    Generate a JSON file from a CSV file containing pharmacogenetic data

    Args:
        csv_file (str): Path to CSV input file
        assembly (str): Reference genome version ("37" or "38").
        categories_path: Path to categories directory

    Returns:
        None
    """
    try:

        # Create dictionary and gene list to store info
        variants_dct = {
            "category": "Pharmacogenetic secondary findings",
            "variants": {}
        }

        # Read CSV file and store it in the dictionary
        with open(csv_file, 'r') as file:
            csv_reader = csv.DictReader(file)
            for row in csv_reader:
                alts = row['Alternative'].split(',')
                if len(alts) > 1:
                    for i in alts:
                        variant = row['Chromosome'] + ':' + row['Position'] + ':' + row['Reference'] + ':' + i
                        variant_info = {
                            'gene_symbol': row['Gene'],
                            'rs': row['dbSNP'],
                        }
                        variants_dct["variants"][variant] = variant_info
                else:
                    variant = row['Chromosome'] + ':' + row['Position'] + ':' + row['Reference'] + ':' + row['Alternative']
                    variant_info = {
                        'gene_symbol': row['Gene'],
                        'rs': row['dbSNP'],
                    }
                    variants_dct["variants"][variant] = variant_info

            # Write JSON file
            out_json = f'{categories_path}FG/fg_risk_genes_GRCh{assembly}.json'
            with open(out_json, 'w') as json_file:
                json.dump(variants_dct, json_file, indent = 4)

        print(f"JSON file '{out_json}' generated successfully.")
    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def generate_bed_from_fg_csv(csv_file, assembly, categories_path):
    """
    Create a BED file from a CSV file with pharmacogenetic data

    Args:
        csv_file (str): Path to CSV input file
        assembly (str): Reference genome version ("37" or "38").
        categories_path: Path to category directory

    Returns:
        None
    """
    bed_data = []

    try:
        with open(csv_file, 'r') as file:
            csv_reader = csv.DictReader(file)

            for row in csv_reader:
                chr_name = row['Chromosome']
                position = int(row['Position']) - 1  # BED uses 0-based coordinates
                ref_allele = row['Reference']
                alt_alleles = row['Alternative'].split(',')  # Split multiple alleles

                for alt_allele in alt_alleles:
                    length_diff = len(alt_allele) - len(ref_allele)
                    start_pos = position
                    end_pos = position + max(len(ref_allele), len(alt_allele))  # Take longer allele length
                    bed_entry = (chr_name, start_pos, end_pos, row['Gene'], ref_allele, alt_allele)
                    bed_data.append(bed_entry)

        # Sort bed_data by chromosome, start position, end position, and allele
        sorted_bed_data = natsorted(bed_data)

        bed_filename = f"{categories_path}FG/fg_risk_genes_GRCh{assembly}.bed"

        with open(bed_filename, "w") as bed_file:
            for entry in sorted_bed_data:
                bed_file.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\t{entry[5]}\n")

        print(f"BED file '{bed_filename}' generated successfully.")
    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")


def build_json_bed_files(category, assembly, categories_path, category_geneset_file):
    """
    Main function: from a CSV file, creates a JSON and a BED files

    Args:
        category (str): category, either PR or RR
        assembly (str): assembly version ("37" or "38").
        categories_path (str): Path to category directory
        category_geneset_file (str): Path to CSV file for the given category

    Returns:
        None
    """

    print("Creating BED and JSON files for " + category.upper() + " catalogue.")

    try:
        if category == 'pr' or category == 'rr':
            # Read CSV and store it in the dictionary
            genes_dct, genes_lst = read_csv(category_geneset_file, category)

            # Write a BED file
            write_bed_file(assembly, genes_lst, category, categories_path)

            # Write a JSON file
            out_json = f"{categories_path}{category.upper()}/{category}_risk_genes.json"
            with open(out_json, 'w') as json_file:
                json.dump(genes_dct, json_file, indent = 4)
                print(f"JSON file '{out_json}' generated successfully.")
        else:  # fg risk category
            generate_json_from_fg_csv(category_geneset_file, assembly, categories_path)
            generate_bed_from_fg_csv(category_geneset_file, assembly, categories_path)

    except Exception as e:
        print(f"An error occurred: {str(e)}")