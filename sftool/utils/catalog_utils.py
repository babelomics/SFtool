
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 19:07:52 2023

@author: Javier Perez Florido, Edurne Urrutia
"""
import csv
import json
from natsort import natsorted
import requests
import vcfpy

def read_csv(in_csv, category):
    """
    Read a CSV file and store information in a dictionary

    Args:
        in_csv (str): Path to CSV file

    Returns:
        dict, list: A dictionary with information from CSV file and a list of gene symbol
    """
    # Create dictionary and gene list to store info
    if category == 'PR':
        cat_str = 'personal'
    elif category == 'RR':
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
    if assembly == "GRCh37":
        server = "https://grch37.rest.ensembl.org"  # GRCh37 (hg19) Ensembl server
    elif assembly == "GRCh38":
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


def write_bed_file(assembly, genes_lst, bed_file, has_chr_prefix):
    """
    Write gene information to a BED format file

    Args:
        assembly (str): reference genome version ("GRCh37" or "GRCh38").
        genes_lst (list): Gene symbol list
        bed_file (str): Path to BED file to be generated
        has_chr_prefix (boolean): Whether VCF file has chr prefix for CHROM value

    Returns:
        None
    """

    chr_prefix = ''
    if has_chr_prefix:
        chr_prefix = 'chr'


    gene_coords = []
    # Some genes change name between assemblies

    for gene in genes_lst:
        gene_query = gene
        if gene == 'MMUT' and assembly == 'GRCh37':
            gene_query = 'MUT'
        elif gene == 'ELP1' and assembly == 'GRCh37':
            gene_query = 'IKBKAP'
        elif gene == 'G6PC1' and assembly == 'GRCh37':
            gene_query = 'G6PC'
        elif gene == 'GBA1' and assembly == 'GRCh37':
            gene_query = 'GBA'

        gene_pos = get_gene_location_ensembl(gene_query, assembly)

        gene_coords.append((gene_pos['Chromosome'], int(gene_pos['Start']), int(gene_pos['End']), gene))
    sorted_coords = natsorted(gene_coords)

    with open(bed_file, "w") as fdw:
        for chrom, start, end, gene in sorted_coords:
            fdw.write(f"{chr_prefix}{chrom}\t{start}\t{end}\t{gene}\n")
        print(f"BED file '{bed_file}' generated successfully.")

def has_chr_prefix_vcf(vcf_file):
    reader = vcfpy.Reader.from_path(vcf_file)
    for record in reader:
        return record.CHROM.startswith("chr")  # Check first variant
    return False  # If no variants are found


def build_catalog_resources(category, assembly, category_geneset_file, bed_file, json_file, vcf_file):
    """
    Main function: from a CSV file, creates a JSON and a BED files

    Args:
        category (str): category, either PR or RR
        assembly (str): assembly version ("GRCh37" or "GRCh38").
        category_geneset_file (str): Path to CSV file for the given category
        bed_file (str): Path to BED file to be generated
        json_file (str): Path to JSON file to be generated
        vcf_file (str): Path to original VCF file

    Returns:
        None
    """

    print("Creating BED and JSON files for " + category + " catalogue.")

    try:

        # Check whether chr prefix is used in VCF file
        has_chr_prefix = has_chr_prefix_vcf(vcf_file)

        # Read CSV and store it in the dictionary
        genes_dct, genes_lst = read_csv(category_geneset_file, category)

        # Write a BED file
        write_bed_file(assembly, genes_lst, bed_file, has_chr_prefix)

        # Write a JSON file
        with open(json_file, 'w') as fdw:
            json.dump(genes_dct, fdw, indent = 4)
            print(f"JSON file '{json_file}' generated successfully.")


    except Exception as e:
        print(f"An error occurred: {str(e)}")