# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 22:45:05 2023

@author: kindi
"""
import csv
import json
from natsort import natsorted


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
        
        print(f"BED file '{bed_filename}' generated and sorted.")
    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def get_json_bed_fg(assembly, categories_path):
    """
    Main function that reads a CSV file with pharmacogenetic info and create a JSON and a BED file

    Args:
        assembly (str): reference genome version ("37" or "38").

    Returns:
        None
    """
    print("Creating BED and JSON files for FG catalogue.")

    try:
        csv_file = f"{categories_path}FG/fg_risk_genes_GRCh{assembly}.csv"
        generate_json_from_fg_csv(csv_file, assembly, categories_path)
        generate_bed_from_fg_csv(csv_file, assembly, categories_path)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
    
    
    

    
    
