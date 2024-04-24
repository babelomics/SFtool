# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 15:12:18 2023

@author: Javier Perez Florido, Edurne Urrutia
"""
import os
import gzip
import csv
import urllib.request
from datetime import datetime
import shutil
from modules.misc.build_json_bed_files import read_csv

def map_review_status(review_status):
    """
    Mapea el estado de revisión (review status) a un nivel de evidencia (evidence level) equivalente.

    Args:
        review_status (str): Estado de revisión en texto.

    Returns:
        int: Nivel de evidencia equivalente en número.
    """
    # Mapear "Review status" a número de estrellas
    mapping = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, conflicting interpretations": 1,
        "criteria provided, single submitter": 1,
        "no assertion for the individual variant": 0,
        "no assertion criteria provided": 0,
        "no assertion provided": 0
    }
    return mapping.get(review_status.lower(), 0)  # Valor predeterminado es 0 si no se encuentra en el mapeo

def run_clinvar(evidence_level, clinvar_db, categories_path, category):
    """
    Run clinvar using the database according to an evidence level

    Args:
        evidence_level (int): Evidence level for variants
        clinvar_db (str): Path to CLINVAR database
        assembly (str): Reference genome version
        categories_path (str): path to directory for categories information
        category (str): either pr or rr

    Returns:
        dict: A dictionary that contains variants from Clinvar and their related information.

    Raises:
        Exception: An exception occurs when an error arises in Clinvar
    """
    try:
        # Select genes for current category

        in_csv = f"{categories_path}{category.upper()}/{category.lower()}_risk_genes.csv"

        # Read CSV and store it in a dictionary
        genes_dct, genes_lst = read_csv(in_csv, category)

        # Read Clinvar database
        clinvar_dct = {}  #  dictionary to store information from CLINVAR

        with open(clinvar_db, "r") as db_file:
            for line in db_file:
                line = line.rstrip()
                if line == "":
                    continue
                fields = line.strip().split("\t")
                gene = fields[2]
                if gene in genes_lst: # Only parse genes for the corresponding category
                    variant = f"{fields[10]}:{fields[15]}:{fields[16]}:{fields[17]}"
                    clinical_significance = fields[3]
                    clinsigsimple = fields[4]
                    rs_id = fields[5]
                    review_status = fields[13]
                    stars = map_review_status(review_status)
                    if stars >= int(evidence_level):
                        clinvar_id = fields[6]
                        clinvar_dct[variant] = {
                            "Gene": gene,
                            "ClinicalSignificance": clinical_significance,
                            "ClinSigSimple": clinsigsimple,
                            "rs": rs_id,
                            "ReviewStatus": '(' + str(stars) + ') ' + review_status,
                            "ClinvarID": clinvar_id
                        }
        return(clinvar_dct)

    except Exception as e:
        print(f"Error when filtering variants: {e}")


def process_clinvar_data(assembly, release_date, clinvar_path):
    """
    Download and process CLINVAR database for a given assembly version
    
    Args:
        assembly (str): Assembly version. Either GRCh37 or GRCh38
        release_date (datetime.datetime): Clinvar release date
        clinvar_path: Path to clinvar directory
    
    Returns:
        str: Output file with processed data
    
    Raises:
        Exception: When an error occurs

    """
    # Columns of interest
    columns_of_interest_names = ["Type", "Name", "GeneSymbol",
                                 "ClinicalSignificance", "ClinSigSimple", "RS# (dbSNP)", "VariationID",
                                 "PhenotypeIDS", "PhenotypeList", "Assembly", 
                                 "Chromosome", "Start", "Stop", "ReviewStatus", 
                                 "SubmitterCategories", "PositionVCF", 
                                 "ReferenceAlleleVCF", "AlternateAlleleVCF"]
    
    # Output file
    output_file = f"{clinvar_path}clinvar_database_{assembly}_{release_date.strftime('%Y%m%d')}.txt"
    
    # Process CLINVAR file
    with gzip.open(f"{clinvar_path}variant_summary.txt.gz", "rt") as gz_file, open(output_file, "w") as output:
        csv_writer = csv.writer(output, delimiter="\t")
        header_line = gz_file.readline().strip()
        header_fields = header_line.split("\t")
        columns_of_interest_positions = [header_fields.index(col) for col in columns_of_interest_names]
        
        # Find the index of the "Assembly" column
        for idx, field in enumerate(header_fields):
            if field == "Assembly":
                assembly_column_index = idx
                break
    
        # Add columns of interest to the output file
        csv_writer.writerow(columns_of_interest_names)
    
        for line in gz_file:
            row = line.strip().split("\t")
            if row[assembly_column_index] == assembly:  # Filter according to genome version (GRCh37 o GRCh38)
                relevant_fields = [row[pos] for pos in columns_of_interest_positions]
                csv_writer.writerow(relevant_fields)
    
    return output_file


def get_clinvar(clinvar_path, assembly):
    """
    Download and process Clinvar database
    
    Args:
        clinvar_path: Path to CLINVAR directory database
    """
    try:        
        # CLINVAR URL
        clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        
        # Open URL
        response = urllib.request.urlopen(clinvar_url)

        # Check whether response is OK (HTTP 200 code)
        if response.status != 200:
            print(f"Error downloading CLINVAR. HTTP code: {response.status}")
            exit(1)
        
        # Open a local file for writing in binary mode
        out_rawfile = f"{clinvar_path}variant_summary.txt.gz"
        with open(out_rawfile, 'wb') as output_file:
            # Copy the response content to the local file
            shutil.copyfileobj(response, output_file)
        print(f"File downloaded to {out_rawfile}")
        
        # Get Clinvar release date
        last_modified = response.headers['Last-Modified']
        if last_modified is None:
            print("Clinvar release date cannot be obtained")
            exit(1)
        
        release_date = datetime.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
        
        # Process CLINVAR file for the assembly
        if assembly == "37":
            clinvar_output_file = process_clinvar_data("GRCh37", release_date, clinvar_path)
            print(f"CLINVAR GRCh37 file is downloaded and processed. Version: {release_date.strftime('%Y%m%d')}")
        else:  # Assembly 38
            clinvar_output_file = process_clinvar_data("GRCh38", release_date, clinvar_path)
            print(f"CLINVAR GRCh38 file is downloaded and processed. Version: {release_date.strftime('%Y%m%d')}")
        
        # Remove donwloaded file
        os.remove(f"{clinvar_path}variant_summary.txt.gz")

        return clinvar_output_file
    
    except Exception as e:
        print(f"Error found: {str(e)}")




def clinvar_manager(clinvar_path, clinvar_ddbb_version, assembly):
    """

    Manage clinvar database: get version contains in the config file or download the latest version

    :param clinvar_path: directory path where Clinvar dabatase is stored
    :param clinvar_ddbb_version: clinvar version. Either a date in the form of YYYYMMDD or 'latest' for downloading the latest one
    :param assembly: reference genome version, either 37 or 38
    :return:
    """


    if clinvar_ddbb_version == "latest": # Download the latest version
        clinvar_db = get_clinvar(clinvar_path, assembly)
    else:  # Use the version contained in the config file
        clinvar_file = os.path.join(clinvar_path, "clinvar_database_GRCh" + str(assembly) + "_" + clinvar_ddbb_version + ".txt")
        if os.path.exists(clinvar_file):
            clinvar_db = clinvar_file
        else:
            print(clinvar_file + " does not exist")
            exit(1)
    return clinvar_db
