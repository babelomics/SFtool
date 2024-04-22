# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 15:12:18 2023

@author: kindi
"""
import os
import gzip
import csv
import urllib.request
from datetime import datetime
import shutil

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
