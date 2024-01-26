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
    Descarga y procesa los datos de CLINVAR para una versión de ensamblaje específica.
    
    Args:
        assembly (str): La versión de ensamblaje para la cual se desean los datos (por ejemplo, 'GRCh37' o 'GRCh38').
        release_date (datetime.datetime): La fecha de lanzamiento de los datos de CLINVAR.
        clinvar_path: Ruta al directorio clinvar.
    
    Returns:
        str: El nombre del archivo de salida que contiene los datos procesados.
    
    Raises:
        Exception: Si ocurre un error durante el procesamiento de los datos.

    """
    # Nombres de las columnas de interés
    columns_of_interest_names = ["Type", "Name", "GeneSymbol",
                                 "ClinicalSignificance", "ClinSigSimple", "RS# (dbSNP)", "VariationID",
                                 "PhenotypeIDS", "PhenotypeList", "Assembly", 
                                 "Chromosome", "Start", "Stop", "ReviewStatus", 
                                 "SubmitterCategories", "PositionVCF", 
                                 "ReferenceAlleleVCF", "AlternateAlleleVCF"]
    
    # Nombre del archivo de salida
    output_file = f"{clinvar_path}clinvar_database_{assembly}_{release_date.strftime('%Y%m%d')}.txt"
    
    # Procesar el archivo CLINVAR
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
    
        # Agregar las columnas de interés al archivo de salida
        csv_writer.writerow(columns_of_interest_names)
    
        for line in gz_file:
            row = line.strip().split("\t")
            if row[assembly_column_index] == assembly:  # Filtro por ensamblaje (GRCh37 o GRCh38)
                relevant_fields = [row[pos] for pos in columns_of_interest_positions]
                csv_writer.writerow(relevant_fields)
    
    return output_file

def get_clinvar(clinvar_path):
    """
    Descarga y procesa los datos de la base de datos CLINVAR.
    
    Args:
        clinvar_path: Ruta al directorio clinvar.
    """
    try:        
        # URL del archivo CLINVAR
        clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        
        # Abrir la URL
        response = urllib.request.urlopen(clinvar_url)

        # Verificar si la respuesta fue exitosa (código de estado HTTP 200)
        if response.status != 200:
            print(f"Error al descargar el archivo CLINVAR. Código de estado HTTP: {response.status}")
            return
        
        # Open a local file for writing in binary mode
        out_rawfile = f"{clinvar_path}variant_summary.txt.gz"
        with open(out_rawfile, 'wb') as output_file:
            # Copy the response content to the local file
            shutil.copyfileobj(response, output_file)
        print(f"File downloaded to {out_rawfile}")
        
        # Obtener la fecha de release del archivo CLINVAR desde la URL
        last_modified = response.headers['Last-Modified']
        if last_modified is None:
            print("No se pudo obtener la fecha de release del archivo CLINVAR.")
            return
        
        release_date = datetime.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
        
        # Procesar el archivo CLINVAR para GRCh37
        grch37_output_file = process_clinvar_data("GRCh37", release_date, clinvar_path)
        print(f"Archivo CLINVAR GRCh37 descargado y procesado. Versión: {release_date.strftime('%Y%m%d')}")
        
        # Procesar el archivo CLINVAR para GRCh38
        grch38_output_file = process_clinvar_data("GRCh38", release_date, clinvar_path)
        print(f"Archivo CLINVAR GRCh38 descargado y procesado. Versión: {release_date.strftime('%Y%m%d')}")
        
        # Eliminar el archivo descargado en formato gz
        os.remove(f"{clinvar_path}variant_summary.txt.gz")
        
        # # Eliminar versiones anteriores si existen
        # for filename in os.listdir(clinvar_path):
        #     if (filename.startswith("clinvar_database_") and filename.endswith(".txt")) and filename != grch37_output_file and filename != grch38_output_file:
        #         os.remove(filename)
        return(grch37_output_file) #si al final eintervar tamb funciona con genoma 38, cambiar esto en función del assembly
    
    except Exception as e:
        print(f"Ocurrió un error: {str(e)}")