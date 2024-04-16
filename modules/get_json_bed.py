# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 19:07:52 2023

@author: kindi
"""

import os
import csv
import json
#from pybedtools import BedTool #DEBERÍA MIRARLO PARA ORDENAR LOS CROMOSOMAS en vez de usar natsorted
from biomart import BiomartServer
from natsort import natsorted

def read_csv(in_csv, category):
    """
    Lee un archivo CSV y almacena la información en un diccionario.
    
    Args:
        in_csv (str): Ruta al archivo CSV.
    
    Returns:
        dict, list: Un diccionario que contiene la información leída y una lista de símbolos de genes.
    """
    # Create dictionary and gene list to store info
    if category == 'pr':
        cat_str = 'personal'
    elif category == 'rr':
        cat_str = 'reproductivo'
    genes_dct = {
    "category": f"Hallazgos secundarios de riesgo {cat_str}",
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
                'variants_to_report': row.get('Variants to Report', '')
                }
            genes_dct["genes"].append(gene_info)
            genes_lst.append(gene_symbol)
    
    return(genes_dct, genes_lst)

def get_gene_pos(gene_symbol, assembly):
    """
    Obtiene la posición de un gen en una versión de ensamblaje específica.
    
    Args:
        gene_symbol (str): Símbolo del gen.
        assembly (str): Versión de ensamblaje ("37" o "38").
    
    Returns:
        dict: Información de posición del gen.
    """
    if assembly == "37":
        server = BiomartServer("http://grch37.ensembl.org/biomart")
    elif assembly == "38":
        server = BiomartServer("http://www.ensembl.org/biomart")
    
    db = server.datasets['hsapiens_gene_ensembl']

    response = db.search({
        'filters': {'external_gene_name': gene_symbol},
        'attributes': ['external_gene_name', 'chromosome_name', 'start_position', 'end_position']
    })
    
    result = {}
    for line in response.iter_lines():
        fields = line.decode('utf-8').split('\t')
        result['Gene_symbol'] = fields[0]
        result['Chromosome'] = fields[1]
        result['Start'] = fields[2]
        result['End'] = fields[3]

    return result

def write_bed_file(assembly, genes_lst, category, categories_path):
    """
    Escribe información de genes en un archivo BED.
    
    Args:
        assembly (str): Versión de ensamblaje ("37" o "38").
        genes_lst (list): Lista de símbolos de genes.
        category (str): Categoría de genes.
        categories_path (str): Ruta al directorio categories.
    
    Returns:
        None
    """
    gene_coords = []
    #### Hay genes que cambian de nombre según el genoma de referencia, pensar como generalizar esto
    for gene in genes_lst:
        if gene == 'MMUT' and assembly == '37':
            gene = 'MUT'
        elif gene == 'ELP1' and assembly == '37':
            gene = 'IKBKAP'
        elif gene == 'G6PC' and assembly == '38':
            gene = 'G6PC1'
        elif gene == 'GBA' and assembly == '38':
            gene = 'GBA1'
        gene_pos = get_gene_pos(gene, assembly)

        if gene_pos['Chromosome'] == 'HG1439_PATCH': #revisar esto, que es una chapuza
            gene_pos['Chromosome'] = 'X'
        elif gene_pos['Chromosome'] == 'HSCHR6_MHC_MCF':
            gene_pos['Chromosome'] = '6'
        
        gene_coords.append((gene_pos['Chromosome'], int(gene_pos['Start']), int(gene_pos['End']), gene))
    sorted_coords = natsorted(gene_coords)
    
    filename = f"{categories_path}{category.upper()}/{category}_risk_genes_GRCh{assembly}.bed"
    with open(filename, "w") as bed_file:
        for chrom, start, end, gene in sorted_coords:
            bed_file.write(f"{chrom}\t{start}\t{end}\t{gene}\n")

def get_json_bed(category, assembly, categories_path):
    """
    Función principal que procesa un archivo CSV y genera archivos JSON y BED.
    
    Args:
        category (str): Categoría de genes.
        assembly (str): Versión de ensamblaje ("37" o "38").
        categories_path (str): Ruta al directorio categories.
        
    Returns:
        None
    """
    
    # CSV infile:
    in_csv = f"{categories_path}{category.upper()}/{category}_risk_genes.csv"
    
    # Read CSV and store it in the dictionary
    genes_dct, genes_lst = read_csv(in_csv, category)
    
    # Write JSON file
    out_json = f"{categories_path}{category.upper()}/{category}_risk_genes.json"
    with open(out_json, 'w') as json_file:
        json.dump(genes_dct, json_file, indent = 4)
    
    # Write BED files   
    write_bed_file(assembly, genes_lst, category, categories_path)