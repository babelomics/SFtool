# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 23:26:14 2023

@author: Edurne Urrutia, Javier Pérez Florido
"""

import json
import pandas as pd
import os
import subprocess


def get_versions_paths(program_arguments, config_data, clinvar_db):
    """
    Get versions of third-party tools from SF tool and program arguments to be shown in the final report
    :param program_arguments: SF tools arguments
    :param config_data: Configuration data
    :param clinvar_db: Clinvar database path
    :return:
    """
    # HPO
    if program_arguments.hpos_file is None:
        hpos_patient_list="Not provided"
    else:
        hpos_list_temp = get_hpos_from_txt(program_arguments.hpos_file)
        hpos_patient_list=",".join(str(hpo) for hpo in hpos_list_temp)

    # Intervar version
    try:
        intervar_file_path = os.path.join(config_data["intervar_path"], "Intervar.py")
        cmd = [intervar_file_path, "--version"]

        intervar_process = subprocess.Popen(cmd, stdout= subprocess.PIPE)
        intervar_out, intervar_err = intervar_process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error running InterVar: {e.output}")

    # Bcftools version
    try:
        bcftools_file_path = os.path.join(config_data["bcftools_path"], "bcftools")
        cmd = [bcftools_file_path, "--version"]

        bcftools_process = subprocess.Popen(cmd, stdout= subprocess.PIPE)
        bcftools_out, bcftools_err = bcftools_process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error running InterVar: {e.output}")



    versions_paths={
        "SF module": "0.1",
        "SF module mode": program_arguments.mode,
        "Input VCF": program_arguments.vcf_file,
        "Temporal dir": config_data["temp_path"],
        "Output dir": config_data["out_path"],
        "HPO list": hpos_patient_list,
        "Human assembly": "hg19" if program_arguments.assembly == 37 else "hg38",
        "Reference genome path":  config_data["reference_genome_37_path"] if program_arguments.assembly == 37 else config_data["reference_genome_38_path"],
        "Clinvar version": "Not used (basic mode)" if program_arguments.mode == "basic" else os.path.splitext(os.path.basename(clinvar_db))[0].split("_")[-1],
        "Clinvar path": clinvar_db,
        "Clinvar Evidence level": str(program_arguments.evidence),
        "Intervar version": str(intervar_out).split(" ")[1] + " " + str(intervar_out).split(" ")[2].split("\\n")[0],
        "Intervar path": config_data["intervar_path"],
        "bcftools version": str(bcftools_out).split(" ")[1].split("\\n")[0],
        "bcftools path": config_data["bcftools_path"],
        "HPO genes to phenotype version": os.path.splitext(os.path.basename(config_data["gene_to_phenotype_file"]))[0].split("_")[-1],
        "HPO genes to phenotype path": config_data["gene_to_phenotype_file"]

    }

    return versions_paths

def combine_variant_and_gene_info(variant_info, gene_info):
    """
    Combina la información de una variante y un gen.

    Args:
        variant_info (dict): Información de la variante.
        gene_info (dict): Información del gen.

    Returns:
        dict: Información combinada de la variante y el gen.
    """
    combined_info = {
        "Gene": variant_info["Gene"],
        "Genotype": variant_info["Genotype"],
        "rs": variant_info.get("rs", ""),
        "IntervarConsequence": variant_info.get("IntervarConsequence", ""),
        "IntervarClassification": variant_info["IntervarClassification"],
        "ClinvarClinicalSignificance": variant_info.get("ClinvarClinicalSignificance", "-"),
        "ReviewStatus": variant_info.get("ReviewStatus", "-"),
        "ClinvarID": variant_info.get("ClinvarID", "-"),
        "Orpha": variant_info.get("Orpha", ""),
        "Phenotype": gene_info["phenotype"],
        "ACMG_version": gene_info.get("ACMG_version", ""),  # Usar get para manejar la falta de 'ACMG_version'
        "OMIM_disorder": gene_info["OMIM_disorder"],
        "inheritance": gene_info["inheritance"],
        "variants_to_report": gene_info.get("variants_to_report", ""),  # Usar get para manejar la falta de 'variants_to_report'
        "related_HPOs": 'NA'
    }
    return combined_info

def check_inheritance(results, category, categories_path):
    """
    Comprueba la herencia de las variantes en función de la categoría de genes y genera un diccionario de variantes
    a informar según las reglas de herencia definidas en el archivo JSON de genes.

    Args:
        results (dict): Un diccionario de resultados de variantes.
        category (str): Categoría de genes para la anotación.

    Returns:
        dict: Un diccionario de variantes a informar siguiendo las reglas de herencia definidas.
    """

    try:
        # Cargar el archivo JSON de categoría de genes
        genes_cat_path = f"{categories_path}{category.upper()}/{category}_risk_genes.json"
        genes_cat = None
        with open(genes_cat_path, "r") as genes_cat_file:
            genes_cat = json.load(genes_cat_file)

        # Crear diccionario de variantes a informar
        reported_variants = {}

        # Recorrer claves del diccionario de resultados combinados
        for variant_key, variant_info in results.items():
            variant_gene = variant_info["Gene"]
            # Obtener el modo de herencia del gen
            for gene in genes_cat['genes']:
                if gene['gene_symbol'] ==  variant_gene:
                    inher = gene["inheritance"]

                    # Si herencia es dominante (AD), semidominante o ligado al X, o si la categoría es RR, se informa la variante
                    if inher in ['AD', 'SD', 'XL'] or category == 'rr':
                        # Combina la información de la variante y el gen
                        combined_info = combine_variant_and_gene_info(variant_info, gene)
                        # Agrega la información combinada a reported_variants
                        reported_variants[variant_key] = combined_info

                    # Si la herencia es recesiva, comprobar genotipo y/o otras variantes en mismo gen
                    elif inher == 'AR':
                        # if gene = 'HFE':
                        #     "variants_to_report": "HFE p.C282Yl\n homozygotes only"
                        # Si la variante está en homocigosis, se reporta
                        if variant_info["Genotype"] == 'hom':
                            # Combina la información de la variante y el gen
                            combined_info = combine_variant_and_gene_info(variant_info, gene)
                            # Agrega la información combinada a reported_variants
                            reported_variants[variant_key] = combined_info

                        # Si la variante está en heterocigosis, sólo se reportará si se encuentra otra variante en el mismo gen
                        elif variant_info["Genotype"] == 'het':
                            # Verifica si hay otra variante en el mismo gen
                            other_variant_in_gene = False
                            for other_variant_key in results:
                                other_variant_info = results[other_variant_key]
                                if (
                                        other_variant_info["Gene"] == variant_gene
                                        and other_variant_key != variant_key
                                ):
                                    other_combined_info = combine_variant_and_gene_info(other_variant_info, gene)
                                    combined_info = combine_variant_and_gene_info(variant_info, gene)

                                    # Agrega la información de ambas variantes a reported_variants
                                    reported_variants[variant_key] = combined_info
                                    reported_variants[other_variant_key] = other_combined_info
        return(reported_variants)
    except Exception as e:
        print(f"Error en check_inheritance: {str(e)}")
        return {}

# def get_hpo_dct(categories_path):

#     hpo_file = f"{categories_path}phenotype_to_genes.txt"
#     hpo_data = {}
#     # Abrir el archivo de texto en modo lectura
#     with open('hpo_file', 'r') as file:
#         # Salta la primera línea si contiene encabezados
#         next(file)

#         # Iterar a través de las líneas del archivo
#         for line in file:
#             # Dividir cada línea en columnas usando tabulaciones o espacios en blanco como delimitadores
#             fields = line.strip().split('\t')

#             # Extraer los valores de cada columna
#             hpo_id, hpo_name, ncbi_gene_id, gene_symbol, disease_id = fields

#             # Comprobar si el HPO ya está en el diccionario, y si no, crea un nuevo diccionario
#             if hpo_id not in hpo_data:
#                 hpo_data[hpo_id] = {
#                     "hpo_name": hpo_name,
#                     "gene_info": []
#                 }

#             # Agrega información del gen y la enfermedad al diccionario correspondiente
#             hpo_data[hpo_id]["gene_info"].append({
#                 "ncbi_gene_id": ncbi_gene_id,
#                 "gene_symbol": gene_symbol,
#                 "disease_id": disease_id
#             })

#     return(hpo_data)

def check_patient_HPO(reported_results, hpos_user, categories_path, gene_to_phenotype_file):
    """
    Comprueba si los HPOs del usuario coinciden con los resultados reportados.

    Args:
        reported_results (dict): Resultados de diagnóstico reportados.
        user_hpos (list): Lista de HPOs proporcionada por el usuario.
        categories_path (str): Ruta al directorio categories.

    Returns:
        list: Lista de resultados coincidentes.
    """
    # Inicializar un diccionario vacío para almacenar los datos
    gene_hpo_dict = {}

    # Abrir el archivo de texto para lectura
    with open(gene_to_phenotype_file, 'r') as file:
        # Leer la primera línea (encabezado) para omitirla
        next(file)

        # Iterar sobre las líneas del archivo
        for line in file:
            # Dividir cada línea en campos separadas por tabulaciones
            fields = line.strip().split('\t')
            gene_symbol = fields[1]
            hpo_id = fields[2]

            # Si el gene_symbol ya está en el diccionario, agregar el hpo_id a la lista existente
            if gene_symbol in gene_hpo_dict:
                gene_hpo_dict[gene_symbol].append(hpo_id)
            # Si el gene_symbol no está en el diccionario, crea una nueva entrada con una lista que contiene el hpo_id
            else:
                gene_hpo_dict[gene_symbol] = [hpo_id]

    for result in reported_results:
        # Comprueba si los hpos asociados al gen coinciden con los hpos del usuario
        gene = reported_results[result]['Gene']
        hpo_results = gene_hpo_dict[gene]

        for hpo in hpo_results:
            if hpo in hpos_user:
                if reported_results[result]['related_HPOs'] == 'NA':
                    reported_results[result]['related_HPOs'] = hpo #habria que mirar que hacer si hay mas de un hpo que coincide
                else:
                    reported_results[result]['related_HPOs'] += ', ' + hpo

    return reported_results

def get_hpos_from_txt(hpos_file):
    """
    Lee un archivo de texto que contiene HPOs y devuelve una lista de HPOs.

    Args:
        hpos_txt (str): Ruta al archivo de texto que contiene los HPOs, uno por línea.

    Returns:
        List[str]: Lista de HPOs leídos del archivo.

    Raises:
        FileNotFoundError: Si el archivo especificado no se encuentra.
    """
    if hpos_file is None:
        return []  # Si no se proporciona un archivo, devolver una lista vacía

    try:
        with open(hpos_file, "r") as file:
            hpos = [line.strip() for line in file.readlines()]
        return hpos

    except FileNotFoundError:
        print(f"El archivo {hpos_file} no se encontró.")
        return []

def generate_report(pr_results, rr_results, fg_results, haplot_results, config_data, args, clinvar_db, categories):
    """
    Write results for all categories to an excel file

    Args:
        pr_results (dict): Results from Personal Risk module
        rr_results (dict): Results from Reproductive Risk module
        fg_results (dict): Results from Pharmacogenetic category module
        out_path (str): Path to output directory
    """
    try:

        categories_path = config_data["categories_path"]
        out_path = config_data["out_path"]
        vcf_file = args.vcf_file
        hpos_file = args.hpos_file
        gene_to_phenotype_file = config_data["gene_to_phenotype_file"]

        # Get versions
        versions_path = get_versions_paths(args, config_data, clinvar_db)

        # Get HPO list
        hpos_user= []
        if versions_path['HPO list'] != "Not provided":
            hpos_user = get_hpos_from_txt(hpos_file)
        outfile = f"{out_path}{vcf_file.split('/')[-1].split('.vcf')[0]}.SF.xlsx"

        with pd.ExcelWriter(outfile) as writer:
            # Write versions and paths
            summary_df = pd.DataFrame.from_dict(versions_path, orient="index")
            summary_df.to_excel(writer, sheet_name= 'Versions and Paths', index=True, header=False)

            for category in categories:
                if category == 'pr':
                    reported_results = check_inheritance(pr_results, category, categories_path)
                    pr_final = check_patient_HPO(reported_results, hpos_user, categories_path, gene_to_phenotype_file)  # pendiente de desarrollar, warning si los términos orpha se corresponden con los hpo del paciente
                    results_df =  pd.DataFrame.from_dict(pr_final, orient='index')
                    results_df.to_excel(writer, sheet_name= category.upper() + ' results', index=True)

                elif category == 'rr':
                    reported_results = check_inheritance(rr_results, category, categories_path)
                    rr_final = check_patient_HPO(reported_results, hpos_user, categories_path, gene_to_phenotype_file)  # pendiente de desarrollar, warning si los términos orpha se corresponden con los hpo del paciente
                    results_df =  pd.DataFrame.from_dict(rr_final, orient='index')
                    results_df.to_excel(writer, sheet_name= category.upper() + ' results', index=True)

                else:
                    fg_df = pd.DataFrame(fg_results)
                    fg_df.to_excel(writer, sheet_name='FG results', index=False)

                    haplot_df = pd.DataFrame(haplot_results)
                    haplot_df.to_excel(writer, sheet_name='FG Diplotype-Phenotype', index=False)

        print(f"Results have been written to '{outfile}'.")

    except Exception as e:
        print(f"Error found when generating final report: {str(e)}")
