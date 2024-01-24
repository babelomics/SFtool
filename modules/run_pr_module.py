# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 22:33:14 2023

@author: kindi
"""
import subprocess
import json
import pysam
import csv
import os

# Definir las funciones para cada módulo y opción
def run_intervar(norm_vcf, category, assembly, intervar_path):
    """
    Ejecuta el programa Intervar para anotar variantes genéticas.
    
    Args:
        norm_vcf (str): Ruta al archivo VCF normalizado.
        category (str): Categoría de genes para la anotación.
        assembly (str): Ensamblaje genómico a utilizar.
        intervar_path (str): Ruta al directorio de InterVar.
    
    Raises:
        subprocess.CalledProcessError: Si hay un error al ejecutar Intervar.
    """
    try:
        # Ruta vcf interseccion y directorio de salida
        input_vcf = f"{norm_vcf.split('normalized')[0]}{category}_intersection.vcf"
        output_file = f"{norm_vcf.split('normalized')[0]}{category}"
        
        if assembly == '37':
            assembly_int = "hg19"
        elif assembly == '38':
            assembly_int = 'hg38'
            
        # Construir el comando para ejecutar Intervar

        intervar_file_path = os.path.join(intervar_path, "Intervar.py")
        cmd = [
            intervar_file_path,
            "-b", assembly_int,    #  no sé si se puede 38 en intervar
            "-i", input_vcf,
            "--input_type", "VCF",
            "-o", output_file
        ]

        # Ejecutar el comando y capturar la salida
        # Cambiar el directorio de trabajo solo para el comando Intervar
        with subprocess.Popen(cmd, stderr=subprocess.STDOUT, text=True, cwd=intervar_path) as process:
        #with subprocess.Popen(cmd, stderr=subprocess.STDOUT, text=True) as process:
            output, _ = process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar InterVar: {e.output}")
        
def parse_intervar_output(norm_vcf, category, mode, assembly):
    """
    Procesa el archivo de salida de InterVar y extrae los campos necesarios.

    Args:
        norm_vcf (str): Ruta al archivo VCF de entrada.
        category (str): Categoría de genes para la anotación.
        mode (str): basic or advanced
        assembly (str): Ensamblaje genómico a utilizar.

    Returns:
        list: Una lista de diccionarios con los campos extraídos.
    """

    if assembly == '37':
        assembly_int = "hg19"
    elif assembly == '38':
        assembly_int = 'hg38'

    intervar_output_file = f"{norm_vcf.split('normalized')[0]}{category}.{assembly_int}_multianno.txt.intervar"
    intervar_results = {}
    
    with open(intervar_output_file, "r") as intervar_file:
        for line in intervar_file:
            # Ignorar líneas que no son variantes
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                # Extraer los campos necesarios y formatearlos como se requiere
                variant = f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
                ref_gene = fields[5]
                intervar_consequence = fields[6] + ',' + fields[7]
                avsnp = fields[9]
                classification = fields[13].split(": ")[1].split(" PVS")[0]
                orpha = fields[32]
                other_info = fields[-1]
                
                # Filtrar solo las variantes patogénicas y posiblemente patogénicas
                if classification in ["Pathogenic", "Likely pathogenic"] or mode == 'advanced':
                
                    # Crear un diccionario con los campos extraídos
                    intervar_results[variant] = {
                        "Gene": ref_gene,
                        "rs": avsnp,
                        "IntervarClassification": classification,
                        "Orpha": orpha,
                        "Genotype": other_info,
                        "IntervarConsequence": intervar_consequence
                    }


    return intervar_results
        
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

def run_clinvar_filtering(evidence_level, clinvar_db, assembly):
    """
    Filtra variantes de la base de datos de CLINVAR según un nivel de evidencia dado.

    Args:
        evidence_level (int): El nivel de evidencia deseado para filtrar las variantes.
        clinvar_db (str): Ruta al archivo de la base de datos de CLINVAR.
        assembly (str):

    Returns:
        dict: Un diccionario que contiene las variantes filtradas de CLINVAR y su información relacionada,
              organizadas por variantes.
    
    Raises:
        Exception: Si ocurre un error durante el filtrado de variantes de CLINVAR.
    """
    try:
        clinvar_path = f"{clinvar_db.split('GRCh')[0]}GRCh{assembly}_{clinvar_db.split('_')[-1]}"

        # Leer la base de datos de CLINVAR
        clinvar_dct = {}  # Un diccionario para almacenar la información de CLINVAR
        
        with open(clinvar_path, "r") as db_file:
            #header = db_file.readline().strip().split("\t")
            #for line in open(clinvar_path):
            for line in db_file:
                line = line.rstrip()
                if line == "":
                    continue
                fields = line.strip().split("\t")
                variant = f"{fields[10]}:{fields[15]}:{fields[16]}:{fields[17]}"
                gene = fields[2]
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
        print(f"Error al filtrar variantes: {e}")
        
        
def combine_results(vcf_norm, category, intervar_results, clinvar_results):
    """
    Combina los resultados de Intervar y ClinVar en una sola línea por variante.

    Args:
        vcf_norm (str): Ruta al archivo VCF normalizado.
        category (str): Categoría de genes para la anotación.
        intervar_results (dict): Resultados de Intervar.
        clinvar_results (dict): Base de datos de ClinVar.

    Returns:
        dict: Un diccionario con los resultados combinados.
    """
    combined_results = {}
    
    # para combinar los resultados de intervar y clinvar, aunque lo ideal sería recorrer
    # las listas, intervar cambia la anotación, eliminando el nt de referencia en las indels
    # por lo que no es posible encontrar esa variante en clinvar (no siempre hay rs disponible)
    # así que hay que recorrer el vcf de interseccion
    


    # Archivo VCF de intersección
    vcf_path = f"{vcf_norm.split('normalized')[0]}{category}_intersection.vcf"
    
    try:
        # Abrir el archivo VCF de intersección
        with open(vcf_path, "r") as vcf_file:
            # Recorrer el archivo línea por línea
            for line in vcf_file:
                fields = line.strip().split("\t")
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]  # Supongo una sola alternativa porque he splitteado los multiallelic sites
        
                # Construye la clave de búsqueda
                variant_key = f"{chrom}:{pos}:{ref}:{alt}"
        
                # Busca la variante en los resultados de Intervar
                # Si es una delecion
                if len(ref) > len(alt):
                    if len(alt) == 1:
                        #variant_int = chrom + ':' + str(int(pos) + 1) + ':' + ref[1:] + ':-'
                        variant_int = f"{chrom}:{str(int(pos) + 1)}:{ref[1:]}:-"
                    else:
                        variant_int = f"{chrom}:{str(int(pos) + 1)}:{ref[1:]}:{alt[1:]}"
                # Si es una inserción
                elif len(ref) < len(alt):
                    if len(ref) == 1:
                        variant_int = f"{chrom}:{pos}:-:{alt[1:]}"
                    else:
                        variant_int = f"{chrom}:{pos}:{ref[1:]}:{alt[1:]}"            
                # Cambio de nt
                else:
                    variant_int = f"{chrom}:{pos}:{ref}:{alt}"
                        
                intervar_info = intervar_results.get(variant_int)
        
                # Busca la variante en los resultados de ClinVar
                clinvar_info = clinvar_results.get(variant_key)

                if clinvar_info is not None:
                    # Combina la información si es "Pathogenic" o "Likely pathogenic" en alguno de los dos
                    if (intervar_info and intervar_info["IntervarClassification"] in ["Pathogenic", "Likely pathogenic"]) or (clinvar_info and clinvar_info["ClinicalSignificance"] in ["Pathogenic", "Likely pathogenic"]) or ((clinvar_info and clinvar_info["ClinicalSignificance"].split(';')[0] == "Conflicting interpretations of pathogenicity") and (clinvar_info["ClinSigSimple"]=="1")):
                            combined_results[variant_key] = {
                                "Gene": clinvar_info["Gene"],
                                "Genotype": intervar_info["Genotype"],
                                "rs": intervar_info["rs"] if intervar_info["rs"] != '.' else clinvar_info["rs"],
                                "IntervarClassification": intervar_info["IntervarClassification"],
                                "ClinvarClinicalSignificance": clinvar_info["ClinicalSignificance"],
                                "ReviewStatus": clinvar_info["ReviewStatus"],
                                "ClinvarID": clinvar_info["ClinvarID"],
                                "Orpha": intervar_info["Orpha"],
                                "IntervarConsequence": intervar_info["IntervarConsequence"]
                            }
                else:
                    # Conserva la info de intervar si no hay info de clinvar
                    if (intervar_info and intervar_info["IntervarClassification"] in ["Pathogenic", "Likely pathogenic"]):
                            combined_results[variant_key] = {
                                "Gene": intervar_info["Gene"],
                                "Genotype": intervar_info["Genotype"],
                                "rs": intervar_info["rs"] if intervar_info["rs"] != '.' else '-',
                                "IntervarClassification": intervar_info["IntervarClassification"],
                                "ClinvarClinicalSignificance": "NA",
                                "ReviewStatus": "NA",
                                "ClinvarID": "NA",
                                "Orpha": intervar_info["Orpha"],
                                "IntervarConsequence": intervar_info["IntervarConsequence"]
                            }

    except Exception as e:
        raise Exception(f"Error al combinar resultados: {e}")
    
    return(combined_results)
    
    
def write_category_results_to_tsv(results, norm_vcf, category):
    """
    Escribe los resultados combinados en un archivo TSV.

    Args:
        results (list): Lista de diccionarios con los resultados combinados.
        norm_vcf (str): Ruta al archivo normalizado.
        category (str): Categoría de genes para la anotación.
    """
    output_tsv = f"{norm_vcf.split('normalized')[0]}{category}_all_results.tsv"
    
    # Abrir el archivo TSV para escritura
    with open(output_tsv, "w", newline="") as tsv_file:
        fieldnames = ["Variant", "Gene", "Genotype", "Intervar Consequence", "rs", "Intervar Classification","Clinvar Clinical Significance", "ReviewStatus", "ClinvarID", "Orpha"]
        writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")
    
        # Escribir el encabezado del archivo TSV
        writer.writeheader()
    
        # Iterar sobre las claves (variantes) en el diccionario result
        for variant, info in results.items():
            # Asegurarse de que el diccionario 'info' tenga todas las claves necesarias
            row = {
                "Variant": variant,
                "Gene": info.get("Gene", ""),
                "Genotype": info.get("Genotype", ""),
                "Intervar Consequence": info.get("IntervarConsequence", "-"),
                "rs": info.get("rs", ""),
                "Intervar Classification": info.get("IntervarClassification", ""),
                "Clinvar Clinical Significance": info.get("ClinvarClinicalSignificance", "-"),
                "ReviewStatus": info.get("ReviewStatus", "-"),
                "ClinvarID": info.get("ClinvarID", "-"),
                "Orpha": info.get("Orpha", "")
            }
    
            # Escribir la fila en el archivo TSV
            writer.writerow(row)    


def run_personal_risk_module(norm_vcf, assembly, mode, evidence_level, clinvar_db, categories_path, intervar_path, temp_path):
    """
    Ejecuta el módulo de riesgo personal según el modo seleccionado.
    
    Args:
        vcf_path (str): Ruta al archivo VCF de entrada.
        assembly (str): Ensamblaje genómico a utilizar.
        mode (str): Modo de ejecución ("basic" o "advanced").
        evidence_level (int): Nivel de evidencia deseado.
        category (str): Categoría de genes para la anotación.
    """
    category = "pr"
    # Intervar is always run
    run_intervar(norm_vcf, category, assembly, intervar_path)
    intervar_results = parse_intervar_output(norm_vcf, category, mode, assembly)
    if mode == "basic":
        category_results = intervar_results
    elif mode == "advanced":
        # Advanced mode: run Clinvar and combine results with Intervar
        clinvar_results = run_clinvar_filtering(evidence_level, clinvar_db, assembly)
        intervar_clinvar_results = combine_results(norm_vcf, category, intervar_results, clinvar_results)
        category_results = intervar_clinvar_results

    # Write results of this category to a file
    write_category_results_to_tsv(category_results, norm_vcf, category)
    return(category_results)




    

