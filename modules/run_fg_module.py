# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 11:21:45 2023

@author: kindi
"""
import subprocess
import json
import pysam
import csv

def annotate_fg_variants(categories_path, norm_vcf, assembly, temp_path):
    """
    Anota variantes genéticas utilizando un archivo JSON de variantes farmacogenéticas.
    
    Args:
        categories_path (str): Ruta al directorio que contiene archivos relacionados con las categorías.
        norm_vcf (str): Ruta al archivo VCF normalizado.
        assembly (str): Ensamblaje genómico a utilizar.
        temp_path (str): Ruta al directorio que contiene archivos intermedios.
    
    Returns:
        list: Una lista de diccionarios que contienen los resultados de la anotación.
    
    Raises:
        FileNotFoundError: Si no se puede encontrar el archivo JSON con variantes farmacogenéticas.
    """

    try:
        # Cargar archivo fg_json con variantes farmacogenéticas
        fg_json_path = f'{categories_path}FG/fg_risk_genes_GRCh{assembly}.json'
        with open(fg_json_path, 'r') as file:
            fg_json = json.load(file)

        annotated_variants = []
    
        # Abrir el archivo VCF
        vcf_int_path = norm_vcf.split('_normalized.vcf')[0] + '_fg_intersection.vcf'
        
        with open(vcf_int_path, 'r') as vcf_file:
        #with pysam.VariantFile(vcf_path) as vcf_file:
            for line in vcf_file:
                fields = line.split('\t')
                variant_key = f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
                #variant_key = f"{record.chrom}:{record.pos}:{record.ref}:{record.alts[0]}"
                variant_info = fg_json['variants'].get(variant_key, None)
    
                if variant_info:
                    annotated_variant = {
                        "Variant": variant_key,
                        "GT": fields[9].split(':')[0],
                        "Gene": variant_info["gene_symbol"],
                        "rs": variant_info["rs"]
                    }
                else:
                    annotated_variant = {
                        "Variant": variant_key,
                        "GT": fields[9].split(':')[0],
                        "Gene": "No encontrado",
                        "rs": "No encontrado",
                    }
    
                annotated_variants.append(annotated_variant)
                
        # Escribir los resultados en un archivo
        result_file = f'{temp_path}/annotated_variants.txt'
        with open(result_file, 'w') as fout:
            json.dump(annotated_variants, fout)
                
        return annotated_variants
    
    except FileNotFoundError as e:
       print(f"Error: {e}")
       return []

def check_gene_variants(variants, gene):
    """
    Verifica la presencia de variantes en un gen específico y recopila información sobre su genotipo.
    
    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        gene (str): El nombre del gen objetivo a buscar en las variantes.
    
    Returns:
        tuple: Una tupla que contiene dos elementos:
            - found (bool): True si se encontró al menos una variante en el gen objetivo, False en caso contrario.
            - variants_gene (dict): Un diccionario que almacena el genotipo de las variantes encontradas en el gen, donde
              la clave es el rsID de la variante y el valor es el genotipo.
    
    Raises:
        TypeError: Si 'variants' no es una lista o si 'gene' no es una cadena de caracteres.
    """
    # Control de errores para los argumentos de entrada
    if not isinstance(variants, list):
        raise TypeError("El argumento 'variants' debe ser una lista de variantes genéticas.")
    if not isinstance(gene, str):
        raise TypeError("El argumento 'gene' debe ser una cadena de caracteres que representa el gen objetivo.")

    # Flag para verificar si se encontró una variante en el gen objetivo
    found = False
    # Crear diccionario para almacenar genotipo de variantes del gene       #igual mejor guardar la variante también, para el caso de indels con mismo rs y diferente nº de repeticiones
    variants_gene = {}
    
    # Iterar a través de las variantes
    for variant in variants:
        if variant['Gene'] == gene:
            found = True
            variants_gene[variant['rs']] = variant['GT']
    
    return found, variants_gene

def assign_cyp2c9_diplotype(variants, diplo_pheno_dct, results):
    """
    Asigna un diplotipo a CYP2C9 basado en las variantes genéticas y almacena los resultados.
    
    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        results (list): Una lista de resultados donde se agregarán los resultados de la asignación.
    
    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de CYP2C9.
    
    Raises:
        TypeError: Si 'variants' no es una lista, si 'diplo_pheno_dct' no es un diccionario o si 'results' no es una lista.
    """
    # Control de errores para los argumentos de entrada
    if not isinstance(variants, list):
        raise TypeError("El argumento 'variants' debe ser una lista de variantes genéticas.")
    if not isinstance(diplo_pheno_dct, dict):
        raise TypeError("El argumento 'diplo_pheno_dct' debe ser un diccionario que contiene información sobre diplotipos y fenotipos genéticos.")
    if not isinstance(results, list):
        raise TypeError("El argumento 'results' debe ser una lista para almacenar los resultados de la asignación de diplotipos.")

    # Gen
    gene = 'CYP2C9'
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'
        
    else:
        if 'rs1799853' in variants_gene.keys():
            if variants_gene['rs1799853'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2'
                else:
                    if 'rs1057910' in variants_gene.keys():
                        diplotype = '*2/*3'
                    else:
                        diplotype = 'NA'
                        print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        elif 'rs1057910' in variants_gene.keys():
            if variants_gene['rs1057910'] == '1/1':
                diplotype = '*3/*3'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3'
                else:
                    diplotype = 'NA'
                    print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        else:
            diplotype = 'NA'
            print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

    # Obtener el fenotipo y AS del diccionario               # esto se repite en todas las assign_diplotype, meter en una función   #contemplar qué hacer cuando no hay diplotipo asignado
    if gene in diplo_pheno_dct and diplotype in diplo_pheno_dct[gene]:
        data = diplo_pheno_dct[gene][diplotype]
        phenotype = data['Phenotype']
        activity_score = data['Activity_Score']
    else:
        phenotype = 'NA'
        activity_score = 'NA'

    # Agregar los resultados a la lista de diccionarios
    results.append({
        'Gene': gene,
        'Diplotipo': diplotype,
        'Phenotype': phenotype,
        'Activity Score': activity_score
    })

    return results              
        
def assign_cyp2c19_diplotype(variants, diplo_pheno_dct, results):
    """
    Asigna un diplotipo a CYP2C19 basado en las variantes genéticas y almacena los resultados.
    
    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        results (list): Una lista de resultados donde se agregarán los resultados de la asignación.
    
    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de CYP2C19.

    """
    # Gen
    gene = 'CYP2C19'
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False or (len(variants_gene) == 1 and 'rs3758581' in variants_gene.keys()):
        diplotype = '*1/*1'
        
    else:
        if 'rs12769205' in variants_gene.keys() and 'rs4244285' in variants_gene.keys():
            if variants_gene['rs12769205'] == '1/1'  and variants_gene['rs4244285'] == '1/1': # en realidad habría que comprobar si es len >2, porque puede que esté en fase con otro alelo
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 2 or (len(variants_gene) == 3 and 'rs3758581' in variants_gene.keys()):
                    diplotype = '*1/*2'
                else:
                    if 'rs4986893' in variants_gene.keys():
                        diplotype = '*2/*3'
                    elif 'rs12248560' in variants_gene.keys():
                        diplotype = '*2/*17'
                    else:
                        print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                        
        elif 'rs4986893' in variants_gene.keys():
            if variants_gene['rs4986893'] == '1/1':
                diplotype = '*3/*3'
            else:
                if len(variants_gene) == 1 or (len(variants_gene) == 2 and 'rs3758581' in variants_gene.keys()):
                    diplotype = '*1/*3'
                elif 'rs12248560' in variants_gene.keys():
                    diplotype = '*3/*17'
                else:
                    print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                    
        # me falta el *17
        elif 'rs12248560' in variants_gene.keys():
            if variants_gene['rs12248560'] == '1/1':
                diplotype = '*17/*17'
            elif len(variants_gene) == 1 or (len(variants_gene) == 2 and 'rs3758581' in variants_gene.keys()):
                diplotype = '*1/*17'
            else:
                print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                
    # Obtener el fenotipo y AS del diccionario
    if gene in diplo_pheno_dct and diplotype in diplo_pheno_dct[gene]:
        data = diplo_pheno_dct[gene][diplotype]
        phenotype = data['Phenotype']
        activity_score = data['Activity_Score']

    # Agregar los resultados a la lista de diccionarios
    results.append({
        'Gene': gene,
        'Diplotipo': diplotype,
        'Phenotype': phenotype,
        'Activity Score': activity_score
    })

    return results   

def assign_dpyd_diplotype(variants, diplo_pheno_dct, results):
    """
    Asigna un diplotipo a DPYD basado en las variantes genéticas y almacena los resultados.
    
    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        results (list): Una lista de resultados donde se agregarán los resultados de la asignación.
    
    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de DPYD.
    """
    # Gen
    gene = 'DPYD'
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'
        
    else:
        # Chequear *2A:
        if 'rs3918290' in variants_gene.keys():
            if variants_gene['rs3918290'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                diplotype = '*2A/*2A'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2A'
                elif len(variants_gene) == 2:
                    if 'rs55886062' in variants_gene.keys():
                        diplotype = '*2A/*13'
                    elif 'rs67376798' in variants_gene.keys():
                        diplotype = '*2A/c.2846A>T'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = '*2a/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *13:
        elif 'rs55886062' in variants_gene.keys():
            if variants_gene['rs55886062'] == '1/1':
                diplotype = '*13/*13'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*13'
                elif len(variants_gene) == 2:
                    if 'rs67376798' in variants_gene.keys():
                        diplotype = '*13/c.2846A>T'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = '*13/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        # Chequear c.2846A>T:
        elif 'rs67376798' in variants_gene.keys():
            if variants_gene['rs67376798'] == '1/1':
                diplotype = 'c.2846A>T/c.2846A>T'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/c.2846A>T'
                elif len(variants_gene) == 2:
                    print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = 'c.2846A>T/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        # Chequear c.2846A>T:
        elif 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
            if variants_gene['rs75017182'] == '1/1' and variants_gene['rs56038477'] == '1/1':
                diplotype = 'HapB3/HapB3'
            else:
                if len(variants_gene) == 2:
                    diplotype = '*1/HapB3'
                elif len(variants_gene) == 3:
                    print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
    # Obtener el fenotipo y AS del diccionario
    if gene in diplo_pheno_dct and diplotype in diplo_pheno_dct[gene]:
        data = diplo_pheno_dct[gene][diplotype]
        phenotype = data['Phenotype']
        activity_score = data['Activity_Score']

    # Agregar los resultados a la lista de diccionarios
    results.append({
        'Gene': gene,
        'Diplotipo': diplotype,
        'Phenotype': phenotype,
        'Activity Score': activity_score
    })

    return results    

def assign_nudt15_diplotype(variants, diplo_pheno_dct, results):
    """
    Asigna un diplotipo a NUDT15 basado en las variantes genéticas y almacena los resultados.
    
    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        results (list): Una lista de resultados donde se agregarán los resultados de la asignación.
    
    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de NUDT15.
    """
    # Gen
    gene = 'NUDT15'
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'
        
    else:
        # Chequear *2A:
        if 'rs116855232' in variants_gene.keys():
            if len(variants_gene) == 1:
                if variants_gene['rs116855232'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                    diplotype = '*3/*3'
                else:
                    diplotype = '*1/*3'
            elif len(variants_gene) == 2:
                if 'rs746071566' in variants_gene.keys():
                    if variants_gene['rs746071566'] == '1/1' and variants_gene['rs116855232'] == '1/1':
                        diplotype = '*2/*2'
                    elif variants_gene['rs746071566'] == '0/1' and variants_gene['rs116855232'] == '1/1': #hbaría que chequear variante si es inserción o del
                        diplotype = '*2/*3'
                    elif variants_gene['rs746071566'] == '0/1' and variants_gene['rs116855232'] == '0/1':
                        diplotype = '*1/*2'
                    else:
                        print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
            else:
                print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        else:
            print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
    # Obtener el fenotipo y AS del diccionario
    if gene in diplo_pheno_dct and diplotype in diplo_pheno_dct[gene]:
        data = diplo_pheno_dct[gene][diplotype]
        phenotype = data['Phenotype']
        activity_score = data['Activity_Score']

    # Agregar los resultados a la lista de diccionarios
    results.append({
        'Gene': gene,
        'Diplotipo': diplotype,
        'Phenotype': phenotype,
        'Activity Score': activity_score
    })

    return results   
 
def assign_tpmt_diplotype(variants, diplo_pheno_dct, results):
    """
    Asigna un diplotipo a TPMT basado en las variantes genéticas y almacena los resultados.
    
    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        results (list): Una lista de resultados donde se agregarán los resultados de la asignación.
    
    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de TPMT.
    """
    # Gen
    gene = 'TPMT'
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'
        
    else:
        # Chequear *2:
        if 'rs1800462' in variants_gene.keys():
            if variants_gene['rs1800462'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2'
                elif len(variants_gene) == 2:
                    if 'rs1800460' in variants_gene.keys():
                        diplotype = '*2/*3B'
                    elif 'rs1142345' in variants_gene.keys():
                        diplotype = '*2/*3C'
                    elif 'rs1800584' in variants_gene.keys():
                        diplotype = '*2/*4'    
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*2/*3A'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *3B:
        elif 'rs1800460' in variants_gene.keys():
            if variants_gene['rs1800460'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*3B/*3B'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3B'
                elif len(variants_gene) == 2:
                    if 'rs1142345' in variants_gene.keys() and (variants_gene['rs1800460'] == '0/1') and (variants_gene['rs1142345'] == '0/1'):
                        diplotype = '*1/*3A o *3B/*3C'
                    elif 'rs1800584' in variants_gene.keys():
                        diplotype = '*3B/*4' 
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*3B'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                    
        # Chequear *3C:
        elif 'rs1142345' in variants_gene.keys():
            if variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*3C/*3C'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3C'
                elif len(variants_gene) == 2:
                    if 'rs1800584' in variants_gene.keys():
                        diplotype = '*3C/*4'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*3B'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *4:
        elif 'rs1800584' in variants_gene.keys():
            if variants_gene['rs1800584'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*4/*4'
            else:
                if len(variants_gene) == 1 and variants_gene['rs1800584'] == '0/1':
                    diplotype = '*1/*4'
                elif len(variants_gene) == 2:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*4'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                    
        # Chequear *3A:  IGUAL REORGANIZO TODO TPMT Y ESTO AL INICIO
        elif 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
            if variants_gene['rs1800460'] == '1/1' and variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3A'
            elif variants_gene['rs1800460'] == '1/1' and variants_gene['rs1142345'] == '0/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3B'
            elif variants_gene['rs1800460'] == '0/1' and variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3C'
            elif variants_gene['rs1800460'] == '0/1' and variants_gene['rs1142345'] == '0/1' and len(variants_gene) == 2:
                diplotype = '*1/*3A o *3A/*3C'
            else:
                print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
    # Obtener el fenotipo y AS del diccionario
    if gene in diplo_pheno_dct and diplotype in diplo_pheno_dct[gene]:
        data = diplo_pheno_dct[gene][diplotype]
        phenotype = data['Phenotype']
        activity_score = data['Activity_Score']

    # Agregar los resultados a la lista de diccionarios
    results.append({
        'Gene': gene,
        'Diplotipo': diplotype,
        'Phenotype': phenotype,
        'Activity Score': activity_score
    })

    return results    
       
def get_diplotype_phenotype_dictionary(categories_path):
    """
    Obtener un diccionario de diplotipos con fenotipos y puntuaciones de actividad.
    
    Args:
        categories_path (str): Ruta al directorio categories.
    
    Returns:
        dict: Un diccionario que asigna diplotipos con su información de fenotipo y puntuación de actividad.
    """
    # Nombre del archivo CSV
    csv_file = f'{categories_path}FG/diplotipo_fenotipo.csv'
    
    # Definir un diccionario para mapear los diplotipos con su fenotipo y Activity Score
    diplotype_data = {}
    
    try:
        # Abrir el archivo CSV y cargar los datos en el diccionario
        with open(csv_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gene = row["GENE"]
                diplotype = row["DIPLOTYPE"]
                phenotype = row["Phenotype"]
                activity_score = row["Activity_Score"]
                # Verificar si el gen ya está en el diccionario
                if gene in diplotype_data:
                    # Si el gen ya está en el diccionario, agregar el nuevo diplotipo
                    diplotype_data[gene][diplotype] = {"Phenotype": phenotype, "Activity_Score": activity_score}
                else:
                    # Si el gen no está en el diccionario, crear una entrada para el gen y agregar el diplotipo
                    diplotype_data[gene] = {diplotype: {"Phenotype": phenotype, "Activity_Score": activity_score}}
    except FileNotFoundError:
        raise FileNotFoundError(f"El archivo CSV no se encuentra en la ruta especificada: {csv_file}")

    return(diplotype_data)

    
def run_pharmacogenomic_risk_module(categories_path, norm_vcf, assembly, temp_path): #sobra category, este modulo es especifico de rr
    """
    Ejecuta el módulo de riesgo farmacogenético.
    
    Args:
        categories_path (str): Ruta al directorio categories.
        norm_vcf (str): Ruta al archivo VCF dnormalizado.
        assembly (str): Ensamblaje genómico a utilizar.
        temp_path (str): Rutal al directorio de archivos intermedios.
        
    Returns:
        list: Una lista de diccionarios que contienen los resultados de los genes procesados.
    """
    # Anotar variantes fg presentes en el vcf
    fg_results = annotate_fg_variants(categories_path, norm_vcf, assembly, temp_path)
    
    # Crear diccionario con asociaciones diplotipo-fenotipo
    diplo_pheno_dct = get_diplotype_phenotype_dictionary(categories_path)
    
    # Asignar los diplotipos de cada gen
    results = []
    results = assign_cyp2c9_diplotype(fg_results, diplo_pheno_dct, results)
    results = assign_cyp2c19_diplotype(fg_results, diplo_pheno_dct, results)
    results = assign_dpyd_diplotype(fg_results, diplo_pheno_dct, results)
    results = assign_nudt15_diplotype(fg_results, diplo_pheno_dct, results)
    results = assign_tpmt_diplotype(fg_results, diplo_pheno_dct, results)
    
    return(fg_results, results)
    

