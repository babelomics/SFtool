# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 11:21:45 2023

@author: jpflorido
"""
from modules.FG.diplotype_scripts.cyp2c9_diplotype_algorithm import assign_cyp2c9_diplotype
from modules.FG.diplotype_scripts.cyp2c19_diplotype_algorithm import assign_cyp2c19_diplotype
from modules.FG.diplotype_scripts.dpyd_diplotype_algorithm import assign_dpyd_diplotype
from modules.FG.diplotype_scripts.nudt15_diplotype_algorithm import assign_nudt15_diplotype
from modules.FG.diplotype_scripts.tpmt_diplotype_algorithm import assign_tpmt_diplotype

from modules.FG.misc_fg import annotate_fg_variants, get_diplotype_phenotype_dictionary

    
def run_pharmacogenomic_risk_module(categories_path, norm_vcf, assembly, temp_path):
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
    fg_annotated_variants = annotate_fg_variants(categories_path, norm_vcf, assembly, temp_path)
    
    # Crear diccionario con asociaciones diplotipo-fenotipo
    diplo_pheno_info = get_diplotype_phenotype_dictionary(categories_path)
    
    # Asignar los diplotipos de cada gen
    aggregated_fg_results = []
    aggregated_fg_results = assign_cyp2c9_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    aggregated_fg_results = assign_cyp2c19_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    aggregated_fg_results = assign_dpyd_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    aggregated_fg_results = assign_nudt15_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    aggregated_fg_results = assign_tpmt_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    
    return(fg_annotated_variants, aggregated_fg_results)
    

