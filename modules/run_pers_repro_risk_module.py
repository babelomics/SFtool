# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 2024

@author: jpflorido
"""
from modules.utils import run_intervar, parse_intervar_output, map_review_status, run_clinvar_filtering, combine_results, write_category_results_to_tsv



def run_pers_repro_risk_module(norm_vcf, assembly, mode, evidence_level, clinvar_db, intervar_path, category):
    """
    Ejecuta el módulo de riesgo personal según el modo seleccionado.

    Args:
        vcf_path (str): Ruta al archivo VCF de entrada.
        assembly (str): Ensamblaje genómico a utilizar.
        mode (str): Modo de ejecución ("basic" o "advanced").
        evidence_level (int): Nivel de evidencia deseado.
        category (str): Categoría de genes para la anotación.
        clinvar_db (str): Ruta al archivo de base de datos de CLINVAR.
    """

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
