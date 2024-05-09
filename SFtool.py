#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Herramienta para el manejo automático de hallazgos secundarios.

Esta herramienta permite a los usuarios analizar archivos VCF para el manejo automático de hallazgos secundarios relacionados con riesgo personal, riesgo reproductivo y farmacogenético.

@Dependencies InterVar and AnnoVar 
@Usage python3 SFtool.py input_file.vcf --mode <Option: 'basic' or 'advanced'> --evidence <integer> --assembly <Option: '37' or '38'>
@Arguments:
    -vcf (str): Ruta al archivo VCF de entrada.
    -outpath (str): Ruta al directorio donde se guardarán los resultados.
    -mode (str): Modo de análisis (básico o avanzado).
    -evidence (int): Nivel de evidencia de ClinVar para el modo avanzado (1-4).#comprobar que lo he puesto de memoria
    -assembly (str): Ensamblaje genómico a utilizar (GRCh37 o GRCh38).

@Author Javier Perez FLorido, Edurne Urrutia Lafuente
@Date 2023/08/01
@email javier.perez.florido.sspa@juntadeandalucia.es, edurlaf@gmail.com
@github https://github.com/babelomics/secondaryfindings
"""

import os
import json

from modules.misc.arguments import arguments
from modules.misc.build_json_bed_files import build_json_bed_files
from modules.misc.clinvar_utils import clinvar_manager
from modules.FG.run_fg_module import run_pharmacogenomic_risk_module
from modules.PR_RR.run_pers_repro_risk_module import run_pers_repro_risk_module
from modules.misc.check_dependencies import check_dependencies
from modules.misc.vcf_utils import normalize_vcf, intersect_vcf_with_bed
from modules.misc.report_utils import generate_report

def main():


    """
    Get the arguments from command line
    """
    [args, categories] = arguments()
    vcf_file = args.vcf_file
    mode = args.mode
    evidence = args.evidence
    assembly = str(args.assembly)


    """
    Read config file
    """
    # Leer el archivo de configuración config.json
    with open(args.config_file, "r") as config_file:
        config_data = json.load(config_file)
    
    # Get values from config file
    categories_path = config_data["categories_path"]
    clinvar_path = config_data["clinvar_path"]
    clinvar_ddbb_version = config_data["clinvar_ddbb_version"]
    temp_path = config_data["temp_path"]
    out_path = config_data["out_path"]
    intervar_path = config_data["intervar_path"]
    bcftools_path = config_data["bcftools_path"]
    diplotype_phenotype_info_file = config_data["diplotype_phenotype_info_file"]
    personal_risk_geneset_file = config_data["personal_risk_geneset_file"]
    reproductive_risk_geneset_file = config_data["reproductive_risk_geneset_file"]
    pharmacogenetic_risk_variant_GRCh37_file = config_data["pharmacogenetic_risk_variant_GRCh37_file"]


    if assembly == '37':
        reference_genome = config_data["reference_genome_37_path"]
    elif assembly == '38':
        reference_genome = config_data["reference_genome_38_path"]

    """ 
    Create several directories: clinvar, temp and final_output directories (if they dont exist)
    """
    for folder in [clinvar_path, temp_path, out_path]:
        if not os.path.exists(folder):
            os.mkdir(folder)

    """
    Check dependencies
    """
    check_dependencies(intervar_path, bcftools_path)


    """
    Generate JSON and BED files 
    """
    # Check whether BED file (and consequently, JSON file) for each category exist. If not, create them
    # Personal risk catalogue
    if not os.path.exists(f"{categories_path}/PR/pr_risk_genes_GRCh{assembly}.bed"):
        build_json_bed_files("pr", assembly, categories_path, personal_risk_geneset_file)
        
    # Reproductive risk catalogue
    if not os.path.exists(f"{categories_path}/RR/rr_risk_genes_GRCh{assembly}.bed"):
        build_json_bed_files("rr", assembly, categories_path, reproductive_risk_geneset_file)
        
    # Pharma risk catalogue
    if not os.path.exists(f"{categories_path}/FG/fg_risk_genes_GRCh{assembly}.bed"):
        build_json_bed_files("fg", assembly, categories_path, pharmacogenetic_risk_variant_GRCh37_file)
        
    
    """
    In advanced mode, check/update clinVar database
    """
    # If "advanced" mode, check whether Clinvar Database exists
    if mode == 'advanced':
        clinvar_db = clinvar_manager(clinvar_path, clinvar_ddbb_version, assembly)
    else:
        clinvar_db = None


    """
    VCF normalization
    """
    norm_vcf_file = normalize_vcf(vcf_file, temp_path, bcftools_path, reference_genome)
    
    """
    Normalized VCF and BED intersection for each category
    """
    input_vcf_files = {}
    for category in categories:
        category_bed_file = os.path.join(categories_path + category.upper(), category + '_risk_genes_GRCh' + assembly + '.bed')
        generated_vcf_file = intersect_vcf_with_bed(norm_vcf_file, category_bed_file, temp_path, category)
        input_vcf_files[category] = generated_vcf_file
    
    """
    Execute modules selected by the user according to categories
    """
    # Run modules selected by user
    if "pr" in categories:
        # Run Personal Risk (PR) module
        pr_results = run_pers_repro_risk_module(input_vcf_files['pr'], assembly, mode, evidence, clinvar_db, intervar_path, 'pr', personal_risk_geneset_file)
    else:
        pr_results = None
        
    if "rr" in categories:
        # Run Reproductive Risk (RR) module
        rr_results = run_pers_repro_risk_module(input_vcf_files['rr'], assembly, mode, evidence, clinvar_db, intervar_path, 'rr', reproductive_risk_geneset_file)
    else:
        rr_results = None        
        
    if "fg" in categories:
        # Run Pharmacogenetic (FG) module
        fg_results, haplot_results = run_pharmacogenomic_risk_module(categories_path, input_vcf_files['fg'], assembly, temp_path, diplotype_phenotype_info_file)
    else:
        fg_results = None    
        haplot_results = None

    """
    Create report
    """
    generate_report(pr_results, rr_results, fg_results, haplot_results, config_data, args, clinvar_db, categories)


    
        
if __name__ == "__main__":
    main()
