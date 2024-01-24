# -*- coding: utf-8 -*-

"""
Herramienta para el manejo automático de hallazgos secundarios.

Esta herramienta permite a los usuarios analizar archivos VCF para el manejo automático de hallazgos secundarios relacionados con riesgo personal, riesgo reproductivo y farmacogenético.

@Dependencies InterVar and AnnoVar 
@Usage python3 secondary_findings.py input_file.vcf --mode <Option: 'basic' or 'advanced'> --evidence <integer> --assembly <Option: '37' or '38'> 
@Arguments:
    -vcf (str): Ruta al archivo VCF de entrada.
    -outpath (str): Ruta al directorio donde se guardarán los resultados.
    -mode (str): Modo de análisis (básico o avanzado).
    -evidence (int): Nivel de evidencia de ClinVar para el modo avanzado (1-4).#comprobar que lo he puesto de memoria
    -assembly (str): Ensamblaje genómico a utilizar (GRCh37 o GRCh38).

@Author Edurne Urrutia Lafuente
@Date 2023/08/01
@email edurlaf@gmail.com
@github github.com/edurlaf
"""

#import sys
import os
import json

from modules.arguments import arguments
from modules.get_json_bed import read_csv, get_gene_pos, write_bed_file, get_json_bed
from modules.get_json_bed_fg import generate_json_from_fg_csv, generate_bed_from_fg_csv, get_json_bed_fg
from modules.get_clinvar import process_clinvar_data, get_clinvar
from modules.normalize_vcf import normalize_vcf
from modules.intersect_vcf_bed import intersect_vcf_with_bed
from modules.run_pr_module import run_intervar, parse_intervar_output, map_review_status, run_clinvar_filtering, combine_results, write_category_results_to_tsv, run_personal_risk_module
from modules.run_rr_module import run_reproductive_risk_module, run_intervar
from modules.run_fg_module import annotate_fg_variants, check_gene_variants, assign_cyp2c9_diplotype, assign_cyp2c19_diplotype, assign_dpyd_diplotype, assign_nudt15_diplotype, assign_tpmt_diplotype, get_diplotype_phenotype_dictionary, run_pharmacogenomic_risk_module
from modules.generate_report import combine_variant_and_gene_info, check_inheritance, check_patient_HPO, get_hpos_from_txt, generate_report

def main():


    """
    Get the arguments
    """
    args = arguments()

    """
    Read config file
    """
    # Leer el archivo de configuración config.json
    with open(args.config_file, "r") as config_file:
        config_data = json.load(config_file)
    
    # Obtener los valores del archivo de configuración
    categories_path = config_data["categories_path"]
    clinvar_path = config_data["clinvar_path"]
    temp_path = config_data["temp_path"]
    out_path = config_data["out_path"]
    intervar_path = config_data["intervar_path"]
    bcftools_path = config_data["bcftools_path"]
    reference_genome_37_path = config_data["reference_genome_37_path"]
    gene_to_phenotype_file = config_data["gene_to_phenotype_file"]

    """ 
    Create clinvar, temp and final_output directories
    """    
    # Comprobar si los directorios 'clinvar', 'temp' y 'final_output' existen y crearlos si no
    for folder in [clinvar_path, temp_path, out_path]:
        if not os.path.exists(folder):
            os.mkdir(folder)       
    

    
    # Argumentos del usuario
    vcf_file = args.vcf_file
    mode = args.mode
    evidence = args.evidence
    assembly = str(args.assembly)
    hpos_file = args.hpos_file
    
    # Obtener la preferencia del usuario para las categorías a analizar (PR, RR, FG)
    categories_usr = input("Elija las categorías a analizar (PR, RR, FG separados por comas): ") ####cambiar al config
    categories = [category.strip().lower() for category in categories_usr.split(",")]
    
    
    """
    Generate JSON and BED files 
    """
    # Comprobar si los archivos JSON existen
    # Catálogo de riesgo personal
    if not os.path.exists(f"{categories_path}/PR/pr_risk_genes.json"):
        print("Generando archivos JSON y BED para riesgo personal.") # mejor dentro de la función get_json
        get_json_bed("pr", assembly, categories_path)
        
    # Catálogo de riesgo reproductivo
    if not os.path.exists(f"{categories_path}/RR/rr_risk_genes.json"):
        print("Generando archivos JSON y BED para riesgo reproductivo.") # mejor dentro de la función get_json
        get_json_bed("rr", assembly, categories_path)
        
    # Catálogo de riesgo farmacogenético
    if not os.path.exists(f"{categories_path}/FG/fg_risk_variants_grch{assembly}.json"):
        print("Generando archivos JSON y BED para riesgo farmacogenético.") # mejor dentro de la función get_json
        get_json_bed_fg(assembly, categories_path)
        
    
    """
    Get ClinVar database
    """
    # Si el modo es avanzado, comprobar si se ha descargado la BD ClinVar
    if mode == 'advanced':
        clinvar_files = [file for file in os.listdir(clinvar_path) if file.startswith("clinvar_database_")]
        
        # Si hay archivos clinvar, seleccionar el más reciente
        if clinvar_files:
            clinvar_files.sort(reverse=True)
            last_clinvar = clinvar_files[0]

            # Obtener la versión del nombre del archivo
            last_version = last_clinvar.split('_')[3].split('.')[0]      
            #print(f"El archivo ClinVar más reciente encontrado es {archivo_mas_reciente}.")
            print(f"La versión actual del archivo ClinVar es {last_version}.")
        
            # Preguntar al usuario si desea actualizar el archivo
            answr = input("¿Deseas actualizarlo? (S/N): ")       
            if answr.lower() == "s":
                # Descargar el archivo actualizado
                clinvar_db = get_clinvar(clinvar_path)
            else:
                print("No se actualizará el archivo ClinVar.")
                clinvar_db = f"{clinvar_path}{last_clinvar}"
        # Si no se encuentran archivos ClinVar, descargarlo y guardarlo
        else:
            print("No se encontraron archivos ClinVar en el directorio.")
            print("El archivo ClinVar se descargará.")
            clinvar_db = get_clinvar(clinvar_path)
    else:
        clinvar_db = None

    """
    Comprobar dependencias
    """
    # Comprobar si InterVar está en el path
    if not os.path.exists(intervar_path):
        print("InterVar no está instalado. Por favor, instálalo para continuar.")

    # Comprobar si bcftools está en el path : TO DO
    # Comprobar si Annovar está en el path: TO DO
    
    """
    Normalizar VCF de entrada
    """
    norm_vcf = normalize_vcf(vcf_file, temp_path, assembly, bcftools_path, reference_genome_37_path)
    
    """
    Realizar la intersección con los archivos BED
    """
    for category in categories:
        intersect_vcf_with_bed(norm_vcf, category, assembly, categories_path)
    
    """
    Ejecutar los módulos que correspondan:
    """
    # Verificar y ejecutar los módulos elegidos por el usuario
    if "pr" in categories:
        # Ejecutar el módulo de riesgo personal (PR)
        print("Ejecutando módulo de riesgo personal...")
        pr_results = run_personal_risk_module(norm_vcf, assembly, mode, evidence, clinvar_db, categories_path, intervar_path, temp_path)
    else:
        pr_results = None
        
    if "rr" in categories:
        # Ejecutar el módulo de riesgo reproductivo (RR)
        print("Ejecutando módulo de riesgo reproductivo...")
        rr_results = run_reproductive_risk_module(norm_vcf, assembly, mode, evidence, clinvar_db, categories_path, intervar_path)
    else:
        rr_results = None        
        
    if "fg" in categories:
        # Ejecutar el módulo farmacogenético (FG)
        print("Ejecutando módulo farmacogenético...")
        fg_results, haplot_results = run_pharmacogenomic_risk_module(categories_path, norm_vcf, assembly, temp_path)
    else:
        fg_results = None    
        haplot_results = None
    # Informar al usuario que los módulos han sido ejecutados
    print("Módulos de análisis completados.")
    
    
    """
    Generar el informe de salida
    """
    out_file = generate_report(pr_results, rr_results, fg_results, haplot_results, categories_path, out_path, categories, vcf_file, hpos_file, gene_to_phenotype_file)
    print("Informe de resultados generado.\n ---Finalizado---")

    
        
if __name__ == "__main__":
    main()
