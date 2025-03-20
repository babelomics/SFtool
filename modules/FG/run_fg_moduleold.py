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
from modules.FG.misc_fg import write_fg_results_to_tsv

    
def run_pharmacogenomic_risk_module(categories_path, norm_vcf, assembly, temp_path, diplotype_phenotype_info_file):
    """
    Pharma risk module execution
    
    Args:
        categories_path (str): Path to categories directory
        norm_vcf (str): Path to normalized and intersected VCF file for this category
        assembly (str): Assembly version.
        temp_path (str): Path to temporary directory
        diplotype_phenotype_info_file (str): Path to diplotype_phenotype info file
        
    Returns:
        list: A list of dictionary with pharma results
    """
    # Get variants from VCF file related to pharma genes and assign them info with rs and genotype
    fg_annotated_variants = annotate_fg_variants(categories_path, norm_vcf, assembly, temp_path)
    
    # Create dictionary with diplotype-phenotype associations
    diplo_pheno_info = get_diplotype_phenotype_dictionary(diplotype_phenotype_info_file)
    
    # Asign diplotype to each pharma gene of interest according to SEFF
    aggregated_fg_results = []
    aggregated_fg_results = assign_cyp2c9_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    aggregated_fg_results = assign_cyp2c19_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    aggregated_fg_results = assign_dpyd_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    aggregated_fg_results = assign_nudt15_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)
    aggregated_fg_results = assign_tpmt_diplotype(fg_annotated_variants, diplo_pheno_info, aggregated_fg_results)

    # Write Pharma results to TSV files
    output_file_variants = f"{norm_vcf.split('norm.FG.vcf.gz')[0]}FG.tsv"
    output_file_diplopheno = f"{norm_vcf.split('norm.FG.vcf.gz')[0]}FG.DiploPheno.tsv"
    write_fg_results_to_tsv(fg_annotated_variants, aggregated_fg_results, output_file_variants, output_file_diplopheno)

    return(fg_annotated_variants, aggregated_fg_results)
    

