# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 2024

@author: jpflorido
"""
from modules.misc.intervar_utils import run_intervar, parse_intervar_output
from modules.misc.clinvar_utils import run_clinvar
from modules.misc.vcf_utils import combine_results
from modules.misc.utils import write_category_results_to_tsv


def run_pers_repro_risk_module(norm_vcf, assembly, mode, evidence_level, clinvar_db, intervar_path, category, category_geneset_file):
    """
    Run Personal Risk or Reproductive Risgk module

    Args:
        vcf_path (str): Path to normalized and intersected VCF file
        assembly (str): Reference genome version
        mode (str): Execution mode ("basic" or "advanced").
        evidence_level (int): Evidence level
        category (str): Gene category for annotation
        clinvar_db (str): Path to CLINVAR database
        category_geneset_file (str): Path to CSV file for the given category
    """

    print("Running " + category.upper() + "risk module")

    # Intervar is always run
    intervar_output_file = run_intervar(norm_vcf, category, assembly, intervar_path)
    intervar_results = parse_intervar_output(intervar_output_file, mode, assembly)
    if mode == "basic":
        category_results = intervar_results
    elif mode == "advanced":
        # Advanced mode: run Clinvar and combine results with Intervar
        clinvar_results = run_clinvar(evidence_level, clinvar_db, category, category_geneset_file)
        intervar_clinvar_results = combine_results(norm_vcf, intervar_results, clinvar_results)
        category_results = intervar_clinvar_results

    # Write results of this category to a file
    output_file = f"{norm_vcf.split('norm.' + category.upper() + '.vcf.gz')[0]}{category.upper()}.SF.tsv"
    write_category_results_to_tsv(category_results, output_file)
    return category_results
