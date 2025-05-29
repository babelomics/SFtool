# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 2024

@author: jpflorido
"""
from modules.misc.clinvar_utils import run_clinvar
from modules.misc.utils import write_category_results_to_tsv, combine_genebe_clinvar_results
from modules.misc.geneBe_utils import run_genebe, parse_genebe_output


def run_pers_repro_risk_module(norm_vcf, assembly, mode, evidence_level, clinvar_db, clinvar_submission, category, category_geneset_file, genebe_path, java_path, genebe_apikey, genebe_username):
    """
    Run Personal Risk or Reproductive Risgk module

    Args:
        vcf_path (str): Path to normalized and intersected VCF file
        assembly (str): Reference genome version
        mode (str): Execution mode ("basic" or "advanced").
        evidence_level (int): Evidence level
        category (str): Gene category for annotation
        clinvar_db (str): Path to CLINVAR database
        clinvar_submission (str): Path to CLINVAR submission summary
        category_geneset_file (str): Path to CSV file for the given category
    """

    print("Running " + category.upper() + "risk module")

    # Run GeneBe
    genebe_output_file = run_genebe(norm_vcf, category, assembly, genebe_path, java_path, genebe_apikey, genebe_username)
    genebe_results = parse_genebe_output(genebe_output_file, mode, category, category_geneset_file)
    if mode == "basic":
        category_results = genebe_results
    elif mode == "advanced":
        # Advanced mode: run Clinvar and combine results with Intervar
        clinvar_results = run_clinvar(evidence_level, clinvar_db, clinvar_submission, category, category_geneset_file)
        genebe_clinvar_results = combine_genebe_clinvar_results(genebe_results, clinvar_results)
        category_results = genebe_clinvar_results

    # Write results of this category to a file
    output_file = f"{norm_vcf.split('norm.' + category.upper() + '.vcf.gz')[0]}{category.upper()}.SF.tsv"
    write_category_results_to_tsv(category_results, output_file)
    return category_results
