# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 2024

@author: jpflorido
"""
from modules.misc.clinvar_utils import run_clinvar
from modules.misc.utils import write_category_results_to_tsv, combine_genebe_clinvar_results
from modules.misc.geneBe_utils import run_genebe, parse_genebe_output
import json
from pathlib import Path


def run_pers_repro_risk_module(mode, category, category_geneset_file, genebe_results_file, clinvar_results_file):
    """
    Run Personal Risk or Reproductive Risgk module

    Args:
        mode (str): Execution mode ("basic" or "advanced").
        category (str): Gene category for annotation
        category_geneset_file (str): Path to CSV file for the given category
        genebe_results_file (str): Path to Genebe annotated VCF file
        clinvar_results_file (str): Path to Clinvar variants for current category
    """

    print("Running " + category.upper() + "risk module")

    # Parse genebe
    genebe_results = parse_genebe_output(genebe_results_file, mode, category, category_geneset_file)
    if mode == "basic":
        category_results = genebe_results
    elif mode == "advanced":
        # Advanced mode: run Clinvar and combine results with Intervar
        with clinvar_results_file.open("r", encoding="utf-8") as fh:
            clinvar_results = json.load(fh)
        genebe_clinvar_results = combine_genebe_clinvar_results(genebe_results, clinvar_results)
        category_results = genebe_clinvar_results

    # Write results of this category to a file
    output_file = clinvar_results_file.parent / f"{category}.SF.csv"
    write_category_results_to_tsv(category_results, str(output_file))
    return category_results
