# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 23:26:14 2023

@author: Edurne Urrutia, Javier Pérez Florido
"""

import json
import pandas as pd
import os
import subprocess


def check_specific_criteria(gene, variant, variant_key, assembly):
    """
    Check specific features for given gene: either specific consequence or specific variant

    :param gene: gene information from current catalogue
    :param variant: variant information under study
    :param variant_key: specific genomic variant under study
    :param assembly: genome assembly version, either 37 or 38
    :return: TRUE if criteria is met, else, FALSE
    """

    meet_criteria = True

    if gene["specific_consequence"] != "":
        variant_consequences = variant["IntervarConsequence"].split(",")
        gene_consequences = gene["specific_consequence"].split(",")

        # If required consequence from gene information is not present in the variant consequence, criteria is not meet (variant is discarded)
        if len(set(variant_consequences) & set(gene_consequences)) == 0:
            meet_criteria = False
    elif gene["specific_variant_GRCh" + str(assembly)] != "":
        specific_variant = gene["specific_variant_GRCh" + str(assembly)].split(',')[0]
        specific_genotype = gene["specific_variant_GRCh" + str(assembly)].split(',')[1]

        if specific_variant != variant_key and specific_genotype.lower() != variant["Genotype"]:
            meet_criteria = False

    return meet_criteria


def get_versions_paths(program_arguments, config_data, clinvar_db):
    """
    Get versions of third-party tools from SF tool and program arguments to be shown in the final report
    :param program_arguments: SF tools arguments
    :param config_data: Configuration data
    :param clinvar_db: Clinvar database path
    :return:
    """
    # HPO
    if program_arguments.hpos_file is None:
        hpos_patient_list = "Not provided"
    else:
        hpos_list_temp = get_hpos_from_txt(program_arguments.hpos_file)
        hpos_patient_list=",".join(str(hpo) for hpo in hpos_list_temp)

    # Intervar version
    try:
        intervar_file_path = os.path.join(config_data["intervar_path"], "Intervar.py")
        cmd = [intervar_file_path, "--version"]

        intervar_process = subprocess.Popen(cmd, stdout= subprocess.PIPE)
        intervar_out, intervar_err = intervar_process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error running InterVar: {e.output}")

    # Bcftools version
    try:
        bcftools_file_path = os.path.join(config_data["bcftools_path"], "bcftools")
        cmd = [bcftools_file_path, "--version"]

        bcftools_process = subprocess.Popen(cmd, stdout= subprocess.PIPE)
        bcftools_out, bcftools_err = bcftools_process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error running InterVar: {e.output}")



    versions_paths={
        "SF module": "0.1",
        "SF module mode": program_arguments.mode,
        "Input VCF": program_arguments.vcf_file,
        "Personal risk catalogue file": config_data["personal_risk_geneset_file"],
        "Reproductive risk catalogue file": config_data["reproductive_risk_geneset_file"],
        "Pharmacogenetic risk catalogue file": config_data["pharmacogenetic_risk_variant_GRCh37_file"],
        "Diplotype-Phenotype file": config_data["diplotype_phenotype_info_file"],
        "Temporal dir": config_data["temp_path"],
        "Output dir": config_data["out_path"],
        "HPO list": hpos_patient_list,
        "Human assembly": "hg19" if program_arguments.assembly == 37 else "hg38",
        "Reference genome path":  config_data["reference_genome_37_path"] if program_arguments.assembly == 37 else config_data["reference_genome_38_path"],
        "Clinvar version": "Not used (basic mode)" if program_arguments.mode == "basic" else os.path.splitext(os.path.basename(clinvar_db))[0].split("_")[-1],
        "Clinvar path": clinvar_db,
        "Clinvar Evidence level": str(program_arguments.evidence),
        "Intervar version": str(intervar_out).split(" ")[1] + " " + str(intervar_out).split(" ")[2].split("\\n")[0],
        "Intervar path": config_data["intervar_path"],
        "bcftools version": str(bcftools_out).split(" ")[1].split("\\n")[0],
        "bcftools path": config_data["bcftools_path"],
        "HPO genes to phenotype version": os.path.splitext(os.path.basename(config_data["gene_to_phenotype_file"]))[0].split("_")[-1],
        "HPO genes to phenotype path": config_data["gene_to_phenotype_file"]

    }

    return versions_paths

def combine_variant_and_gene_info(variant_info, gene_info):
    """
    Combine information of variant and gene

    Args:
        variant_info (dict): Variant info
        gene_info (dict): Gene info

    Returns:
        dict: Merged information
    """
    combined_info = {
        "Gene": variant_info["Gene"],
        "Genotype": variant_info["Genotype"],
        "rs": variant_info.get("rs", ""),
        "IntervarConsequence": variant_info.get("IntervarConsequence", ""),
        "IntervarClassification": variant_info["IntervarClassification"],
        "ClinvarClinicalSignificance": variant_info.get("ClinvarClinicalSignificance", "-"),
        "ReviewStatus": variant_info.get("ReviewStatus", "-"),
        "ClinvarID": variant_info.get("ClinvarID", "-"),
        "Orpha": variant_info.get("Orpha", ""),
        "Phenotype": gene_info["phenotype"],
        "ACMG_version": gene_info.get("ACMG_version", ""),  # Usar get para manejar la falta de 'ACMG_version'
        "OMIM_disorder": gene_info["OMIM_disorder"],
        "inheritance": gene_info["inheritance"],
        "variants_to_report": gene_info.get("variants_to_report", ""),  # Usar get para manejar la falta de 'variants_to_report'
        "related_HPOs_for_sample": 'NA'
    }
    return combined_info

def check_inheritance(results, category, categories_path, assembly):
    """
    Check inheritance for the set of variants according to the category (personal risk or reproductive risk) and generate a dictionary with variants to be
    informed according to inheritance specified for genes in the corresponding catalogue

    Args:
        results (dict): A dictionary with results from PR or RR module
        category (str): category: either PR or RR
        categories_path (str): Path to categories directory
        assembly (str): genome reference assembly. Either 37 or 38

    Returns:
        dict: A dictionary with the variants to be informed according to the inheritance defined for the set of genes in the catalogue
    """

    try:
        # Load JSON file for the given category. This file contains the inheritance mode for each gene
        genes_cat_path = f"{categories_path}{category.upper()}/{category}_risk_genes.json"
        genes_cat = None
        with open(genes_cat_path, "r") as genes_cat_file:
            genes_cat = json.load(genes_cat_file)

        # Create a dictionary with the set of variants to be informed
        reported_variants = {}

        # Move through the set of combined results
        for variant_key, variant_info in results.items():
            variant_gene = variant_info["Gene"]
            # Get inheritance mode for the given gene
            for gene in genes_cat['genes']:
                if gene['gene_symbol'] ==  variant_gene:
                    inher = gene["inheritance"]

                    # Check for specific consequence or genomic variant
                    if check_specific_criteria(gene, variant_info, variant_key, assembly):
                        # For Personal Risk module, report variant if autosomic Dominant (AD), SemiDominant (SD) and X-linked (XL) inheritance mode. For Reproductive Risk category, report variant in any case
                        if inher in ['AD', 'SD', 'XL'] or category == 'rr':
                            # Combina la información de la variante y el gen
                            combined_info = combine_variant_and_gene_info(variant_info, gene)
                            # Agrega la información combinada a reported_variants
                            reported_variants[variant_key] = combined_info

                        # For Autosomic Recessive, check genotype and/or other variants in the same gene
                        elif inher == 'AR':
                            # For a variant in HOM, report variant
                            if variant_info["Genotype"] == 'hom':
                                # Combine information for gene and variant
                                combined_info = combine_variant_and_gene_info(variant_info, gene)
                                # Add merged information to the dictionary
                                reported_variants[variant_key] = combined_info

                            # For a variant in HET, only report if there is another variant in the same gene
                            elif variant_info["Genotype"] == 'het':
                                # Look for another variant in the same gene
                                other_variant_in_gene = False
                                for other_variant_key in results:
                                    other_variant_info = results[other_variant_key]
                                    if other_variant_info["Gene"] == variant_gene and other_variant_key != variant_key:
                                        other_combined_info = combine_variant_and_gene_info(other_variant_info, gene)
                                        combined_info = combine_variant_and_gene_info(variant_info, gene)

                                        # Add merged information of the two variants of the gene to the dictionary
                                        reported_variants[variant_key] = combined_info
                                        reported_variants[other_variant_key] = other_combined_info
        return reported_variants
    except Exception as e:
        print(f"Error in check_inheritance function: {str(e)}")
        return {}

def check_patient_HPO(reported_results, hpos_user, gene_to_phenotype_file):
    """
    Check if the set of HPOs provided for the sample must be reported according the list of HPOs described for genes

    Args:
        reported_results (dict): Results from PR or RR module
        user_hpos (list): HPO list provided for the sample
        gene_to_phenotype_file (str): Gene to phenotype file containing correspondence between HPO, Genes and OMIM terms
                                      (https://hpo.jax.org/app/data/annotations, genes to phenotype, https://github.com/obophenotype/human-phenotype-ontology/releases/)

    Returns:
        list: List of matching results.
    """
    gene_hpo_dict = {}

    # Process file with information of HPO, Genes and OMIM terms
    with open(gene_to_phenotype_file, 'r') as file:
        next(file)
        for line in file:
            fields = line.strip().split('\t')
            gene_symbol = fields[1]
            hpo_id = fields[2]

            # If gene symbol is in the dictionary, add the new HPO term
            if gene_symbol in gene_hpo_dict:
                gene_hpo_dict[gene_symbol].append(hpo_id)
            else:  # If gene symbol is not in the dictionary, create a new entry
                gene_hpo_dict[gene_symbol] = [hpo_id]

    for result in reported_results:
        # Check wether HPOs for a given gene according to gene_to_phenotype_file are in the list of HPOs for a sample
        gene = reported_results[result]['Gene']
        hpo_results = gene_hpo_dict[gene]

        for hpo in hpo_results:
            if hpo in hpos_user:
                if reported_results[result]['related_HPOs_for_sample'] == 'NA':
                    reported_results[result]['related_HPOs_for_sample'] = hpo
                else:
                    if hpo not in reported_results[result]['related_HPOs_for_sample'].split(','):
                        reported_results[result]['related_HPOs_for_sample'] += ',' + hpo

    return reported_results

def get_hpos_from_txt(hpos_file):
    """
    Read a txt file that contains a list of HPOs for a given patient (one HPO per line) and returns a list of such HPOs

    Args:
        hpos_file (str): Path to file that contains the list of HPOs

    Returns:
        List[str]: HPO list

    Raises:
        FileNotFoundError: if file is not found
    """
    if hpos_file is None:
        return []  # If there is no file, returns an empty list

    try:
        with open(hpos_file, "r") as file:
            hpos = [line.strip() for line in file.readlines()]
        return hpos

    except FileNotFoundError:
        print(f"File {hpos_file} not found.")
        return []

def generate_report(pr_results, rr_results, fg_results, haplot_results, config_data, args, clinvar_db, categories):
    """
    Write results for all categories to an excel file

    Args:
        pr_results (dict): Results from Personal Risk module
        rr_results (dict): Results from Reproductive Risk module
        fg_results (dict): Results from Pharmacogenetic category module
        out_path (str): Path to output directory
    """
    try:

        categories_path = config_data["categories_path"]
        out_path = config_data["out_path"]
        vcf_file = args.vcf_file
        hpos_file = args.hpos_file
        gene_to_phenotype_file = config_data["gene_to_phenotype_file"]
        assembly = args.assembly

        # Get versions
        versions_path = get_versions_paths(args, config_data, clinvar_db)

        # Get HPO list
        hpos_user= []
        if versions_path['HPO list'] != "Not provided":
            hpos_user = versions_path['HPO list'].split(',')
        outfile = f"{out_path}{vcf_file.split('/')[-1].split('.vcf')[0]}.SF.xlsx"

        with pd.ExcelWriter(outfile) as writer:
            # Write versions and paths
            summary_df = pd.DataFrame.from_dict(versions_path, orient="index")
            summary_df.to_excel(writer, sheet_name= 'Versions and Paths', index=True, header=False)

            for category in categories:
                if category == 'pr':
                    reported_results = check_inheritance(pr_results, category, categories_path, assembly)
                    pr_final = check_patient_HPO(reported_results, hpos_user, gene_to_phenotype_file)
                    results_df = pd.DataFrame.from_dict(pr_final, orient='index')
                    results_df.to_excel(writer, sheet_name= category.upper() + ' results', index=True)
                elif category == 'rr':
                    reported_results = check_inheritance(rr_results, category, categories_path, assembly)
                    rr_final = check_patient_HPO(reported_results, hpos_user, gene_to_phenotype_file)
                    results_df = pd.DataFrame.from_dict(rr_final, orient='index')
                    results_df.to_excel(writer, sheet_name= category.upper() + ' results', index=True)
                else:
                    fg_df = pd.DataFrame(fg_results)
                    fg_df.to_excel(writer, sheet_name='FG results', index=False)

                    haplot_df = pd.DataFrame(haplot_results)
                    haplot_df.to_excel(writer, sheet_name='FG Diplotype-Phenotype', index=False)

        print(f"Results have been written to '{outfile}'.")

    except Exception as e:
        print(f"Error found when generating final report: {str(e)}")
