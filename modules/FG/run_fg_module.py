# -*- coding: utf-8 -*-
"""
Created on Wed March  19 17:21:45 2025

@author: jpflorido

"""

import subprocess
import os
import json
import pandas as pd

def run_pharmCAT_vcf_preprocessor(vcf_input, python_path, pharmCAT_path, bcftools_path, htslib_path):
    """
    Run pharmCAT's VCF preprocessor: https://pharmcat.org/using/VCF-Preprocessor/

    :param vcf_input:
    :param python_path:
    :param pharmCAT_path:
    :return:
    """

    preprocessed_vcf = ''

    try:
        # Run pharmCAT's VCF preprocessor script
        vcf_preprocessor_command = [python_path, pharmCAT_path + "/pharmcat_vcf_preprocessor.py", "--path-to-bcftools", bcftools_path + "bcftools", "--path-to-bgzip", htslib_path + "bgzip", "-vcf", vcf_input]

        with subprocess.Popen(vcf_preprocessor_command, stderr=subprocess.STDOUT, text=True, cwd=pharmCAT_path) as process:
            output, _ = process.communicate()

        preprocessed_vcf = vcf_input.split(".vcf.gz")[0] + ".preprocessed.vcf.bgz"
        if os.path.exists(preprocessed_vcf):
            return preprocessed_vcf
        else:
            print("pharmCAT's preprocessed file does not exist. Exiting")
            exit(-1)

    except subprocess.CalledProcessError as e:
        print(f"Error when running pharmCAT's VCF preprocessor script: {e.output}")


def run_pharmCAT(preprocessed_vcf, pharmCAT_path, java_path, out_path):
    """
    Run pharmCAT

    :param preprocessed_vcf:
    :param out_path:
    :return:
    """

    try:
        pharmCAT_command = [java_path, "-jar", pharmCAT_path + "/pharmcat.jar", "-vcf", preprocessed_vcf, "--output-dir", out_path]
        with subprocess.Popen(pharmCAT_command, stderr=subprocess.STDOUT, text=True, cwd=pharmCAT_path) as process:
            output, _ = process.communicate()

        file_name_prefix = os.path.basename(preprocessed_vcf).split(".preprocessed.vcf.bgz")[0]

        report_file = os.path.join(out_path, file_name_prefix + ".report.html")
        phenotype_file = os.path.join(out_path, file_name_prefix + ".phenotype.json")

        if os.path.exists(report_file) and os.path.exists(phenotype_file):
            return [report_file, phenotype_file]
        else:
            print("pharmCAT report could not be generated. Exiting")
            exit(-1)
    except subprocess.CalledProcessError as e:
        print(f"Error when running pharmCAT: {e.output}")


def parse_pharmCAT(phenotype_file):
    """
    Parse pharmCAT JSON file and return a data frame with basic information of genes, diplotypes and phenotypes

    :param phenotype_file:
    :return:
    """

    with open(phenotype_file, "r") as fd:
        data = json.load(fd)

    # Extract the required information

    diplotype_data_cpic = []
    diplotype_data_dpwg = []

    gene_reports = data.get("geneReports", {}).get("CPIC", {})

    for gene, details in gene_reports.items():
        recommendation_diplotypes = details.get("recommendationDiplotypes", [])

        for rec in recommendation_diplotypes:
            phenotypes = ",".join(map(str, rec.get("phenotypes", "N/A")))
            label = rec.get("label", "N/A")

            if label != "Unknown/Unknown" and label != "Unknown":
                diplotype_data_cpic.append({"Gene": gene, "Genotype": label, "Phenotype": phenotypes, "Source": "CPIC"})
            else:
                diplotype_data_cpic.append({"Gene": gene, "Genotype": "Not determined", "Phenotype": "Not determined", "Source": "-"})

    gene_reports = data.get("geneReports", {}).get("DPWG", {})

    for gene, details in gene_reports.items():
        recommendation_diplotypes = details.get("recommendationDiplotypes", [])

        for rec in recommendation_diplotypes:
            phenotypes = ",".join(map(str, rec.get("phenotypes", "N/A")))
            label = rec.get("label", "N/A")

            if label != "Unknown/Unknown" and label != "Unknown":
                diplotype_data_dpwg.append({"Gene": gene, "Genotype": label, "Phenotype": phenotypes, "Source": "DPWG"})


    diplotype_data = pd.concat([pd.DataFrame(diplotype_data_cpic), pd.DataFrame(diplotype_data_dpwg)], ignore_index=True)

    # Sorting to keep identical genes together while maintaining order
    diplotype_data = diplotype_data.sort_values(by=["Gene"], key=lambda x: x.map({gene: i for i, gene in enumerate(diplotype_data["Gene"].unique())}))

    # Reset index for a clean structure
    diplotype_data.reset_index(drop=True, inplace=True)


    return diplotype_data


def run_pharmacogenomic_risk_module(vcf_input, python_path, pharmCAT_path, bcftools_path, htslib_path, java_path, out_path):
    """
    :param vcf_input:
    :return:
    """

    # 1. Run pharmCAT's preprocessor script (https://pharmcat.org/using/VCF-Preprocessor/)
    preprocessed_vcf = run_pharmCAT_vcf_preprocessor(vcf_input, python_path, pharmCAT_path, bcftools_path, htslib_path)

    # 2. Run pharmCAT
    [report_file, phenotype_file] = run_pharmCAT(preprocessed_vcf, pharmCAT_path, java_path, out_path)

    # 3. Besides pharmCAT results, generate a short report with genes, diplotypes and phenotype
    haplot_results = parse_pharmCAT(phenotype_file)

    return [report_file, haplot_results]