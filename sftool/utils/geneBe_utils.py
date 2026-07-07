# -*- coding: utf-8 -*-
"""

@author: jpflorido
"""
import subprocess
import os
import gzip
import io
import vcfpy
from sftool.utils.catalog_utils import read_csv
from pathlib import Path
import re

def run_genebe(norm_vcf, category, assembly, genebe_path, java_path, api_key, username, tmp_dir):
    """
    Run GeneBe for annotate variants

    :param norm_vcf: Path to normalized file
    :param category: Gene category for annotation
    :param assembly: Reference genome version
    :param genebe_path: Path to geneBe annotator
    :param java_path: Path to Java
    :param api_key: Api key for annotating using GeneBe
    :param username: User name for annotating using GeneBe
    :param tmp_dir: temporary dir where output file will be saved
    :return: Annotated VCF file
    """

    try:
        # Path to VCF intersected and output directory

        norm_vcf = Path(norm_vcf)

        category_tmp_dir = Path(tmp_dir) / category.upper()
        category_tmp_dir.mkdir(parents=True, exist_ok=True)

        basename = norm_vcf.name.replace(
            f".{category.upper()}.vcf.gz",
            f".{category.upper()}.geneBe.vcf.gz"
        )

        genebe_output_file = category_tmp_dir / basename


        if assembly == 'GRCh37':
            assembly_int = "hg19"
        elif assembly == 'GRCh38':
            assembly_int = 'hg38'

        # Build command to run GeneBe
        cmd = [java_path,
               "-jar",
               genebe_path,
               "vcf",  "annotate",
               "--input-vcf", norm_vcf,
               "--output-vcf", genebe_output_file,
               "--genome", assembly_int,
               "--api-key", api_key,
               "--username", username
               ]

        # Run command and get output
        with subprocess.Popen(cmd, stderr=subprocess.STDOUT, text=True, cwd=os.path.dirname(genebe_path)) as process:
            output, _ = process.communicate()


        return genebe_output_file

    except subprocess.CalledProcessError as e:
        print(f"Error when running Genebe: {e.output}")

def parse_genebe_output(genebe_output_vcf_file, variant_classification_sources, category, category_geneset_file):
    """

    :param genebe_output_vcf_file: VCF annotated by GeneBe
    :param variant_classification_sources: list of variant classification sources
    :param category: pr or rr
    :param category_geneset_file: Path to CSV file for the given category
    :return:
    """

    try:

        # Get the list of genes for the current category
        genes_dct, genes_lst = read_csv(category_geneset_file, category)

        # Read VCF file
        genebe_results = {}

        with gzip.open(genebe_output_vcf_file, "rb") as f:
            text_stream = io.TextIOWrapper(f, encoding="utf-8", errors="replace")  # Convert to text
            vcf_reader = vcfpy.Reader(text_stream)
            for variant_record in vcf_reader:
                chrom = str(variant_record.CHROM)
                pos = str(variant_record.POS)
                ref = str(variant_record.REF)
                alt = str(variant_record.ALT[0].value)
                variant = chrom + ':' + pos + ':' + ref + ':' + alt

                if 'gene_symbol_base' in variant_record.INFO: # There are entries in the VCF file whose ALT is * (avoid those entries which have no gene annotation)

                    # A variant might overlap with more than a single gene. If so, get the information for the gene of interest (contained in the category list)
                    genes_info = [item.split("|") for item in variant_record.INFO['acmg_by_gene_base']]

                    for i, gene_values in enumerate(genes_info, 1):
                        if gene_values[0] in genes_lst:
                            ref_gene = gene_values[0]
                            transcript = gene_values[2]
                            variant_consequence = gene_values[3]
                            acmg_criteria = gene_values[7]
                            classification = gene_values[8]
                            hgvsc = gene_values[9]
                            hgvsp = gene_values[10]

                            genotype = variant_record.calls[0].data['GT'] # A single sample in the VCF is assumed
                            rs = variant_record.INFO.get('dbsnp_base','.')

                            # Get only pathogenic and likely pathogenic variants or add them all if clinvar in variant_classification_sources
                            if classification in ["Pathogenic", "Likely_pathogenic"] or 'clinvar' in variant_classification_sources:
                                # Create a dictionary with interesting fields

                                if variant not in genebe_results:
                                    genebe_results[variant] = []

                                # There might be variants that are annotated to more than one gene (p.e. CYP21A2 – TNXB in RR category)
                                genebe_results[variant].append({
                                    "Gene": ref_gene,
                                    "rs": rs,
                                    "GeneBeClassification": classification,
                                    "Genotype": genotype,
                                    "Consequence": variant_consequence,
                                    "Transcript": transcript,
                                    "ACMG_criteria": acmg_criteria,
                                    "HGVSC": hgvsc,
                                    "HGVSP": hgvsp,
                                    "VCFSampleFormat": "; ".join(
                                        f"{key}: {', '.join(map(str, value)) if isinstance(value, list) else value}"
                                        for key, value in variant_record.calls[0].data.items()
                                    )
                                })

        return genebe_results

    except Exception as e:
        raise Exception(f"Error when parsing Genebe annotated VCF file: {e}")

def get_clinvar_main_gene(variant_name):
    """
    Extracts the gene in the first parentheses of Clinvar's VariantName,
    immediately after the transcript accession.
    """
    match = re.search(r'^[^(]+\(([^)]+)\)', variant_name)
    if match:
        return match.group(1)
    return None



def combine_genebe_clinvar_results(genebe_results, clinvar_results):

    combined_results = {}

    for variant_key in genebe_results.keys():
        for genebe_info in genebe_results.get(variant_key):  # iterate over each entry of genebe

            if variant_key.startswith('chr'):
                clinvar_info = clinvar_results.get(variant_key.removeprefix('chr'))  # Clinvar database does not have chr prefix from CHROM
            else:
                clinvar_info = clinvar_results.get(variant_key)

            if clinvar_info is not None and (get_clinvar_main_gene(clinvar_info["VariantName"]) == genebe_info["Gene"]):  # The gene in both annotation datasets must be the same. Important to check for variants that overlap genes
                clinvar_clinical_significance_tmp = list(map(str.strip,re.split(';|,|/',clinvar_info["ClinicalSignificance"])))
                clinvar_clinical_significance = list(map(str.lower,clinvar_clinical_significance_tmp))

                if (genebe_info and genebe_info["GeneBeClassification"] in ["Pathogenic", "Likely_pathogenic"]) or \
                        ("pathogenic" in clinvar_clinical_significance) or \
                        ("likely pathogenic" in clinvar_clinical_significance) or \
                        (("conflicting classifications of pathogenicity" in clinvar_clinical_significance) and (clinvar_info["ClinSigSimple"]=="1")):

                    if variant_key not in combined_results:
                        combined_results[variant_key] = []

                    combined_results[variant_key].append({
                        "Gene": genebe_info["Gene"],
                        "Genotype": genebe_info["Genotype"],
                        "rs": genebe_info["rs"] if genebe_info["rs"] != '.' else clinvar_info["rs"],
                        "Transcript": genebe_info["Transcript"],
                        "HGVSC": genebe_info["HGVSC"],
                        "HGVSP": genebe_info["HGVSP"],
                        "GeneBeClassification": genebe_info["GeneBeClassification"],
                        "ACMG_criteria": genebe_info["ACMG_criteria"],
                        "ClinvarClinicalSignificance": clinvar_info["ClinicalSignificance"],
                        "ClinvarSummary": clinvar_info["ClinSigSummary"],
                        "ReviewStatus": clinvar_info["ReviewStatus"],
                        "ClinvarID": clinvar_info["ClinvarID"],
                        "Orpha": ",".join(re.findall(r'Orphanet:(\d+)', clinvar_info["PhenotypeIDS"])),
                        "OMIM": ",".join(re.findall(r'OMIM:\s*([^,|]+)', clinvar_info["PhenotypeIDS"])),
                        "Consequence": genebe_info["Consequence"],
                        "VCFSampleFormat": genebe_info["VCFSampleFormat"]
                    })
            else:
                # If there is no info in Clinvar, get info from GeneBe
                if (genebe_info and genebe_info["GeneBeClassification"] in ["Pathogenic", "Likely_pathogenic"]):
                    if variant_key not in combined_results:
                        combined_results[variant_key] = []

                    combined_results[variant_key].append({
                        "Gene": genebe_info["Gene"],
                        "Genotype": genebe_info["Genotype"],
                        "rs": genebe_info["rs"] if genebe_info["rs"] != '.' else '-',
                        "Transcript": genebe_info["Transcript"],
                        "HGVSC": genebe_info["HGVSC"],
                        "HGVSP": genebe_info["HGVSP"],
                        "GeneBeClassification": genebe_info["GeneBeClassification"],
                        "ACMG_criteria": genebe_info["ACMG_criteria"],
                        "ClinvarClinicalSignificance": "NA",
                        "ClinvarSummary": "NA",
                        "ReviewStatus": "NA",
                        "ClinvarID": "NA",
                        "Orpha": "NA",
                        "OMIM": "NA",
                        "Consequence": genebe_info["Consequence"],
                        "VCFSampleFormat": genebe_info["VCFSampleFormat"]
                    })

    return combined_results