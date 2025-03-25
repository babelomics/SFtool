

import csv
import re

def write_category_results_to_tsv(results, output_tsv):
    """
    Write combined results from Intervar and Clinvar to a TSV file

    Args:
        results (list): List of dictionaries with combined results
        output_tsv (str): Output TSV file
    """

    with open(output_tsv, "w", newline="") as tsv_file:
        fieldnames = ["Variant", "Gene", "Genotype", "Consequence", "rs", "Transcript", "HGVSC", "HGVSP", "GeneBe ACMG Classification", "GeneBe ACMG criteria", "Clinvar Clinical Significance", "ReviewStatus", "ClinvarID", "Orpha"]
        writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")

        writer.writeheader()

        for variant, info in results.items():
            # Be sure that info dictionary has all necessary keys
            row = {
                "Variant": variant,
                "Gene": info.get("Gene", ""),
                "Genotype": info.get("Genotype", ""),
                "Consequence": info.get("Consequence", "-"),
                "rs": info.get("rs", ""),
                "Transcript": info.get("Transcript",""),
                "HGVSC": info.get("HGVSC",""),
                "HGVSP": info.get("HGVSP",""),
                "GeneBe ACMG Classification": info.get("GeneBeClassification", ""),
                "GeneBe ACMG criteria": info.get("ACMG_criteria",""),
                "Clinvar Clinical Significance": info.get("ClinvarClinicalSignificance", "-"),
                "ReviewStatus": info.get("ReviewStatus", "-"),
                "ClinvarID": info.get("ClinvarID", "-"),
                "Orpha": info.get("Orpha", "")
            }

            writer.writerow(row)



def combine_genebe_clinvar_results(genebe_results, clinvar_results):

    combined_results = {}

    for variant_key in genebe_results.keys():
        genebe_info= genebe_results.get(variant_key)

        if variant_key.startswith('chr'):
            clinvar_info = clinvar_results.get(variant_key.removeprefix('chr')) # Clinvar database does not have chr prefix from CHROM
        else:
            clinvar_info = clinvar_results.get(variant_key)

        if clinvar_info is not None:
            clinvar_clinical_significance_tmp = list(map(str.strip,re.split(';|,|/',clinvar_info["ClinicalSignificance"])))
            clinvar_clinical_significance = list(map(str.lower,clinvar_clinical_significance_tmp))

            if (genebe_info and genebe_info["GeneBeClassification"] in ["Pathogenic", "Likely_pathogenic"]) or \
                    ("pathogenic" in clinvar_clinical_significance) or \
                    ("likely pathogenic" in clinvar_clinical_significance) or \
                    (("conflicting classifications of pathogenicity" in clinvar_clinical_significance) and (clinvar_info["ClinSigSimple"]=="1")):

                combined_results[variant_key] = {
                    "Gene": genebe_info["Gene"],
                    "Genotype": genebe_info["Genotype"],
                    "rs": genebe_info["rs"] if genebe_info["rs"] != '.' else clinvar_info["rs"],
                    "Transcript": genebe_info["Transcript"],
                    "HGVSC": genebe_info["HGVSC"],
                    "HGVSP": genebe_info["HGVSP"],
                    "GeneBeClassification": genebe_info["GeneBeClassification"],
                    "ACMG_criteria": genebe_info["ACMG_criteria"],
                    "ClinvarClinicalSignificance": clinvar_info["ClinicalSignificance"],
                    "ReviewStatus": clinvar_info["ReviewStatus"],
                    "ClinvarID": clinvar_info["ClinvarID"],
                    "Orpha": ",".join(re.findall(r'Orphanet:(\d+)', clinvar_info["PhenotypeIDS"])),
                    "Consequence": genebe_info["Consequence"]
                }
        else:
            # If there is no info in Clinvar, get info from GeneBe
            if (genebe_info and genebe_info["GeneBeClassification"] in ["Pathogenic", "Likely_pathogenic"]):
                combined_results[variant_key] = {
                    "Gene": genebe_info["Gene"],
                    "Genotype": genebe_info["Genotype"],
                    "rs": genebe_info["rs"] if genebe_info["rs"] != '.' else '-',
                    "Transcript": genebe_info["Transcript"],
                    "HGVSC": genebe_info["HGVSC"],
                    "HGVSP": genebe_info["HGVSP"],
                    "GeneBeClassification": genebe_info["GeneBeClassification"],
                    "ACMG_criteria": genebe_info["ACMG_criteria"],
                    "ClinvarClinicalSignificance": "NA",
                    "ReviewStatus": "NA",
                    "ClinvarID": "NA",
                    "Orpha": "NA",
                    "Consequence": genebe_info["Consequence"]
                }

    return combined_results