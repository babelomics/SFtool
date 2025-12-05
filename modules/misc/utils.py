

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
        fieldnames = ["Variant", "Gene", "Genotype", "Consequence", "rs", "Transcript", "HGVSC", "HGVSP", "GeneBe ACMG Classification", "GeneBe ACMG criteria", "Clinvar Clinical Significance", "ReviewStatus", "ClinvarSummary", "ClinvarID", "Orpha"]
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
                "ClinvarSummary": info.get("ClinvarSummary", "-"),
                "ClinvarID": info.get("ClinvarID", "-"),
                "Orpha": info.get("Orpha", "")
            }

            writer.writerow(row)


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