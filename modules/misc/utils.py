

import csv

def write_category_results_to_tsv(results, output_tsv):
    """
    Write combined results from Intervar and Clinvar to a TSV file

    Args:
        results (list): List of dictionaries with combined results
        output_tsv (str): Output TSV file
    """

    with open(output_tsv, "w", newline="") as tsv_file:
        fieldnames = ["Variant", "Gene", "Genotype", "Intervar Consequence", "rs", "Intervar Classification","Clinvar Clinical Significance", "ReviewStatus", "ClinvarID", "Orpha"]
        writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")

        writer.writeheader()

        for variant, info in results.items():
            # Be sure that info dictionary has all necessary keys
            row = {
                "Variant": variant,
                "Gene": info.get("Gene", ""),
                "Genotype": info.get("Genotype", ""),
                "Intervar Consequence": info.get("IntervarConsequence", "-"),
                "rs": info.get("rs", ""),
                "Intervar Classification": info.get("IntervarClassification", ""),
                "Clinvar Clinical Significance": info.get("ClinvarClinicalSignificance", "-"),
                "ReviewStatus": info.get("ReviewStatus", "-"),
                "ClinvarID": info.get("ClinvarID", "-"),
                "Orpha": info.get("Orpha", "")
            }

            writer.writerow(row)



