

import csv

def write_category_results_to_tsv(results, norm_vcf, category):
    """
    Escribe los resultados combinados en un archivo TSV.

    Args:
        results (list): Lista de diccionarios con los resultados combinados.
        norm_vcf (str): Ruta al archivo normalizado.
        category (str): Categoría de genes para la anotación.
    """
    output_tsv = f"{norm_vcf.split('normalized')[0]}{category}_all_results.tsv"

    # Abrir el archivo TSV para escritura
    with open(output_tsv, "w", newline="") as tsv_file:
        fieldnames = ["Variant", "Gene", "Genotype", "Intervar Consequence", "rs", "Intervar Classification","Clinvar Clinical Significance", "ReviewStatus", "ClinvarID", "Orpha"]
        writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")

        # Escribir el encabezado del archivo TSV
        writer.writeheader()

        # Iterar sobre las claves (variantes) en el diccionario result
        for variant, info in results.items():
            # Asegurarse de que el diccionario 'info' tenga todas las claves necesarias
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

            # Escribir la fila en el archivo TSV
            writer.writerow(row)



