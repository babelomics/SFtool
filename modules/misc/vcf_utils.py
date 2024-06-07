# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 22:15:15 2023

@author: Edurne Urrutia, Javier Perez Florido
"""
import os
import subprocess
from pybedtools import BedTool
import vcfpy
import re

def normalize_vcf(input_vcf_path, temp_path, bcftools_path, reference_genome_path):
    """
    VCF normalization

    Args:
        input_vcf_path (str): La ruta al archivo VCF de entrada que se va a normalizar.
        temp_path (str): La ruta al directorio temporal donde se guardarán los archivos intermedios.
        bcftools_path (str): path a bcftools
        reference_genome_path (str): path to reference genome

    Returns:
        str: La ruta del archivo VCF normalizado. Este archivo se encuentra en el directorio temporal.
    """
    # split multiallelic (-m -) y left-alignment.
    try:
        print("Normalizing " + input_vcf_path + " file...")

        just_filename = os.path.basename(input_vcf_path)

        if not (os.path.exists(input_vcf_path + ".csi") or os.path.exists(input_vcf_path + ".tbi")):
            # Index VCF file if not indexed
            index_command = [bcftools_path + "bcftools", "index", input_vcf_path]
            subprocess.run(index_command, check=True, capture_output=True)

        # Output files
        output_vcf_path = os.path.join(temp_path, just_filename.split(".vcf.gz")[0] + ".tmp.vcf.gz")
        output2_vcf_path = os.path.join(temp_path, just_filename.split(".vcf.gz")[0] + ".norm.vcf.gz")

        # bcftools normalization command

        bcftools_command = [bcftools_path + "bcftools", "norm", "-O", "z", "-m", "-any", "--check-ref", "w",  "-f", reference_genome_path, "-o", output_vcf_path, input_vcf_path]

        # Normalize with bcftools
        subprocess.run(bcftools_command, check=True)

        # Remove duplicates with bcftools
        rm_dup_command = [bcftools_path + "bcftools", "norm", "--rm-dup", "none", "-Oz", "-o", output2_vcf_path, output_vcf_path]
        subprocess.run(rm_dup_command, check=True)

        print("bcftools normalization completed.")

        os.remove(output_vcf_path)
        return(output2_vcf_path)

    except Exception as e:
        print(f"Error given by bcftools normalization: {e}")


def intersect_vcf_with_bed(vcf_norm_file, category_bed_file, temp_path, category):
    """
    Intersection of VCF and BED file. Variants are saved to a new VCF file

    Args:
        vcf_norm_file (str): Normalized VCF file.
        category_bed_file (str): BED file with positions of interest for the given category
        temp_path (str): temporal directory where VCF files are generated
        category (str): Category: PR, RR or FG.

    Returns:
        None

    Raises:
        Exception: An error is rised if an error occurrs with the intersection
    """
    try:

        print("Intersecting VCF file with BED file ( " + category.upper() + " category)...")
        just_filename = os.path.basename(vcf_norm_file)
        output_vcf_path = os.path.join(temp_path, just_filename.split(".vcf.gz")[0] + "." + category.upper() + ".vcf.gz")

        # Load VCF and BED files using Python's BedTools
        vcf = BedTool(vcf_norm_file)
        bed = BedTool(category_bed_file)

        # Intersection of BED and VCF files
        intersected_variants = vcf.intersect(bed, u=True, header=True)

        # Save intersected VCF file
        intersected_variants.saveas(output_vcf_path)

        return output_vcf_path

        print(f"Intersection completed. VCF file saved to {output_vcf_path}")
    except Exception as e:
        print(f"Error during the VCF intersection {e}")

def combine_results(vcf_norm, intervar_results, clinvar_results):
    """
    Combine results from Intervar and Clinvar into a single line per variant

    Args:
        vcf_norm (str): Path to normalized and intersected VCF file
        intervar_results (dict): Intervar results
        clinvar_results (dict): Clinvar results

    Returns:
        dict: A diccionary with combined results
    """
    combined_results = {}

    # Since Intervar removes reference nucleotide in indels, the best way to combine results
    # is to parse VCF file and search for a variant in Intervar and Clinvar results


    try:
        # Read VCF file
        vcf_reader = vcfpy.Reader.from_path(vcf_norm)
        for variant_record in vcf_reader:
            chrom = str(variant_record.CHROM)
            pos = str(variant_record.POS)
            ref = str(variant_record.REF)
            alt = str(variant_record.ALT[0].value)
            variant_key = chrom + ':' + pos + ':' + ref + ':' + alt

            if len(ref) > len(alt): # Deletion
                if len(alt) == 1:
                    variant_int = f"{chrom}:{str(int(pos) + 1)}:{ref[1:]}:-"
                else:
                    variant_int = f"{chrom}:{str(int(pos) + 1)}:{ref[1:]}:{alt[1:]}"
            elif len(ref) < len(alt): # Insertion
                if len(ref) == 1:
                    variant_int = f"{chrom}:{pos}:-:{alt[1:]}"
                else:
                    variant_int = f"{chrom}:{pos}:{ref[1:]}:{alt[1:]}"
            else:
                variant_int = f"{chrom}:{pos}:{ref}:{alt}"


            intervar_info = intervar_results.get(variant_int)

            # Search variant in Clinvar dictionary
            clinvar_info = clinvar_results.get(variant_key)


            if clinvar_info is not None:

                clinvar_clinical_significance_tmp = list(map(str.strip,re.split(';|,|/',clinvar_info["ClinicalSignificance"])))
                clinvar_clinical_significance = list(map(str.lower,clinvar_clinical_significance_tmp))

                if (intervar_info and intervar_info["IntervarClassification"] in ["Pathogenic", "Likely pathogenic"]) or \
                        ("pathogenic" in clinvar_clinical_significance) or \
                        ("likely pathogenic" in clinvar_clinical_significance) or \
                        (("conflicting interpretations of pathogenicity" in clinvar_clinical_significance) and (clinvar_info["ClinSigSimple"]=="1")):

                    combined_results[variant_key] = {
                        "Gene": clinvar_info["Gene"],
                        "Genotype": intervar_info["Genotype"],
                        "rs": intervar_info["rs"] if intervar_info["rs"] != '.' else clinvar_info["rs"],
                        "IntervarClassification": intervar_info["IntervarClassification"],
                        "ClinvarClinicalSignificance": clinvar_info["ClinicalSignificance"],
                        "ReviewStatus": clinvar_info["ReviewStatus"],
                        "ClinvarID": clinvar_info["ClinvarID"],
                        "Orpha": intervar_info["Orpha"],
                        "IntervarConsequence": intervar_info["IntervarConsequence"]
                    }
            else:
                # If there is no info in Clinvar, get info from Intervar
                if (intervar_info and intervar_info["IntervarClassification"] in ["Pathogenic", "Likely pathogenic"]):
                    combined_results[variant_key] = {
                        "Gene": intervar_info["Gene"],
                        "Genotype": intervar_info["Genotype"],
                        "rs": intervar_info["rs"] if intervar_info["rs"] != '.' else '-',
                        "IntervarClassification": intervar_info["IntervarClassification"],
                        "ClinvarClinicalSignificance": "NA",
                        "ReviewStatus": "NA",
                        "ClinvarID": "NA",
                        "Orpha": intervar_info["Orpha"],
                        "IntervarConsequence": intervar_info["IntervarConsequence"]
                    }

    except Exception as e:
        raise Exception(f"Error al combinar resultados: {e}")

    return combined_results

