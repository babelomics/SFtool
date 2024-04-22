# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 22:15:15 2023

@author: Edurne Urrutia, Javier Perez Florido
"""
import os
import subprocess
from pybedtools import BedTool


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

        just_filename = os.path.basename(input_vcf_path)

        if not (os.path.exists(input_vcf_path + ".csi") or os.path.exists(input_vcf_path + ".tbi")):
            # Index VCF file if not indexed
            index_command = [bcftools_path + "bcftools", "index", input_vcf_path]
            subprocess.run(index_command, check=True)

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

