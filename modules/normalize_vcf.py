# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 22:15:15 2023

@author: kindi
"""
import os
import gzip
import subprocess

def normalize_vcf(input_vcf_path, temp_path, bcftools_path, reference_genome_path):
    """
    VCF normalization
    
    Args:
        input_vcf_path (str): La ruta al archivo VCF de entrada que se va a normalizar.
        temp_path (str): La ruta al directorio temporal donde se guardarán los archivos intermedios.
        bcftools_path (str): path a bcftools
        reference_genome_37_path (str): path to reference genome (hs37d5)
    
    Returns:
        str: La ruta del archivo VCF normalizado. Este archivo se encuentra en el directorio temporal.
    """
    # split multiallelic (-m -) y left-alignment.
    try:

        filename, file_extension = os.path.splitext(input_vcf_path)
        just_filename = os.path.basename(input_vcf_path)


        compressed = False
        if file_extension == ".gz":
            compressed = True
            if not (os.path.exists(input_vcf_path + ".csi") or os.path.exists(input_vcf_path + ".tbi")):
                # Index VCF file if not indexed
                index_command = [bcftools_path + "bcftools", "index", input_vcf_path]
                subprocess.run(index_command, check=True)

        # Output files
        if compressed:
            output_vcf_path = os.path.join(temp_path, just_filename.split(".vcf.gz")[0] + ".tmp.vcf.gz")
            output2_vcf_path = os.path.join(temp_path, just_filename.split(".vcf.gz")[0] + ".norm.vcf.gz")
        else:
            output_vcf_path = os.path.join(temp_path, just_filename.split(".vcf")[0] + ".tmp.vcf")
            output2_vcf_path = os.path.join(temp_path, just_filename.split(".vcf")[0] + ".norm.vcf")

        
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

