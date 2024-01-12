# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 22:15:15 2023

@author: kindi
"""
import os
import gzip
import subprocess

def normalize_vcf(input_vcf_path, temp_path, assembly, bcftools_path, reference_genome_37_path):
    """
    Normaliza un archivo VCF de entrada utilizando bcftools.
    
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
        # Ruta del archivo de salida
        output_vcf_path = f"{temp_path}{input_vcf_path.split('/')[-1].split('.vcf')[0]}_normalized_int.vcf"
        output2_vcf_path = f"{temp_path}{input_vcf_path.split('/')[-1].split('.vcf')[0]}_normalized.vcf" 
        
        # Verificar si el archivo VCF está comprimido
        if input_vcf_path.endswith(".gz"):
            # Verificar si el archivo VCF necesita ser indexado
            if not (os.path.exists(input_vcf_path + ".csi") or os.path.exists(input_vcf_path + ".tbi")):
                # Indexar el archivo VCF si no está indexado
                index_command = [bcftools_path + "bcftools", "index", input_vcf_path]
                subprocess.run(index_command, check=True)
        
        # Comando para normalizar con bcftools
        if assembly == '37':
            bcftools_command = [bcftools_path + "bcftools", "norm", "-O", "v", "-m", "-any", "--check-ref", "w", "-f", reference_genome_37_path, "-o", output_vcf_path, input_vcf_path]
        elif assembly == '38':
            bcftools_command = [bcftools_path + "bcftools", "norm", "-O", "v", "-m", "-any", "--check-ref", "w", "-f", "./Homo_sapiens.GRCh38.dna.primary_assembly.fa", "-o", output_vcf_path, input_vcf_path]

        # Ejecutar el comando utilizando subprocess
        subprocess.run(bcftools_command, check=True)
        
        # Eliminar duplicados con bcftools
        rm_dup_command = [bcftools_path + "bcftools", "norm", "--rm-dup", "none", "-Oz", "-o", output2_vcf_path, output_vcf_path]
        subprocess.run(rm_dup_command, check=True)
        
        print("Normalización con bcftools completada.")
        
        return(output2_vcf_path)

    except Exception as e:
        print(f"Error durante la normalización con bcftools: {e}")

