# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 23:13:23 2023

@author: kindi
"""

from pybedtools import BedTool
import os

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



