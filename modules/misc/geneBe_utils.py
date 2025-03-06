# -*- coding: utf-8 -*-
"""

@author: jpflorido
"""
import subprocess
import os
import vcfpy
import gzip
import io

def run_genebe(norm_vcf, category, assembly, genebe_path, java_path, api_key, username):
    """
    Run GeneBe for annotate variants

    :param norm_vcf: Path to normalized file
    :param category: Gene category for annotation
    :param assembly: Reference genome version
    :param genebe_path: Path to geneBe annotator
    :param java_path: Path to Java
    :param api_key: Api key for annotating using GeneBe
    :param username: User name for annotating using GeneBe
    :return: Annotated VCF file
    """

    try:
        # Path to VCF intersected and output directory
        genebe_output_file = f"{norm_vcf.split('norm.' + category.upper() + '.vcf.gz')[0]}{category.upper()}{'.geneBe.vcf.gz'}"

        if assembly == '37':
            assembly_int = "hg19"
        elif assembly == '38':
            assembly_int = 'hg38'

        # Build command to run GeneBe
        genebe_file_path = os.path.join(genebe_path, "GeneBeClient.jar")
        cmd = [java_path,
               "-jar",
               genebe_file_path,
               "vcf",  "annotate",
               "--input-vcf", norm_vcf,
               "--output-vcf", genebe_output_file,
               "--genome", assembly_int,
               "--api-key", api_key,
               "--username", username
               ]

        # Run command and get output
        with subprocess.Popen(cmd, stderr=subprocess.STDOUT, text=True, cwd=genebe_path) as process:
            output, _ = process.communicate()


        return genebe_output_file

    except subprocess.CalledProcessError as e:
        print(f"Error when running Genebe: {e.output}")

def parse_genebe_output(genebe_output_vcf_file, mode):
    """
    Parse VCF file annotated by GeneBe

    :param genebe_output_vcf_file:
    :param mode:
    :return:
    """

    try:
        # Read VCF file

        genebe_results = {}
        #vcf_reader = vcfpy.Reader.from_path(genebe_output_vcf_file)
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
                    ref_gene = variant_record.INFO['gene_symbol_base']
                    genotype = variant_record.INFO['zygosity'][0]

                    if 'dbsnp_base' in variant_record.INFO:
                        rs = variant_record.INFO['dbsnp_base']
                    else:
                        rs = '.'
                    variant_consequence = variant_record.INFO['effect_base']
                    classification = variant_record.INFO['acmg_classification_base']

                    # Get only pathogenic and likely pathogenic variants or add them all if advanced (Clinvar) mode
                    if classification in ["Pathogenic", "Likely pathogenic"] or mode == 'advanced':
                        # Create a dictionary with interesting fields
                        genebe_results[variant] = {
                            "Gene": ref_gene,
                            "rs": rs,
                            "GeneBeClassification": classification,
                            "Genotype": genotype,
                            "Consequence": variant_consequence
                        }

        return genebe_results

    except Exception as e:
        raise Exception(f"Error when parsing Genebe annotated VCF file: {e}")
