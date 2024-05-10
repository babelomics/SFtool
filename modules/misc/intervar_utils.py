# -*- coding: utf-8 -*-
"""

@author: jpflorido
"""
import subprocess
import os

def run_intervar(norm_vcf, category, assembly, intervar_path):
    """
    Run Intervar to annotate variants

    Args:
        norm_vcf (str): Path to normalized file
        category (str): Gene category for annotation
        assembly (str): Reference genome version
        intervar_path (str): Path to InterVar directory

    Raises:
        subprocess.CalledProcessError: If an error arises when running Intervar
    """
    try:
        # Path to VCF intersected and output directory
        output_file = f"{norm_vcf.split('norm.' + category.upper() + '.vcf.gz')[0]}{category.upper()}"

        if assembly == '37':
            assembly_int = "hg19"
        elif assembly == '38':
            assembly_int = 'hg38'

        # Build command to run Intervar

        intervar_file_path = os.path.join(intervar_path, "Intervar.py")
        cmd = ["python3",
            intervar_file_path,
            "-d", intervar_path + '/humandb/',
            "-b", assembly_int,
            "-i", norm_vcf,
            "--input_type", "VCF",
            "-o", output_file
        ]

        # Run command and get output
        with subprocess.Popen(cmd, stderr=subprocess.STDOUT, text=True, cwd=intervar_path) as process:
            output, _ = process.communicate()

        intervar_output_file = f"{norm_vcf.split('norm.' + category.upper() + '.vcf.gz')[0]}{category.upper()}.{assembly_int}_multianno.txt.intervar"
        return intervar_output_file

    except subprocess.CalledProcessError as e:
        print(f"Error when running Intervar: {e.output}")

def parse_intervar_output(intervar_output_file, mode, assembly):
    """
    Parse Intervar output file and get interesting fields

    Args:
        intervar_output_file (str): Path to intervar output file
        mode (str): basic or advanced
        assembly (str): Reference genome version

    Returns:
        list: A list of dictionaries with interesting fields
    """

    if assembly == '37':
        assembly_int = "hg19"
    elif assembly == '38':
        assembly_int = 'hg38'

    intervar_results = {}

    with open(intervar_output_file, "r") as intervar_file:
        for line in intervar_file:
            # Ignore lines that are not variants
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                # Get interesting fields and format them
                variant = f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
                ref_gene = fields[5]
                intervar_consequence = fields[6] + ',' + fields[7]
                avsnp = fields[9]
                classification = fields[13].split(": ")[1].split(" PVS")[0]
                orpha = fields[32]
                other_info = fields[-1]

                # Get only pathogenic and likely pathogenic variants
                if classification in ["Pathogenic", "Likely pathogenic"] or mode == 'advanced':

                    # Create a dictionary with interesting fields
                    intervar_results[variant] = {
                        "Gene": ref_gene,
                        "rs": avsnp,
                        "IntervarClassification": classification,
                        "Orpha": orpha,
                        "Genotype": other_info,
                        "IntervarConsequence": intervar_consequence
                    }

    return intervar_results