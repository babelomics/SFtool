# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 22:15:15 2023

@author: Edurne Urrutia, Javier Perez Florido
"""
import os
import subprocess
from pybedtools import BedTool
import logging
from pathlib import Path

from sftool.utils.errors import ValidationError

logger = logging.getLogger(__name__)

#####################################################################################
# Preprocessing functions
#####################################################################################


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
            #index_command = [bcftools_path + "bcftools", "index", input_vcf_path]
            index_command = [bcftools_path, "index", input_vcf_path]
            subprocess.run(index_command, check=True, capture_output=True)

        # Output files
        output_vcf_path = os.path.join(temp_path, just_filename.split(".vcf.gz")[0] + ".tmp.vcf.gz")
        output2_vcf_path = os.path.join(temp_path, just_filename.split(".vcf.gz")[0] + ".tmp2.vcf.gz")
        output3_vcf_path = os.path.join(temp_path, just_filename.split(".vcf.gz")[0] + ".norm.vcf.gz")

        # bcftools normalization command
        bcftools_command = [bcftools_path, "norm", "-O", "z", "-m", "-any", "--check-ref", "w",  "-f", reference_genome_path, "-o", output_vcf_path, input_vcf_path]

        # Normalize with bcftools
        subprocess.run(bcftools_command, check=True)

        # Remove duplicates with bcftools
        rm_dup_command = [bcftools_path, "norm", "--rm-dup", "none", "-Oz", "-o", output2_vcf_path, output_vcf_path]
        subprocess.run(rm_dup_command, check=True)

        print("bcftools normalization completed.")

        # Remove non-variant sites (genotypes with 0/0)
        rm_nonvariantsites_command = [bcftools_path, "view", "-e", 'ALT="*" || GT="0/0"', "-Oz", "-o", output3_vcf_path, output2_vcf_path]
        subprocess.run(rm_nonvariantsites_command, check=True)

        print("bcftools filtering non variant sites completed.")

        os.remove(output_vcf_path)
        os.remove(output2_vcf_path)
        return(output3_vcf_path)

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
        output_dir = os.path.join(str(temp_path), category)
        os.makedirs(output_dir, exist_ok=True)
        output_vcf_path = os.path.join(output_dir, just_filename.split(".vcf.gz")[0] + "." + category.upper() + ".vcf.gz")

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

# ==========================================================
# Generic helpers
# ==========================================================

def run_command(cmd: list[str]) -> str:
    """
    Execute a shell command and return stdout.
    Raise ValidationError on failure.
    """
    try:
        return subprocess.check_output(
            cmd,
            text=True,
            stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise ValidationError(
            f"Command failed:\n{' '.join(cmd)}\n\n{e.stderr}"
        ) from e


def get_vcf_contigs(vcf_path: Path) -> dict[str, int]:
    """
    Return contigs declared in the VCF header.

    Returns
    -------
    dict
        {contig_name: contig_length}
    """
    header = run_command(["bcftools", "view", "-h", str(vcf_path)])

    contigs = {}

    for line in header.splitlines():
        if not line.startswith("##contig=<"):
            continue

        content = line.removeprefix("##contig=<").removesuffix(">")

        fields = dict(
            item.split("=", 1)
            for item in content.split(",")
            if "=" in item
        )

        if "ID" in fields and "length" in fields:
            contigs[fields["ID"]] = int(fields["length"])

    return contigs





def get_vcf_positions(vcf_path: Path) -> set[tuple[str, int]]:
    """
    Return all (chromosome, position) pairs in a VCF.
    """
    output = run_command([
        "bcftools",
        "query",
        "-f",
        "%CHROM\t%POS\n",
        str(vcf_path)
    ])

    return {
        (chrom, int(pos))
        for chrom, pos in (
            line.split("\t")
            for line in output.splitlines()
            if line
        )
    }


# ==========================================================
# Validation helpers
# ==========================================================


def validate_chr_prefix(vcf_path: Path):
    """
    Validate that all VCF contigs use the 'chr' prefix.
    """

    contigs = get_vcf_contigs(vcf_path)

    invalid = [
        chrom
        for chrom in contigs
        if not chrom.startswith("chr")
    ]

    if invalid:
        raise ValidationError(
            "VCF chromosome names must use the 'chr' prefix. "
            f"Examples: {', '.join(invalid[:10])}"
        )


def check_vcf_positions_present(
        input_vcf: Path,
        required_vcf: Path,
        output_file: Path | None = None
)-> set[tuple[str, int]]:
    """
    Validate that all positions present in required_vcf
    also exist in input_vcf.
    """

    required = get_vcf_positions(required_vcf)
    observed = get_vcf_positions(input_vcf)

    missing = required - observed

    if missing:
        examples = sorted(missing)[:10]

        logger.warning(
            "%d PharmCAT positions are missing from %s. "
            "PharmCAT may still run successfully. "
            "First missing positions: %s",
            len(missing),
            input_vcf,
            ", ".join(f"{chrom}:{pos}" for chrom, pos in examples)
        )

        if output_file:
            output_file.parent.mkdir(parents=True, exist_ok=True)

            with output_file.open("w") as fh:
                fh.write("CHROM\tPOS\n")

                for chrom, pos in sorted(missing):
                    fh.write(f"{chrom}\t{pos}\n")

    return missing