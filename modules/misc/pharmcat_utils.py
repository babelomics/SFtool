import subprocess
import os
import shutil
import json
from collections import defaultdict

def pharmCAT_vcf_preprocessor(vcf_input, python_path, pharmCAT_path, bgzip_path, htslib_path, tmp_dir):
    """
    Run pharmCAT's VCF preprocessor: https://pharmcat.org/using/VCF-Preprocessor/

    :param vcf_input:
    :param python_path:
    :param pharmCAT_path:
    :return:
    """

    preprocessed_vcf = ''

    try:
        # Run pharmCAT's VCF preprocessor script
        vcf_preprocessor_command = [python_path, os.path.dirname(pharmCAT_path) + "/pharmcat_vcf_preprocessor", "--path-to-bcftools", bgzip_path , "--path-to-bgzip", htslib_path , "-vcf", vcf_input]
        with subprocess.Popen(vcf_preprocessor_command, stderr=subprocess.STDOUT, text=True, cwd=os.path.dirname(pharmCAT_path)) as process:
            output, _ = process.communicate()

        preprocessed_vcf = vcf_input.split(".vcf.gz")[0] + ".preprocessed.vcf.bgz"
        if os.path.exists(preprocessed_vcf):
            file_name = os.path.basename(preprocessed_vcf)
            pgx_output_dir = os.path.join(str(tmp_dir), "PGx")
            os.makedirs(pgx_output_dir, exist_ok=True)
            shutil.move(preprocessed_vcf, os.path.join(pgx_output_dir, file_name))
            return os.path.join(pgx_output_dir, file_name)
        else:
            print("pharmCAT's preprocessed file does not exist. Exiting")
            exit(-1)

    except subprocess.CalledProcessError as e:
        print(f"Error when running pharmCAT's VCF preprocessor script: {e.output}")


def run_pharmCAT(preprocessed_vcf, pharmCAT_path, java_path, out_path):
    """
    Run pharmCAT

    :param preprocessed_vcf:
    :param out_path:
    :return:
    """

    try:
        pgx_output_dir = os.path.join(str(out_path),'PGx')
        os.makedirs(pgx_output_dir, exist_ok=True)

        pharmCAT_command = [java_path, "-jar", pharmCAT_path, "-vcf", preprocessed_vcf, "--output-dir", pgx_output_dir]
        with subprocess.Popen(pharmCAT_command, stderr=subprocess.STDOUT, text=True, cwd=os.path.dirname(pharmCAT_path)) as process:
            output, _ = process.communicate()

        file_name_prefix = os.path.basename(preprocessed_vcf).split(".preprocessed.vcf.bgz")[0]

        report_file = os.path.join(pgx_output_dir, file_name_prefix + ".report.html")
        phenotype_file = os.path.join(pgx_output_dir, file_name_prefix + ".phenotype.json")

        if os.path.exists(report_file) and os.path.exists(phenotype_file):
            return [report_file, phenotype_file]
        else:
            print("pharmCAT report could not be generated. Exiting")
            exit(-1)
    except subprocess.CalledProcessError as e:
        print(f"Error when running pharmCAT: {e.output}")


def pharmCAT_collection(phenotype_file):
    """
    Parse pharmCAT JSON file and return a data frame with basic information of genes, diplotypes and phenotypes

    :param phenotype_file:
    :return:
    """

    with open(phenotype_file, "r") as fd:
        data = json.load(fd)

    # Initialize data structure
    pgx_variants = defaultdict(list)

    gene_reports = data.get("geneReports", {})

    for gene, details in gene_reports.items():
        recommendation_diplotypes = details.get(
            "recommendationDiplotypes",
            []
        )

        grouped = defaultdict(list)

        for rec in recommendation_diplotypes:
            label = rec.get("label", "N/A")
            phenotypes = rec.get("phenotypes", [])

            if label in {"Unknown", "Unknown/Unknown"}:
                genotype = "Not determined"
                phenotype = "Not determined"
            else:
                genotype = label
                phenotype = (
                    ",".join(map(str, phenotypes))
                    if phenotypes
                    else "N/A"
                )

            grouped[phenotype].append(genotype)

        for phenotype, genotypes in grouped.items():
            pgx_variants[gene].append({
                "genotype": "; ".join(sorted(set(genotypes))),
                "phenotype": phenotype
            })

    return dict(pgx_variants)

