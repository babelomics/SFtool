import subprocess
import os

def pharmCAT_vcf_preprocessor(vcf_input, python_path, pharmCAT_path, bgzip_path, htslib_path):
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
        vcf_preprocessor_command = [python_path, os.path.dirname(pharmCAT_path) + "/pharmcat_vcf_preprocessor.py", "--path-to-bcftools", bgzip_path , "--path-to-bgzip", htslib_path , "-vcf", vcf_input]
        with subprocess.Popen(vcf_preprocessor_command, stderr=subprocess.STDOUT, text=True, cwd=os.path.dirname(pharmCAT_path)) as process:
            output, _ = process.communicate()

        preprocessed_vcf = vcf_input.split(".vcf.gz")[0] + ".preprocessed.vcf.bgz"
        if os.path.exists(preprocessed_vcf):
            return preprocessed_vcf
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
        pharmCAT_command = [java_path, "-jar", pharmCAT_path, "-vcf", preprocessed_vcf, "--output-dir", str(out_path)]
        with subprocess.Popen(pharmCAT_command, stderr=subprocess.STDOUT, text=True, cwd=os.path.dirname(pharmCAT_path)) as process:
            output, _ = process.communicate()

        file_name_prefix = os.path.basename(preprocessed_vcf).split(".preprocessed.vcf.bgz")[0]

        report_file = os.path.join(out_path, file_name_prefix + ".report.html")
        phenotype_file = os.path.join(out_path, file_name_prefix + ".phenotype.json")

        if os.path.exists(report_file) and os.path.exists(phenotype_file):
            return [report_file, phenotype_file]
        else:
            print("pharmCAT report could not be generated. Exiting")
            exit(-1)
    except subprocess.CalledProcessError as e:
        print(f"Error when running pharmCAT: {e.output}")
