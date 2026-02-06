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