

from modules.generate_report import get_hpos_from_txt
import os
import subprocess


def get_versions_paths(program_arguments, config_data, clinvar_db):

    # HPO
    if program_arguments.hpos_file is None:
        hpos_patient_list="Not provided"
    else:
        hpos_list_temp = get_hpos_from_txt(program_arguments.hpos_file)
        hpos_patient_list=",".join(str(hpo) for hpo in hpos_list_temp)

    # Intervar version
    try:
        intervar_file_path = os.path.join(config_data["intervar_path"], "Intervar.py")
        cmd = [intervar_file_path, "--version"]

        intervar_process = subprocess.Popen(cmd, stdout= subprocess.PIPE)
        intervar_out, intervar_err = intervar_process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error running InterVar: {e.output}")

    # Bcftools version
    try:
        bcftools_file_path = os.path.join(config_data["bcftools_path"], "bcftools")
        cmd = [bcftools_file_path, "--version"]

        bcftools_process = subprocess.Popen(cmd, stdout= subprocess.PIPE)
        bcftools_out, bcftools_err = bcftools_process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error running InterVar: {e.output}")



    versions_paths={
        "SF module": "0.1",
        "SF module mode": program_arguments.mode,
        "Input VCF": program_arguments.vcf_file,
        "Temporal dir": config_data["temp_path"],
        "Output dir": config_data["out_path"],
        "HPO list": hpos_patient_list,
        "Human assembly": "hg19" if program_arguments.assembly == 37 else "hg38",
        "Reference genome path":  config_data["reference_genome_37_path"] if program_arguments.assembly == 37 else config_data["reference_genome_38_path"],
        "Clinvar version": "Not used (basic mode)" if program_arguments.mode == "basic" else os.path.splitext(os.path.basename(clinvar_db))[0].split("_")[-1],
        "Clinvar path": clinvar_db,
        "Clinvar Evidence level": str(program_arguments.evidence),
        "Intervar version": str(intervar_out).split(" ")[1] + " " + str(intervar_out).split(" ")[2].split("\\n")[0],
        "Intervar path": config_data["intervar_path"],
        "bcftools version": str(bcftools_out).split(" ")[1].split("\\n")[0],
        "bcftools path": config_data["bcftools_path"],
        "HPO genes to phenotype version": os.path.splitext(os.path.basename(config_data["gene_to_phenotype_file"]))[0].split("_")[-1],
        "HPO genes to phenotype path": config_data["gene_to_phenotype_file"]

    }

    return versions_paths

