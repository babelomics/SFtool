import subprocess
import os
import shutil
import json
from collections import defaultdict
from pathlib import Path


def pharmCAT_vcf_preprocessor(
        *,
        vcf_input: Path | str,
        python_path: Path | str,
        pharmcat_path: Path | str,
        bcftools_path: Path | str,
        bgzip_path: Path | str,
        tmp_dir: Path | str,
) -> Path:
    """
    Run the PharmCAT VCF preprocessor (https://pharmcat.org/using/VCF-Preprocessor/) and return the generated VCF.

    Executable paths are provided explicitly by the runtime
    configuration.
    """

    vcf_input = Path(vcf_input).resolve()
    pharmcat_path = Path(pharmcat_path).resolve()
    tmp_dir = Path(tmp_dir).resolve()

    if not vcf_input.is_file():
        raise ResourceOperationError(
            f"PharmCAT input VCF does not exist: "
            f"{vcf_input}"
        )

    preprocessor_path = (
            pharmcat_path.parent
            / "pharmcat_vcf_preprocessor"
    )

    if not preprocessor_path.is_file():
        raise ResourceOperationError(
            "PharmCAT VCF preprocessor script does not "
            f"exist: {preprocessor_path}"
        )

    command = [
        str(python_path),
        str(preprocessor_path),
        "--path-to-bcftools",
        str(bcftools_path),
        "--path-to-bgzip",
        str(bgzip_path),
        "-vcf",
        str(vcf_input),
    ]

    result = subprocess.run(
        command,
        cwd=pharmcat_path.parent,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        raise ResourceOperationError(
            "PharmCAT VCF preprocessing failed "
            f"for {vcf_input}:\n{result.stdout}"
        )

    input_name = vcf_input.name

    if input_name.endswith(".vcf.gz"):
        output_name = (
                input_name[:-len(".vcf.gz")]
                + ".preprocessed.vcf.bgz"
        )
    elif input_name.endswith(".vcf.bgz"):
        output_name = (
                input_name[:-len(".vcf.bgz")]
                + ".preprocessed.vcf.bgz"
        )
    elif input_name.endswith(".vcf"):
        output_name = (
                input_name[:-len(".vcf")]
                + ".preprocessed.vcf.bgz"
        )
    else:
        raise ResourceOperationError(
            f"Unsupported PharmCAT input VCF name: "
            f"{input_name}"
        )

    generated_vcf = (
            vcf_input.parent
            / output_name
    )

    if not generated_vcf.is_file():
        raise ResourceOperationError(
            "PharmCAT VCF preprocessing completed but "
            "the expected output was not generated: "
            f"{generated_vcf}"
        )

    pgx_tmp_dir = tmp_dir / "PGx"
    pgx_tmp_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    destination = (
            pgx_tmp_dir
            / output_name
    )

    if destination.exists():
        destination.unlink()

    shutil.move(
        str(generated_vcf),
        str(destination),
    )

    return destination


def run_pharmCAT(
        *,
        preprocessed_vcf: Path | str,
        pharmCAT_path: Path | str,
        java_path: Path | str,
        out_path: Path | str,
) -> tuple[Path, Path]:
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

