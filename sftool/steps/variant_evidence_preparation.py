
from sftool.core.context import ExecutionContext
from sftool.utils.geneBe_utils import run_genebe
from sftool.utils.pharmcat_utils import run_pharmCAT


def run(ctx: ExecutionContext) -> None:
    """
    Variant evidence preparation (variant annotation) with GeneBe and PharmCAT
    """
    # GeneBe annotation and PharmCAT execution: specific for each sample
    for sample in ctx.samples:
        categories = sample.categories
        for category in categories:
            if category == 'PR' or category == 'RR':
                vcf_file = sample.vcf_outputs["intersected"][category]
                genebe_output_file = run_genebe(
                    norm_vcf=vcf_file,
                    category=category,
                    assembly=ctx.assembly,
                    genebe_path=ctx.config.paths.genebe,
                    java_path=ctx.config.paths.java,
                    api_key=ctx.config.genebe_credentials.api_key,
                    username=ctx.config.genebe_credentials.username,
                    tmp_dir=ctx.tmp_dir
                )
                sample.vcf_outputs["genebe_annotated"][category] = genebe_output_file
            elif category == 'PGx' and ctx.assembly == 'GRCh38': # Run PharmCAT
                [pharmCAT_report_file, pharmCAT_phenotype_file] = run_pharmCAT(
                    preprocessed_vcf=sample.vcf_outputs["PGx_preprocessed"],
                    pharmCAT_path=ctx.config.paths.pharmCAT,
                    java_path=ctx.config.paths.java,
                    out_path=ctx.run_dir,
                )
                sample.reports["PGx"] = pharmCAT_report_file
                sample.results["PGx"] = pharmCAT_phenotype_file
