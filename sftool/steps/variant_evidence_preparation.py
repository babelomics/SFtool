
from sftool.core.context import ExecutionContext
from sftool.utils.geneBe_utils import run_genebe
from sftool.utils.pharmcat_utils import run_pharmCAT


def run(ctx: ExecutionContext) -> None:
    """
    Variant annotation with GeneBe and PharmCAT
    """

    genebe_apikey = ctx.config.genebe_credentials.api_key
    genebe_username = ctx.config.genebe_credentials.username
    genebe_path = ctx.config.paths.genebe
    java_path = ctx.config.paths.java
    pharmCAT_path = ctx.config.paths.pharmCAT
    assembly = ctx.assembly
    tmp_dir = ctx.tmp_dir

    unique_categories = sorted(
        {cat for sample in ctx.samples for cat in sample.categories}
    )

    # 1. GeneBe annotation and PharmCAT execution: specific for each sample
    for sample in ctx.samples:
        categories = sample.categories
        for category in categories:
            if category == 'PR' or category == 'RR':
                vcf_file = sample.vcf_outputs["intersected"][category]
                genebe_output_file = run_genebe(vcf_file, category, assembly, genebe_path, java_path, genebe_apikey, genebe_username, tmp_dir)
                sample.vcf_outputs["genebe_annotated"][category] = genebe_output_file
            elif category == 'PGx' and assembly == 'GRCh38': # Run PharmCAT
                [pharmCAT_report_file, pharmCAT_phenotype_file] = run_pharmCAT(sample.vcf_outputs["PGx_preprocessed"], pharmCAT_path, java_path, ctx.run_dir)
                sample.reports["PGx"] = pharmCAT_report_file
                sample.results["PGx"] = pharmCAT_phenotype_file
