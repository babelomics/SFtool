
from sftool.core.context import ExecutionContext
from sftool.utils.geneBe_utils import run_genebe
from sftool.utils.clinvar_utils import run_clinvar
from pathlib import Path
from sftool.utils.pharmcat_utils import run_pharmCAT


def run(ctx: ExecutionContext) -> None:
    """
    Variant annotation with GeneBe and Clinvar
    """

    genebe_apikey = ctx.config.genebe_credentials.api_key
    genebe_username = ctx.config.genebe_credentials.username
    genebe_path = ctx.config.paths.genebe
    java_path = ctx.config.paths.java
    pharmCAT_path = ctx.config.paths.pharmCAT
    assembly = ctx.assembly
    variant_classification_sources = ctx.variant_classification_sources
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

    # 2. Clinvar variant selection according to gene catalogs
    if 'clinvar' in variant_classification_sources and ('PR' in unique_categories or 'RR' in unique_categories):
        clinvar_evidence = ctx.clinvar_evidence
        clinvar_db = ctx.outputs["clinvar"]["clinvar_db"]
        clinvar_submission = ctx.outputs["clinvar"]["clinvar_summary_db"]

        for category in unique_categories:
            if category == 'PR' or category == 'RR':
                if category == 'PR':
                    category_geneset_file = ctx.config.catalogs.personal_risk_geneset
                elif category == 'RR':
                    category_geneset_file = ctx.config.catalogs.reproductive_risk_geneset

                category_tmp_dir = Path(tmp_dir) / category.upper()
                category_tmp_dir.mkdir(parents=True, exist_ok=True)

                clinvar_output_file = category_tmp_dir / f"clinvar.{category}.json"

                run_clinvar(clinvar_evidence, clinvar_db, clinvar_submission, category, category_geneset_file, clinvar_output_file)
                json_category = category + '_json'
                ctx.outputs["clinvar"][json_category] = clinvar_output_file