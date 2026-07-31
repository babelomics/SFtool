# steps/sample_preprocessing_setup.py

from sftool.core.context import ExecutionContext
from sftool.utils.vcf_utils import normalize_vcf, intersect_vcf_with_bed, check_vcf_positions_present
from sftool.utils.pharmcat_utils import pharmCAT_vcf_preprocessor

CATALOG_CATEGORIES = {"PR", "RR"}


def run(ctx: ExecutionContext) -> None:
    """
    Sample preprocessing (VCF normalization and intersection with corresponding BED file)
    """

    for sample in ctx.samples:
        categories = sample.categories
        # 1. VCF normalization

        if 'PR' in categories or 'RR' in categories or sample.variant_confirmation_request is not None:
            norm_vcf_file = normalize_vcf(
                input_vcf_path=str(sample.vcf.resolve()),
                temp_path=ctx.tmp_dir,
                bcftools_path=ctx.config.paths.bcftools,
                reference_genome_path=ctx.resources.reference_genome
            )
            sample.vcf_outputs["normalized"] = norm_vcf_file

        # 2. Normalized VCF file intersection with BED file for each category

        for category in categories:
            if category in CATALOG_CATEGORIES:
                sample.vcf_outputs["intersected"][category] = intersect_vcf_with_bed(norm_vcf_file, ctx.resources.selected_catalog_bed(category), ctx.tmp_dir, category)
            elif category == 'PGx' and ctx.assembly == 'GRCh38':
                # 1. Run pharmCAT's preprocessor script (https://pharmcat.org/using/VCF-Preprocessor/)
                pgx_vcf = sample.pgx_vcf or sample.vcf

                sample.vcf_outputs["PGx_preprocessed"] = (
                    pharmCAT_vcf_preprocessor(
                        vcf_input=pgx_vcf,
                        python_path=ctx.config.paths.python,
                        pharmcat_path=ctx.config.paths.pharmCAT,
                        bcftools_path=ctx.config.paths.bcftools,
                        bgzip_path=ctx.config.paths.bgzip,
                        tmp_dir=ctx.tmp_dir,
                    )
                )







