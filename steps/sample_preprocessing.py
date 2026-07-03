# steps/sample_preprocessing_setup.py

from sftool.context import ExecutionContext
from sftool.misc.vcf_utils import normalize_vcf, intersect_vcf_with_bed
from sftool.misc.pharmcat_utils import pharmCAT_vcf_preprocessor


def run(ctx: ExecutionContext) -> None:
    """
    Sample preprocessing (VCF normalization and intersection with corresponding BED file)
    """

    for sample in ctx.samples:
        categories = sample.categories
        # 1. VCF normalization

        if 'PR' in categories or 'RR' in categories:
            norm_vcf_file = normalize_vcf(
                str(sample.vcf.resolve()),
                ctx.tmp_dir,
                ctx.config.paths.bcftools,
                ctx.config.references.genomes[ctx.assembly]
            )
            sample.vcf_outputs["normalized"] = norm_vcf_file

        # 2. Normalized VCF file intersection with BED file for each category

        for category in categories:
            if category == 'PR' or category == 'RR':
                sample.vcf_outputs["intersected"][category] = intersect_vcf_with_bed(norm_vcf_file, ctx.outputs["catalogs"]["bed_files"][category], ctx.tmp_dir, category)
            elif category == 'PGx' and ctx.assembly == 'GRCh38':
                # 1. Run pharmCAT's preprocessor script (https://pharmcat.org/using/VCF-Preprocessor/)
                pgx_vcf = sample.pgx_vcf or sample.vcf
                sample.vcf_outputs["PGx_preprocessed"] = pharmCAT_vcf_preprocessor(str(pgx_vcf.resolve()), ctx.config.paths.python, ctx.config.paths.pharmCAT, ctx.config.paths.bcftools, ctx.config.paths.bgzip, ctx.tmp_dir)








