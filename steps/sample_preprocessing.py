# steps/sample_preprocessing_setup.py

from modules.context import ExecutionContext
from modules.misc.vcf_utils import normalize_vcf, intersect_vcf_with_bed
from pathlib import Path


def run(ctx: ExecutionContext) -> None:
    """
    Sample preprocessing (VCF normalization and intersection with corresponding BED file)
    """

    for sample in ctx.samples:
        # 1. VCF normalization
        norm_vcf_file = normalize_vcf(
            str(sample.vcf.resolve()),
            ctx.tmp_dir,
            ctx.config.paths.bcftools,
            ctx.config.references.genomes[ctx.assembly]
        )

        sample.vcf_outputs["normalized"] = norm_vcf_file

        # 2. Normalized VCF file intersection with BED file for each category
        categories = sample.categories
        for category in categories:
            if category == 'PR' or category == 'RR':
                generated_vcf_file = intersect_vcf_with_bed(norm_vcf_file, ctx.outputs["catalogs"]["bed_files"][category], ctx.tmp_dir, category)
                sample.vcf_outputs["intersected"][category] = generated_vcf_file







