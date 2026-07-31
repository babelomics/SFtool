from sftool.core.context import ExecutionContext
from sftool.utils.errors import ValidationError
from sftool.utils.vcf_utils import vcf_uses_chr_prefix


CATALOG_CATEGORIES = {"PR", "RR"}


def run(ctx: ExecutionContext) -> None:
    """
    Select installed catalog resources for the current execution.

    All sample VCFs participating in PR or RR analysis must use the same
    chromosome naming convention.
    """

    categories = sorted({
        category
        for sample in ctx.samples
        for category in sample.categories
        if category in CATALOG_CATEGORIES
    })

    if not categories:
        return

    relevant_samples = [
        sample
        for sample in ctx.samples
        if any(
            category in CATALOG_CATEGORIES
            for category in sample.categories
        )
    ]

    uses_chr_prefix = _resolve_execution_chr_prefix(
        ctx,
        relevant_samples,
    )

    for category in categories:
        ctx.resources.select_catalog_bed(
            category,
            uses_chr_prefix=uses_chr_prefix,
        )


def _resolve_execution_chr_prefix(
        ctx: ExecutionContext,
        samples: list,
) -> bool:
    conventions = {
        sample.sample_id: vcf_uses_chr_prefix(
            sample.vcf,
            bcftools_path=ctx.config.paths.bcftools,
        )
        for sample in samples
    }

    distinct_conventions = set(conventions.values())

    if len(distinct_conventions) != 1:
        prefixed_samples = sorted(
            sample_id
            for sample_id, uses_prefix in conventions.items()
            if uses_prefix
        )

        non_prefixed_samples = sorted(
            sample_id
            for sample_id, uses_prefix in conventions.items()
            if not uses_prefix
        )

        raise ValidationError(
            "All sample VCFs in one execution must use the same "
            "chromosome naming convention. "
            f"Samples using the 'chr' prefix: "
            f"{', '.join(prefixed_samples) or 'none'}. "
            f"Samples without the 'chr' prefix: "
            f"{', '.join(non_prefixed_samples) or 'none'}."
        )

    return next(iter(distinct_conventions))