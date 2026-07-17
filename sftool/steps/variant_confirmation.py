from sftool.core.context import ExecutionContext
from sftool.variant_confirmation.candidate_vcf import (
    CandidateVcfNormalizer,
    CandidateVcfWriter,
)
from sftool.variant_confirmation.converter import GeneBeVariantConverter
from sftool.variant_confirmation.matcher import CandidateMatcher
from sftool.variant_confirmation.parser import VariantRepresentationParser
from sftool.variant_confirmation.result_builder import (
    VariantConfirmationResultBuilder,
)
from sftool.variant_confirmation.utils import (
    get_variant_confirmation_output_dir,
    has_detected_candidates,
    require_normalized_patient_vcf,
    run_variant_confirmation_genebe,
    store_variant_confirmation_outputs,
    write_conversion_json,
)


def run(ctx: ExecutionContext) -> None:
    """
    Run the Variant Confirmation workflow for each configured sample.

    This step assumes that:
    - ``variant_confirmation`` is enabled in ``ctx.modes``;
    - sample preprocessing has already generated the normalized patient VCF.
    """

    parser = VariantRepresentationParser()
    converter = GeneBeVariantConverter()
    writer = CandidateVcfWriter()

    normalizer = CandidateVcfNormalizer(
        bcftools_path=ctx.config.paths.bcftools,
    )

    matcher = CandidateMatcher(
        bcftools_path=ctx.config.paths.bcftools,
    )

    result_builder = VariantConfirmationResultBuilder()

    for sample in ctx.samples:
        request = sample.variant_confirmation_request

        if request is None:
            continue

        output_dir = get_variant_confirmation_output_dir(
            ctx=ctx,
            sample=sample,
        )

        normalized_patient_vcf = require_normalized_patient_vcf(
            sample
        )

        # 1. Detect representation type.
        parsed_request = parser.parse(
            request
        )

        # 2. Convert the representation into genomic candidate(s).
        candidates = converter.convert(
            request=parsed_request,
            assembly=ctx.assembly,
        )

        conversion_json = write_conversion_json(
            request=parsed_request,
            candidates=candidates,
            output_dir=output_dir,
        )

        # 3. Generate and normalize the candidate VCF.
        raw_candidate_vcf = writer.write_to_directory(
            candidates=candidates,
            output_directory=output_dir,
        )

        normalized_candidate_vcf = normalizer.normalize_to_directory(
            raw_vcf_path=raw_candidate_vcf,
            candidates=candidates,
            reference_fasta_path=(
                ctx.config.references.genomes[ctx.assembly]
            ),
            output_directory=output_dir,
        )

        # 4. Match normalized candidates against the normalized patient VCF.
        matching_output = matcher.match_to_directory(
            candidates=candidates,
            patient_vcf_path=normalized_patient_vcf,
            output_directory=output_dir,
        )

        # 5. Annotate only when at least one candidate was detected.
        genebe_annotated_vcf = None

        if has_detected_candidates(matching_output):
            genebe_annotated_vcf = run_variant_confirmation_genebe(
                ctx=ctx,
                input_vcf=matching_output.vcf_path,
                output_dir=output_dir,
            )

        # 6. Build and serialize the structured result.
        structured_output = result_builder.build_to_directory(
            request=parsed_request,
            candidates=candidates,
            matching_output=matching_output,
            annotated_vcf_path=genebe_annotated_vcf,
            output_directory=output_dir,
        )

        # 7. Store paths and in-memory result in SampleContext.
        store_variant_confirmation_outputs(
            sample=sample,
            conversion_json=conversion_json,
            raw_candidate_vcf=raw_candidate_vcf,
            normalized_candidate_vcf=normalized_candidate_vcf,
            matching_vcf=matching_output.vcf_path,
            genebe_annotated_vcf=genebe_annotated_vcf,
            result_json=structured_output.json_path,
            result=structured_output.result,
        )
