from pathlib import Path
from types import SimpleNamespace

from sftool.steps import variant_confirmation
from sftool.variant_confirmation.models import (
    VariantCandidate,
    VariantConfirmationRequest,
    VariantConfirmationResult,
    VariantMatch,
)


class FakeMatchingOutput:
    """
    Minimal matching output used by the step integration tests.
    """

    def __init__(self, matches, vcf_path):
        self.matches = matches
        self.vcf_path = vcf_path


class FakeStructuredOutput:
    """
    Minimal structured result output used by the step integration tests.
    """

    def __init__(self, result, json_path):
        self.result = result
        self.json_path = json_path


def build_context(tmp_path, with_request=True):
    """
    Build the minimum ExecutionContext-like object required by the step.
    """
    normalized_patient_vcf = (
            tmp_path
            / "patient.normalized.vcf.gz"
    )
    normalized_patient_vcf.write_bytes(
        b"patient-vcf"
    )

    reference_fasta = (
            tmp_path
            / "GRCh38.fa"
    )
    reference_fasta.write_text(
        ">2\nA\n",
        encoding="utf-8",
    )

    sample = SimpleNamespace(
        sample_id="sample_1",
        variant_confirmation_request=(
            VariantConfirmationRequest(
                "NM_000251.3:c.2030C>A"
            )
            if with_request
            else None
        ),
        vcf_outputs={
            "normalized": normalized_patient_vcf,
            "variant_confirmation": {
                "raw": None,
                "normalized": None,
                "genebe_annotated": None,
                "matches": None,
            },
        },
        results={
            "variant_confirmation": None,
        },
    )

    ctx = SimpleNamespace(
        assembly="GRCh38",
        run_dir=tmp_path / "run",
        tmp_dir=tmp_path / "tmp",
        samples=[sample],
        config=SimpleNamespace(
            paths=SimpleNamespace(
                bcftools="/usr/bin/bcftools",
                genebe="/opt/genebe/GeneBeClient.jar",
                java="/usr/bin/java",
            ),
            references=SimpleNamespace(
                genomes={
                    "GRCh38": reference_fasta,
                }
            ),
            genebe_credentials=SimpleNamespace(
                api_key="test-api-key",
                username="test-user",
            ),
        ),
    )

    ctx.run_dir.mkdir()
    ctx.tmp_dir.mkdir()

    return ctx, sample


def install_step_mocks(
        monkeypatch,
        tmp_path,
        *,
        candidate_found=True,
):
    """
    Replace external workflow components while retaining the real step.
    """
    calls = []

    class FakeParser:
        def parse(self, request):
            calls.append("parse")
            request.set_representation_type(
                "hgvsc"
            )
            return request

    class FakeConverter:
        def convert(self, request, assembly):
            calls.append("convert")

            return [
                VariantCandidate(
                    candidate_id="candidate_1",
                    chromosome="2",
                    position=47476367,
                    reference="G",
                    alternate="T",
                    assembly=assembly,
                )
            ]

    class FakeWriter:
        def write_to_directory(
                self,
                candidates,
                output_directory,
        ):
            calls.append("write")

            output_path = (
                    Path(output_directory)
                    / "diagnostic_candidates.raw.vcf"
            )
            output_path.write_text(
                "##fileformat=VCFv4.2\n",
                encoding="utf-8",
            )

            return output_path

    class FakeNormalizer:
        def __init__(self, bcftools_path):
            calls.append("normalizer_init")

        def normalize_to_directory(
                self,
                raw_vcf_path,
                candidates,
                reference_fasta_path,
                output_directory,
        ):
            calls.append("normalize")

            candidates[0].set_normalized_coordinates(
                chromosome="2",
                position=47476367,
                reference="G",
                alternate="T",
            )

            output_path = (
                    Path(output_directory)
                    / "diagnostic_candidates.normalized.vcf.gz"
            )
            output_path.write_bytes(
                b"normalized-candidates"
            )

            return output_path

    class FakeMatcher:
        def __init__(self, bcftools_path):
            calls.append("matcher_init")

        def match_to_directory(
                self,
                candidates,
                patient_vcf_path,
                output_directory,
        ):
            calls.append("match")

            matching_vcf = (
                    Path(output_directory)
                    / "diagnostic_candidates.matches.vcf.gz"
            )
            matching_vcf.write_bytes(
                b"matching-vcf"
            )

            variant_match = VariantMatch(
                candidate_id="candidate_1",
                chromosome="2",
                position=47476367,
                reference="G",
                alternate="T",
                found=candidate_found,
            )

            return FakeMatchingOutput(
                matches=[variant_match],
                vcf_path=matching_vcf,
            )

    class FakeResultBuilder:
        def build_to_directory(
                self,
                request,
                candidates,
                matching_output,
                annotated_vcf_path,
                output_directory,
        ):
            calls.append("build_result")

            result = VariantConfirmationResult(
                request=request
            )

            for candidate in candidates:
                result.add_candidate(
                    candidate
                )

            for variant_match in matching_output.matches:
                result.add_match(
                    variant_match
                )

            result_path = (
                    Path(output_directory)
                    / "diagnostic_variant_results.json"
            )
            result_path.write_text(
                "{}\n",
                encoding="utf-8",
            )

            return FakeStructuredOutput(
                result=result,
                json_path=result_path,
            )

    def fake_get_output_dir(ctx, sample):
        calls.append("get_output_dir")

        output_dir = (
                Path(ctx.run_dir)
                / sample.sample_id
                / "variant_confirmation"
        )
        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        return output_dir

    def fake_require_normalized_patient_vcf(sample):
        calls.append("require_patient_vcf")
        return Path(
            sample.vcf_outputs["normalized"]
        )

    def fake_write_conversion_json(
            request,
            candidates,
            output_dir,
    ):
        calls.append("write_conversion_json")

        output_path = (
                Path(output_dir)
                / "conversion.json"
        )
        output_path.write_text(
            "{}\n",
            encoding="utf-8",
        )

        return output_path

    def fake_has_detected_candidates(matching_output):
        calls.append("has_detected_candidates")

        return any(
            match.found
            for match in matching_output.matches
        )

    def fake_run_variant_confirmation_genebe(
            ctx,
            input_vcf,
            output_dir,
    ):
        calls.append("run_genebe")

        output_path = (
                Path(output_dir)
                / "diagnostic_candidates.genebe.vcf.gz"
        )
        output_path.write_bytes(
            b"annotated-vcf"
        )

        return output_path

    def fake_store_outputs(
            sample,
            conversion_json,
            raw_candidate_vcf,
            normalized_candidate_vcf,
            matching_vcf,
            genebe_annotated_vcf,
            result_json,
            result,
    ):
        calls.append("store_outputs")

        sample.vcf_outputs[
            "variant_confirmation"
        ].update(
            {
                "conversion": conversion_json,
                "raw": raw_candidate_vcf,
                "normalized": normalized_candidate_vcf,
                "matches": matching_vcf,
                "genebe_annotated": genebe_annotated_vcf,
                "result_json": result_json,
            }
        )

        sample.results[
            "variant_confirmation"
        ] = result

    monkeypatch.setattr(
        variant_confirmation,
        "VariantRepresentationParser",
        FakeParser,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "GeneBeVariantConverter",
        FakeConverter,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "CandidateVcfWriter",
        FakeWriter,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "CandidateVcfNormalizer",
        FakeNormalizer,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "CandidateMatcher",
        FakeMatcher,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "VariantConfirmationResultBuilder",
        FakeResultBuilder,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "get_variant_confirmation_output_dir",
        fake_get_output_dir,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "require_normalized_patient_vcf",
        fake_require_normalized_patient_vcf,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "write_conversion_json",
        fake_write_conversion_json,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "has_detected_candidates",
        fake_has_detected_candidates,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "run_variant_confirmation_genebe",
        fake_run_variant_confirmation_genebe,
    )
    monkeypatch.setattr(
        variant_confirmation,
        "store_variant_confirmation_outputs",
        fake_store_outputs,
    )

    return calls


def test_run_executes_complete_workflow_for_configured_sample(
        tmp_path,
        monkeypatch,
):
    """
    A configured sample must execute every Variant Confirmation stage in
    order and store all generated outputs in SampleContext.
    """
    ctx, sample = build_context(
        tmp_path
    )

    calls = install_step_mocks(
        monkeypatch,
        tmp_path,
        candidate_found=True,
    )

    variant_confirmation.run(
        ctx
    )

    assert calls == [
        "normalizer_init",
        "matcher_init",
        "get_output_dir",
        "require_patient_vcf",
        "parse",
        "convert",
        "write_conversion_json",
        "write",
        "normalize",
        "match",
        "has_detected_candidates",
        "run_genebe",
        "build_result",
        "store_outputs",
    ]

    outputs = sample.vcf_outputs[
        "variant_confirmation"
    ]

    assert outputs["conversion"].name == "conversion.json"
    assert outputs["raw"].name == "diagnostic_candidates.raw.vcf"
    assert (
            outputs["normalized"].name
            == "diagnostic_candidates.normalized.vcf.gz"
    )
    assert (
            outputs["matches"].name
            == "diagnostic_candidates.matches.vcf.gz"
    )
    assert (
            outputs["genebe_annotated"].name
            == "diagnostic_candidates.genebe.vcf.gz"
    )
    assert (
            outputs["result_json"].name
            == "diagnostic_variant_results.json"
    )
    assert sample.results["variant_confirmation"] is not None


def test_run_skips_sample_without_variant_confirmation_request(
        tmp_path,
        monkeypatch,
):
    """
    A sample without a variant_confirmation request must be ignored by the
    step without generating outputs.
    """
    ctx, sample = build_context(
        tmp_path,
        with_request=False,
    )

    calls = install_step_mocks(
        monkeypatch,
        tmp_path,
    )

    variant_confirmation.run(
        ctx
    )

    assert calls == [
        "normalizer_init",
        "matcher_init",
    ]
    assert sample.results["variant_confirmation"] is None
    assert (
            sample.vcf_outputs["variant_confirmation"]["raw"]
            is None
    )


def test_run_skips_genebe_when_candidate_is_not_detected(
        tmp_path,
        monkeypatch,
):
    """
    GeneBe annotation must not run when none of the converted candidates is
    present in the normalized patient VCF.
    """
    ctx, sample = build_context(
        tmp_path
    )

    calls = install_step_mocks(
        monkeypatch,
        tmp_path,
        candidate_found=False,
    )

    variant_confirmation.run(
        ctx
    )

    assert "run_genebe" not in calls
    assert (
            sample.vcf_outputs[
                "variant_confirmation"
            ]["genebe_annotated"]
            is None
    )
    assert sample.results["variant_confirmation"] is not None
