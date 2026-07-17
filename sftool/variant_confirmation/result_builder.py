from __future__ import annotations

import json
import os
import re
import tempfile

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import vcfpy

from sftool.variant_confirmation.matcher import CandidateMatchingOutput
from sftool.variant_confirmation.models import (
    VariantCandidate,
    VariantConfirmationRequest,
    VariantConfirmationResult,
    VariantMatch,
)


class ResultBuilderError(RuntimeError):
    """
    Base exception for structured diagnostic variant result failures.
    """

    pass


class InvalidResultBuilderInputError(ResultBuilderError):
    """
    Raised when requests, candidates or matching results are inconsistent.
    """

    pass


class AnnotationParsingError(ResultBuilderError):
    """
    Raised when a GeneBe-annotated VCF cannot be parsed consistently.
    """

    pass


class ResultSerializationError(ResultBuilderError):
    """
    Raised when the structured result JSON cannot be written.
    """

    pass


class StructuredResultOutput:
    """
    Structured result generated for one sample.

    Contains the in-memory VariantConfirmationResult and the path to the
    serialized JSON file.
    """

    def __init__(
            self,
            result: VariantConfirmationResult,
            json_path: Path,
    ):
        if not isinstance(result, VariantConfirmationResult):
            raise TypeError(
                "result must be a VariantConfirmationResult"
            )

        if not isinstance(json_path, Path):
            raise TypeError(
                "json_path must be a pathlib.Path"
            )

        self.result = result
        self.json_path = json_path

    def to_dict(self) -> dict:
        return {
            "result": self.result.to_dict(),
            "json_path": str(self.json_path),
        }

    def __repr__(self) -> str:
        return (
            "StructuredResultOutput("
            f"status={self.result.get_status()!r}, "
            f"json_path={str(self.json_path)!r})"
        )


class VariantConfirmationResultBuilder:
    """
    Build and serialize the structured diagnostic variant confirmation result.

    This component combines an original request, normalized candidates,
    patient VCF matching results and optional GeneBe annotations.

    It does not run conversion, normalization, matching or GeneBe annotation.
    """

    DEFAULT_FILENAME = "diagnostic_variant_results.json"

    AMBIGUOUS_CONVERSION_WARNING = (
        "The input representation generated multiple genomic candidates"
    )

    NOT_FOUND_WARNING = (
        "None of the diagnostic variant candidates was found in the "
        "normalized patient VCF"
    )

    MISSING_GENOTYPE_WARNING = (
        "A detected diagnostic variant has no genotype information"
    )

    FILTER_WARNING = (
        "A detected diagnostic variant does not pass the patient VCF filters"
    )

    MISSING_ANNOTATION_WARNING = (
        "No GeneBe annotation was found for a detected diagnostic variant"
    )

    MULTIPLE_DETECTED_WARNING = (
        "More than one genomic candidate was detected for the input "
        "representation"
    )

    def build(
            self,
            request: VariantConfirmationRequest,
            candidates: Sequence[VariantCandidate],
            matching_output: CandidateMatchingOutput,
            annotated_vcf_path: str | Path | None = None,
    ) -> VariantConfirmationResult:
        """
        Build a complete structured result for one sample.
        """
        validated_request = self._validate_request(
            request
        )

        validated_candidates = self._validate_candidates(
            candidates
        )

        validated_matches = self._validate_matching_output(
            candidates=validated_candidates,
            matching_output=matching_output,
        )

        result = VariantConfirmationResult(
            request=validated_request
        )

        for candidate in validated_candidates:
            result.add_candidate(
                candidate
            )

        for variant_match in validated_matches:
            result.add_match(
                variant_match
            )

        if annotated_vcf_path is not None:
            annotations_by_key = self._read_annotations(
                annotated_vcf_path
            )

            self._attach_annotations(
                candidates=validated_candidates,
                matches=validated_matches,
                annotations_by_key=annotations_by_key,
            )

        self._add_result_warnings(
            result
        )

        result.validate_complete()

        return result

    def build_to_directory(
            self,
            request: VariantConfirmationRequest,
            candidates: Sequence[VariantCandidate],
            matching_output: CandidateMatchingOutput,
            output_directory: str | Path,
            annotated_vcf_path: str | Path | None = None,
            filename: str = DEFAULT_FILENAME,
    ) -> StructuredResultOutput:
        """
        Build the result and write it using the standard workflow filename.
        """
        output_directory = self._validate_output_directory(
            output_directory
        )

        filename = self._validate_filename(
            filename
        )

        result = self.build(
            request=request,
            candidates=candidates,
            matching_output=matching_output,
            annotated_vcf_path=annotated_vcf_path,
        )

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        json_path = (
                output_directory
                / filename
        )

        self._write_json_atomic(
            result=result,
            output_path=json_path,
        )

        return StructuredResultOutput(
            result=result,
            json_path=json_path,
        )

    @staticmethod
    def _validate_request(
            request: VariantConfirmationRequest,
    ) -> VariantConfirmationRequest:
        if not isinstance(request, VariantConfirmationRequest):
            raise InvalidResultBuilderInputError(
                "request must be a VariantConfirmationRequest"
            )

        return request

    @staticmethod
    def _validate_candidates(
            candidates: Sequence[VariantCandidate],
    ) -> List[VariantCandidate]:
        if isinstance(candidates, (str, bytes)):
            raise InvalidResultBuilderInputError(
                "candidates must be a sequence of VariantCandidate objects"
            )

        if not isinstance(candidates, Sequence):
            raise InvalidResultBuilderInputError(
                "candidates must be a sequence of VariantCandidate objects"
            )

        candidates = list(
            candidates
        )

        if not candidates:
            raise InvalidResultBuilderInputError(
                "At least one diagnostic variant candidate is required"
            )

        candidate_ids = set()
        matching_keys = set()

        for index, candidate in enumerate(
                candidates,
                start=1,
        ):
            if not isinstance(candidate, VariantCandidate):
                raise InvalidResultBuilderInputError(
                    f"Candidate {index} must be a VariantCandidate, "
                    f"got {type(candidate).__name__}"
                )

            if not candidate.is_normalized():
                raise InvalidResultBuilderInputError(
                    f"Candidate {candidate.candidate_id!r} "
                    "has not been normalized"
                )

            if candidate.candidate_id in candidate_ids:
                raise InvalidResultBuilderInputError(
                    "Candidate IDs must be unique; duplicated ID: "
                    f"{candidate.candidate_id}"
                )

            matching_key = candidate.get_matching_key()

            if matching_key in matching_keys:
                raise InvalidResultBuilderInputError(
                    "Normalized diagnostic candidates must be unique; "
                    "duplicated variant: "
                    f"{candidate.get_normalized_variant()}"
                )

            candidate_ids.add(
                candidate.candidate_id
            )

            matching_keys.add(
                matching_key
            )

        return candidates

    @staticmethod
    def _validate_matching_output(
            candidates: Sequence[VariantCandidate],
            matching_output: CandidateMatchingOutput,
    ) -> List[VariantMatch]:
        if not isinstance(matching_output, CandidateMatchingOutput):
            raise InvalidResultBuilderInputError(
                "matching_output must be a CandidateMatchingOutput"
            )

        matches = list(
            matching_output.matches
        )

        candidate_by_id = {
            candidate.candidate_id: candidate
            for candidate in candidates
        }

        matches_by_id: Dict[str, VariantMatch] = {}

        for index, variant_match in enumerate(
                matches,
                start=1,
        ):
            if not isinstance(variant_match, VariantMatch):
                raise InvalidResultBuilderInputError(
                    f"Match {index} must be a VariantMatch, "
                    f"got {type(variant_match).__name__}"
                )

            candidate_id = variant_match.candidate_id

            if candidate_id not in candidate_by_id:
                raise InvalidResultBuilderInputError(
                    "Match references an unknown candidate: "
                    f"{candidate_id}"
                )

            if candidate_id in matches_by_id:
                raise InvalidResultBuilderInputError(
                    "A match already exists for candidate: "
                    f"{candidate_id}"
                )

            candidate = candidate_by_id[candidate_id]

            if VariantConfirmationResultBuilder._get_match_key(
                    variant_match
            ) != candidate.get_matching_key():
                raise InvalidResultBuilderInputError(
                    "Match coordinates do not correspond to candidate "
                    f"{candidate_id}: {variant_match.get_variant()}"
                )

            matches_by_id[candidate_id] = variant_match

        candidates_without_match = [
            candidate.candidate_id
            for candidate in candidates
            if candidate.candidate_id not in matches_by_id
        ]

        if candidates_without_match:
            raise InvalidResultBuilderInputError(
                "Candidates without matching result: "
                + ", ".join(candidates_without_match)
            )

        return [
            matches_by_id[candidate.candidate_id]
            for candidate in candidates
        ]

    @staticmethod
    def _get_match_key(
            variant_match: VariantMatch,
    ) -> Tuple[str, int, str, str]:
        return (
            variant_match.chromosome,
            variant_match.position,
            variant_match.reference,
            variant_match.alternate,
        )

    def _read_annotations(
            self,
            annotated_vcf_path: str | Path,
    ) -> Dict[Tuple[str, int, str, str], List[dict]]:
        annotated_vcf_path = self._validate_annotated_vcf_path(
            annotated_vcf_path
        )

        annotations_by_key: Dict[
            Tuple[str, int, str, str],
            List[dict],
        ] = {}

        try:
            reader = vcfpy.Reader.from_path(
                str(annotated_vcf_path)
            )

            try:
                for record in reader:
                    alternate_values = [
                        alternate.value
                        for alternate in record.ALT
                        if alternate.value not in (None, ".", "*")
                    ]

                    for alternate_value in alternate_values:
                        record_key = self._get_record_key(
                            record=record,
                            alternate=str(alternate_value),
                        )

                        annotations = self._parse_record_annotations(
                            record
                        )

                        if annotations:
                            annotations_by_key.setdefault(
                                record_key,
                                [],
                            ).extend(
                                annotations
                            )

            finally:
                reader.close()

        except AnnotationParsingError:
            raise

        except Exception as exc:
            raise AnnotationParsingError(
                "Could not parse GeneBe-annotated VCF "
                f"{annotated_vcf_path}: {exc}"
            ) from exc

        return annotations_by_key

    def _parse_record_annotations(
            self,
            record: vcfpy.Record,
    ) -> List[dict]:
        raw_annotations = record.INFO.get(
            "acmg_by_gene_base"
        )

        raw_annotations = self._as_list(
            raw_annotations
        )

        if not raw_annotations:
            gene_symbols = self._as_list(
                record.INFO.get(
                    "gene_symbol_base"
                )
            )

            return [
                self._build_minimal_annotation(
                    gene_symbol=gene_symbol,
                    record=record,
                )
                for gene_symbol in gene_symbols
                if self._normalize_optional_string(gene_symbol) is not None
            ]

        annotations = []

        for raw_annotation in raw_annotations:
            if raw_annotation is None:
                continue

            raw_annotation = str(
                raw_annotation
            ).strip()

            if not raw_annotation or raw_annotation == ".":
                continue

            annotations.append(
                self._parse_genebe_annotation(
                    raw_annotation=raw_annotation,
                    record=record,
                )
            )

        return annotations

    def _parse_genebe_annotation(
            self,
            raw_annotation: str,
            record: vcfpy.Record,
    ) -> dict:
        """
        Parse the acmg_by_gene_base layout already used by SFtool.

        Expected positions:
          0 gene
          2 transcript
          3 consequence
          7 ACMG criteria
          8 ACMG classification
          9 HGVSc
         10 HGVSp
        """
        values = raw_annotation.split(
            "|"
        )

        if len(values) < 11:
            raise AnnotationParsingError(
                "Invalid acmg_by_gene_base entry at "
                f"{record.CHROM}:{record.POS}: "
                f"expected at least 11 fields, got {len(values)}"
            )

        return {
            "gene": self._normalize_optional_string(
                values[0]
            ),
            "transcript": self._normalize_optional_string(
                values[2]
            ),
            "consequence": self._normalize_optional_string(
                values[3]
            ),
            "hgvsc": self._normalize_optional_string(
                values[9]
            ),
            "hgvsp": self._normalize_optional_string(
                values[10]
            ),
            "dbsnp": self._get_first_info_value(
                record.INFO,
                "dbsnp_base",
            ),
            "acmg_classification": self._normalize_optional_string(
                values[8]
            ),
            "acmg_criteria": self._split_acmg_criteria(
                values[7]
            ),
        }

    def _build_minimal_annotation(
            self,
            gene_symbol: Any,
            record: vcfpy.Record,
    ) -> dict:
        return {
            "gene": self._normalize_optional_string(
                gene_symbol
            ),
            "transcript": None,
            "consequence": None,
            "hgvsc": None,
            "hgvsp": None,
            "dbsnp": self._get_first_info_value(
                record.INFO,
                "dbsnp_base",
            ),
            "acmg_classification": None,
            "acmg_criteria": [],
        }


    def _attach_annotations(
            self,
            candidates: Sequence[VariantCandidate],
            matches: Sequence[VariantMatch],
            annotations_by_key: Dict[
                Tuple[str, int, str, str],
                List[dict],
            ],
    ) -> None:
        candidate_by_id = {
            candidate.candidate_id: candidate
            for candidate in candidates
        }

        for variant_match in matches:
            if not variant_match.found:
                continue

            candidate = candidate_by_id[
                variant_match.candidate_id
            ]

            matching_key = candidate.get_matching_key()

            for annotation in annotations_by_key.get(
                    matching_key,
                    [],
            ):
                variant_match.add_annotation(
                    annotation
                )

    def _add_result_warnings(
            self,
            result: VariantConfirmationResult,
    ) -> None:
        detected_matches = result.get_detected_matches()

        if result.is_ambiguous():
            result.add_warning(
                self.AMBIGUOUS_CONVERSION_WARNING
            )

        if not detected_matches:
            result.add_warning(
                self.NOT_FOUND_WARNING
            )

        if len(detected_matches) > 1:
            result.add_warning(
                self.MULTIPLE_DETECTED_WARNING
            )

        if any(
                variant_match.genotype in (None, "", ".")
                for variant_match in detected_matches
        ):
            result.add_warning(
                self.MISSING_GENOTYPE_WARNING
            )

        if any(
                variant_match.passed_filter() is False
                for variant_match in detected_matches
        ):
            result.add_warning(
                self.FILTER_WARNING
            )

        if any(
                not variant_match.annotations
                for variant_match in detected_matches
        ):
            result.add_warning(
                self.MISSING_ANNOTATION_WARNING
            )

    @staticmethod
    def _get_record_key(
            record: vcfpy.Record,
            alternate: str,
    ) -> Tuple[str, int, str, str]:
        return (
            str(record.CHROM),
            int(record.POS),
            str(record.REF).upper(),
            str(alternate).upper(),
        )

    @staticmethod
    def _as_list(
            value: Any,
    ) -> List[Any]:
        if value is None:
            return []

        if isinstance(value, list):
            return value

        if isinstance(value, tuple):
            return list(
                value
            )

        return [
            value
        ]

    @staticmethod
    def _normalize_optional_string(
            value: Any,
    ) -> Optional[str]:
        if value is None:
            return None

        value = str(
            value
        ).strip()

        if value in ("", ".", "-", "NA"):
            return None

        return value

    @classmethod
    def _get_first_info_value(
            cls,
            info: Dict[str, Any],
            key: str,
    ) -> Optional[str]:
        values = cls._as_list(
            info.get(
                key
            )
        )

        for value in values:
            normalized_value = cls._normalize_optional_string(
                value
            )

            if normalized_value is not None:
                return normalized_value

        return None

    @classmethod
    def _get_first_available_info_value(
            cls,
            info: Dict[str, Any],
            keys: Sequence[str],
    ) -> Optional[str]:
        for key in keys:
            value = cls._get_first_info_value(
                info,
                key,
            )

            if value is not None:
                return value

        return None

    @staticmethod
    def _split_acmg_criteria(
            raw_criteria: Any,
    ) -> List[str]:
        raw_criteria = VariantConfirmationResultBuilder \
            ._normalize_optional_string(
            raw_criteria
        )

        if raw_criteria is None:
            return []

        criteria = [
            criterion.strip()
            for criterion in re.split(
                r"[,&;]+",
                raw_criteria,
            )
            if criterion.strip()
        ]

        return list(
            dict.fromkeys(
                criteria
            )
        )

    @staticmethod
    def _validate_annotated_vcf_path(
            annotated_vcf_path: str | Path,
    ) -> Path:
        if not isinstance(annotated_vcf_path, (str, Path)):
            raise TypeError(
                "annotated_vcf_path must be a string or pathlib.Path"
            )

        annotated_vcf_path = Path(
            annotated_vcf_path
        )

        if not annotated_vcf_path.is_file():
            raise FileNotFoundError(
                "GeneBe-annotated VCF not found: "
                f"{annotated_vcf_path}"
            )

        if not (
                annotated_vcf_path.name.endswith(".vcf")
                or annotated_vcf_path.name.endswith(".vcf.gz")
        ):
            raise AnnotationParsingError(
                "GeneBe-annotated VCF must use the .vcf or .vcf.gz extension"
            )

        return annotated_vcf_path

    @staticmethod
    def _validate_output_directory(
            output_directory: str | Path,
    ) -> Path:
        if not isinstance(output_directory, (str, Path)):
            raise TypeError(
                "output_directory must be a string or pathlib.Path"
            )

        return Path(
            output_directory
        )

    @staticmethod
    def _validate_filename(
            filename: str,
    ) -> str:
        if not isinstance(filename, str):
            raise TypeError(
                "filename must be a string"
            )

        filename = filename.strip()

        if not filename:
            raise ValueError(
                "filename must be a non-empty string"
            )

        if Path(filename).name != filename:
            raise ValueError(
                "filename must not contain directory components"
            )

        if not filename.endswith(".json"):
            raise ValueError(
                "Structured result filename must use the .json extension"
            )

        return filename

    @staticmethod
    def _write_json_atomic(
            result: VariantConfirmationResult,
            output_path: Path,
    ) -> None:
        temporary_path = None

        try:
            with tempfile.NamedTemporaryFile(
                    mode="w",
                    encoding="utf-8",
                    suffix=".tmp",
                    prefix=f".{output_path.name}.",
                    dir=output_path.parent,
                    delete=False,
            ) as output_handle:
                temporary_path = Path(
                    output_handle.name
                )

                json.dump(
                    result.to_dict(),
                    output_handle,
                    indent=2,
                    ensure_ascii=False,
                )

                output_handle.write(
                    "\n"
                )

            os.replace(
                temporary_path,
                output_path,
            )

        except Exception as exc:
            if temporary_path is not None:
                try:
                    temporary_path.unlink(
                        missing_ok=True
                    )
                except OSError:
                    pass

            raise ResultSerializationError(
                "Could not write structured diagnostic variant result "
                f"{output_path}: {exc}"
            ) from exc