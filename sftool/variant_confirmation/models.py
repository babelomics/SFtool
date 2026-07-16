from __future__ import annotations

from typing import Any, Dict, List, Optional

class VariantConfirmationRequest:
    """
    Diagnostic variant representation provided by the user.

    The representation type is initially unknown and will be assigned
    by the representation parser.
    """

    SUPPORTED_REPRESENTATION_TYPES = {
        "genomic",
        "hgvsc",
        "hgvsg",
        "hgvsp",
    }

    def __init__(
            self,
            variant: str,
            representation_type: Optional[str] = None,
    ):
        self.variant = self._validate_variant(variant)
        self.representation_type = self._validate_representation_type(
            representation_type
        )

    @staticmethod
    def _validate_variant(variant: str) -> str:
        if not isinstance(variant, str):
            raise TypeError(
                "variant must be a string, "
                f"got {type(variant).__name__}"
            )

        variant = variant.strip()

        if not variant:
            raise ValueError("variant must be a non-empty string")

        return variant

    @classmethod
    def _validate_representation_type(
            cls,
            representation_type: Optional[str],
    ) -> Optional[str]:
        if representation_type is None:
            return None

        if representation_type not in cls.SUPPORTED_REPRESENTATION_TYPES:
            raise ValueError(
                "Unsupported variant representation type: "
                f"{representation_type}"
            )

        return representation_type

    @classmethod
    def from_dict(cls, data: dict) -> "VariantConfirmationRequest":
        if not isinstance(data, dict):
            raise TypeError(
                "VariantConfirmationRequest input must be a dictionary"
            )

        if set(data.keys()) != {"variant"}:
            raise ValueError(
                "VariantConfirmationRequest must contain "
                "exactly one field: 'variant'"
            )

        return cls(variant=data["variant"])

    def set_representation_type(self, representation_type: str):
        self.representation_type = self._validate_representation_type(
            representation_type
        )

    def to_dict(self) -> dict:
        return {
            "variant": self.variant,
            "representation_type": self.representation_type,
        }

    def __repr__(self) -> str:
        return (
            "VariantConfirmationRequest("
            f"variant={self.variant!r}, "
            f"representation_type={self.representation_type!r})"
        )



class VariantCandidate:
    """
    Genomic candidate generated from a diagnostic variant representation.

    Original coordinates correspond to the GeneBe conversion result.
    Normalized coordinates are populated after VCF normalization.
    """

    SUPPORTED_ASSEMBLIES = {
        "GRCh37",
        "GRCh38",
    }

    def __init__(
            self,
            candidate_id: str,
            chromosome: str,
            position: int,
            reference: str,
            alternate: str,
            assembly: str,
            conversion_warnings: Optional[List[str]] = None,
    ):
        self.candidate_id = self._validate_non_empty_string(
            candidate_id,
            "candidate_id",
        )

        self.chromosome = self._validate_non_empty_string(
            chromosome,
            "chromosome",
        )

        self.position = self._validate_position(
            position,
            "position",
        )

        self.reference = self._validate_allele(
            reference,
            "reference",
        )

        self.alternate = self._validate_allele(
            alternate,
            "alternate",
        )

        self.assembly = self._validate_assembly(assembly)

        self.normalized_chromosome: Optional[str] = None
        self.normalized_position: Optional[int] = None
        self.normalized_reference: Optional[str] = None
        self.normalized_alternate: Optional[str] = None

        self.conversion_warnings: List[str] = (
            list(conversion_warnings)
            if conversion_warnings
            else []
        )

    @staticmethod
    def _validate_non_empty_string(
            value: str,
            field_name: str,
    ) -> str:
        if not isinstance(value, str):
            raise TypeError(
                f"{field_name} must be a string, "
                f"got {type(value).__name__}"
            )

        value = value.strip()

        if not value:
            raise ValueError(
                f"{field_name} must be a non-empty string"
            )

        return value

    @staticmethod
    def _validate_position(
            position: int,
            field_name: str,
    ) -> int:
        if (
                not isinstance(position, int)
                or isinstance(position, bool)
                or position < 1
        ):
            raise ValueError(
                f"{field_name} must be a positive integer"
            )

        return position

    @classmethod
    def _validate_allele(
            cls,
            allele: str,
            field_name: str,
    ) -> str:
        allele = cls._validate_non_empty_string(
            allele,
            field_name,
        ).upper()

        if allele == ".":
            raise ValueError(
                f"{field_name} must contain a valid allele"
            )

        return allele

    @classmethod
    def _validate_assembly(cls, assembly: str) -> str:
        if assembly not in cls.SUPPORTED_ASSEMBLIES:
            raise ValueError(
                f"Unsupported reference genome: {assembly}"
            )

        return assembly

    def set_normalized_coordinates(
            self,
            chromosome: str,
            position: int,
            reference: str,
            alternate: str,
    ):
        """
        Assign the coordinates obtained after candidate VCF normalization.
        """
        self.normalized_chromosome = self._validate_non_empty_string(
            chromosome,
            "normalized_chromosome",
        )

        self.normalized_position = self._validate_position(
            position,
            "normalized_position",
        )

        self.normalized_reference = self._validate_allele(
            reference,
            "normalized_reference",
        )

        self.normalized_alternate = self._validate_allele(
            alternate,
            "normalized_alternate",
        )

    def add_warning(self, warning: str):
        warning = self._validate_non_empty_string(
            warning,
            "warning",
        )

        if warning not in self.conversion_warnings:
            self.conversion_warnings.append(warning)

    def is_normalized(self) -> bool:
        return all(
            value is not None
            for value in (
                self.normalized_chromosome,
                self.normalized_position,
                self.normalized_reference,
                self.normalized_alternate,
            )
        )

    def get_genomic_variant(self) -> str:
        return (
            f"{self.chromosome}:"
            f"{self.position}:"
            f"{self.reference}:"
            f"{self.alternate}"
        )

    def get_normalized_variant(self) -> Optional[str]:
        if not self.is_normalized():
            return None

        return (
            f"{self.normalized_chromosome}:"
            f"{self.normalized_position}:"
            f"{self.normalized_reference}:"
            f"{self.normalized_alternate}"
        )

    def get_matching_key(self) -> tuple:
        """
        Return CHROM, POS, REF and ALT used for patient VCF matching.

        Candidates must be normalized before matching.
        """
        if not self.is_normalized():
            raise ValueError(
                f"Candidate {self.candidate_id} has not been normalized"
            )

        return (
            self.normalized_chromosome,
            self.normalized_position,
            self.normalized_reference,
            self.normalized_alternate,
        )

    def to_dict(self) -> dict:
        return {
            "candidate_id": self.candidate_id,
            "assembly": self.assembly,
            "chromosome": self.chromosome,
            "position": self.position,
            "reference": self.reference,
            "alternate": self.alternate,
            "genomic_variant": self.get_genomic_variant(),
            "normalized_chromosome": self.normalized_chromosome,
            "normalized_position": self.normalized_position,
            "normalized_reference": self.normalized_reference,
            "normalized_alternate": self.normalized_alternate,
            "normalized_variant": self.get_normalized_variant(),
            "conversion_warnings": self.conversion_warnings.copy(),
        }

    def __repr__(self) -> str:
        return (
            "VariantCandidate("
            f"candidate_id={self.candidate_id!r}, "
            f"variant={self.get_genomic_variant()!r}, "
            f"assembly={self.assembly!r})"
        )


class VariantMatch:
    """
    Result of matching one normalized candidate against the normalized
    patient VCF.
    """

    def __init__(
            self,
            candidate_id: str,
            chromosome: str,
            position: int,
            reference: str,
            alternate: str,
            found: bool = False,
    ):
        self.candidate_id = VariantCandidate._validate_non_empty_string(
            candidate_id,
            "candidate_id",
        )

        self.chromosome = VariantCandidate._validate_non_empty_string(
            chromosome,
            "chromosome",
        )

        self.position = VariantCandidate._validate_position(
            position,
            "position",
        )

        self.reference = VariantCandidate._validate_allele(
            reference,
            "reference",
        )

        self.alternate = VariantCandidate._validate_allele(
            alternate,
            "alternate",
        )

        if not isinstance(found, bool):
            raise TypeError("found must be a boolean")

        self.found = found

        self.genotype: Optional[str] = None
        self.quality: Optional[float] = None
        self.filters: List[str] = []
        self.sample_format: Dict[str, Any] = {}

        # A variant may have annotations for multiple genes/transcripts.
        self.annotations: List[dict] = []

        self.warnings: List[str] = []

    def set_match_data(
            self,
            genotype: Optional[str],
            quality: Optional[float],
            filters: Optional[List[str]],
            sample_format: Optional[Dict[str, Any]],
    ):
        """
        Populate information extracted from the matching patient VCF record.
        """
        self.found = True
        self.genotype = genotype
        self.quality = quality
        self.filters = list(filters) if filters else []
        self.sample_format = (
            dict(sample_format)
            if sample_format
            else {}
        )

    def add_annotation(self, annotation: dict):
        if not isinstance(annotation, dict):
            raise TypeError("annotation must be a dictionary")

        self.annotations.append(dict(annotation))

    def add_warning(self, warning: str):
        warning = VariantCandidate._validate_non_empty_string(
            warning,
            "warning",
        )

        if warning not in self.warnings:
            self.warnings.append(warning)

    def get_variant(self) -> str:
        return (
            f"{self.chromosome}:"
            f"{self.position}:"
            f"{self.reference}:"
            f"{self.alternate}"
        )

    def passed_filter(self) -> Optional[bool]:
        if not self.found:
            return None

        if not self.filters:
            return True

        return self.filters == ["PASS"]

    def to_dict(self) -> dict:
        return {
            "candidate_id": self.candidate_id,
            "found": self.found,
            "chromosome": self.chromosome,
            "position": self.position,
            "reference": self.reference,
            "alternate": self.alternate,
            "variant": self.get_variant(),
            "genotype": self.genotype,
            "quality": self.quality,
            "filters": self.filters.copy(),
            "passed_filter": self.passed_filter(),
            "sample_format": dict(self.sample_format),
            "annotations": [
                dict(annotation)
                for annotation in self.annotations
            ],
            "warnings": self.warnings.copy(),
        }

    def __repr__(self) -> str:
        return (
            "VariantMatch("
            f"candidate_id={self.candidate_id!r}, "
            f"variant={self.get_variant()!r}, "
            f"found={self.found})"
        )



class VariantConfirmationResult:
    """
    Complete variant confirmation result for one sample.

    Maintains traceability from the original representation to all
    generated candidates and their patient VCF matching results.
    """

    def __init__(
            self,
            request: VariantConfirmationRequest,
    ):
        if not isinstance(request, VariantConfirmationRequest):
            raise TypeError(
                "request must be a VariantConfirmationRequest"
            )

        self.request = request
        self.candidates: List[VariantCandidate] = []
        self.matches: List[VariantMatch] = []
        self.warnings: List[str] = []

    def add_candidate(self, candidate: VariantCandidate):
        if not isinstance(candidate, VariantCandidate):
            raise TypeError(
                "candidate must be a VariantCandidate"
            )

        if self.get_candidate(candidate.candidate_id) is not None:
            raise ValueError(
                "Duplicated candidate_id: "
                f"{candidate.candidate_id}"
            )

        self.candidates.append(candidate)

    def add_match(self, match: VariantMatch):
        if not isinstance(match, VariantMatch):
            raise TypeError(
                "match must be a VariantMatch"
            )

        candidate = self.get_candidate(match.candidate_id)

        if candidate is None:
            raise ValueError(
                "Match references an unknown candidate: "
                f"{match.candidate_id}"
            )

        if self.get_match(match.candidate_id) is not None:
            raise ValueError(
                "A match already exists for candidate: "
                f"{match.candidate_id}"
            )

        self.matches.append(match)

    def get_candidate(
            self,
            candidate_id: str,
    ) -> Optional[VariantCandidate]:
        for candidate in self.candidates:
            if candidate.candidate_id == candidate_id:
                return candidate

        return None

    def get_match(
            self,
            candidate_id: str,
    ) -> Optional[VariantMatch]:
        for match in self.matches:
            if match.candidate_id == candidate_id:
                return match

        return None

    def add_warning(self, warning: str):
        warning = VariantCandidate._validate_non_empty_string(
            warning,
            "warning",
        )

        if warning not in self.warnings:
            self.warnings.append(warning)

    def get_detected_matches(self) -> List[VariantMatch]:
        return [
            match
            for match in self.matches
            if match.found
        ]

    def is_found(self) -> bool:
        return len(self.get_detected_matches()) > 0

    def is_ambiguous(self) -> bool:
        return len(self.candidates) > 1

    def get_status(self) -> str:
        if not self.is_found():
            return "not_found"

        if self.is_ambiguous():
            return "confirmed_ambiguous"

        return "confirmed"

    def is_complete(self) -> bool:
        """
        Every generated candidate must have a matching result, including
        candidates not detected in the patient VCF.
        """
        if not self.candidates:
            return False

        candidate_ids = {
            candidate.candidate_id
            for candidate in self.candidates
        }

        match_candidate_ids = {
            match.candidate_id
            for match in self.matches
        }

        return candidate_ids == match_candidate_ids

    def validate_complete(self):
        if not self.candidates:
            raise ValueError(
                "VariantConfirmationResult contains no candidates"
            )

        candidates_without_match = [
            candidate.candidate_id
            for candidate in self.candidates
            if self.get_match(candidate.candidate_id) is None
        ]

        if candidates_without_match:
            raise ValueError(
                "Candidates without matching result: "
                + ", ".join(candidates_without_match)
            )

    def to_dict(self) -> dict:
        self.validate_complete()

        detected_matches = self.get_detected_matches()

        return {
            "request": self.request.to_dict(),
            "status": self.get_status(),
            "found": self.is_found(),
            "ambiguous": self.is_ambiguous(),
            "candidate_count": len(self.candidates),
            "match_count": len(detected_matches),
            "candidates": [
                candidate.to_dict()
                for candidate in self.candidates
            ],
            "matches": [
                match.to_dict()
                for match in self.matches
            ],
            "warnings": self.warnings.copy(),
        }

    def __repr__(self) -> str:
        return (
            "VariantConfirmationResult("
            f"variant={self.request.variant!r}, "
            f"candidates={len(self.candidates)}, "
            f"matches={len(self.matches)}, "
            f"status={self.get_status()!r})"
        )
