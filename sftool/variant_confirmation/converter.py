# sftool/variant_confirmation/converter.py

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import requests

from sftool.variant_confirmation.models import (
    VariantCandidate,
    VariantConfirmationRequest,
)


class VariantConversionError(RuntimeError):
    """
    Base exception for blocking diagnostic variant conversion failures.
    """

    pass


class GeneBeRequestError(VariantConversionError):
    """
    Raised when the GeneBe service cannot be reached or returns an HTTP error.
    """

    pass


class GeneBeResponseError(VariantConversionError):
    """
    Raised when GeneBe returns malformed or unexpected response data.
    """

    pass


class GeneBeNoCandidatesError(VariantConversionError):
    """
    Raised when GeneBe cannot generate genomic candidates for a representation.
    """

    pass


class GeneBeAssemblyMismatchError(VariantConversionError):
    """
    Raised when GeneBe returns coordinates for a different genome assembly.
    """

    pass


class GeneBeVariantConverter:
    """
    Convert one diagnostic variant representation into genomic candidates
    using the public GeneBe conversion API.

    The converter performs coordinate conversion only. It does not generate
    VCF files, normalize variants, match variants against the patient VCF or
    annotate variants.
    """

    API_URL = "https://api.genebe.net/cloud/api-public/v1/convert"

    ASSEMBLY_TO_GENEBE = {
        "GRCh37": "hg19",
        "GRCh38": "hg38",
    }

    GENEBE_TO_ASSEMBLY = {
        "hg19": "GRCh37",
        "hg38": "GRCh38",
    }

    def __init__(
            self,
            session: Optional[requests.Session] = None,
            timeout: float = 30.0,
    ):
        """
        Parameters
        ----------
        session
            Optional requests-compatible session. Dependency injection allows
            the HTTP client to be mocked during unit tests.
        timeout
            Maximum number of seconds allowed for the HTTP request.
        """
        if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
            raise TypeError("timeout must be a positive number")

        if timeout <= 0:
            raise ValueError("timeout must be greater than zero")

        self.session = session or requests.Session()
        self.timeout = float(timeout)

    def convert(
            self,
            request: VariantConfirmationRequest,
            assembly: str,
    ) -> List[VariantCandidate]:
        """
        Convert a parsed diagnostic variant representation into one or more
        genomic candidates.

        Parameters
        ----------
        request
            Variant confirmation request whose representation type has
            already been assigned by VariantRepresentationParser.
        assembly
            SFtool reference genome: GRCh37 or GRCh38.

        Returns
        -------
        List[VariantCandidate]
            Unique genomic candidates in the order returned by GeneBe.

        Raises
        ------
        TypeError
            If request is not a VariantConfirmationRequest.
        ValueError
            If request has not been parsed or assembly is unsupported.
        GeneBeRequestError
            If the external service cannot be reached or returns an HTTP error.
        GeneBeResponseError
            If the service response is malformed.
        GeneBeNoCandidatesError
            If no genomic candidates are returned.
        GeneBeAssemblyMismatchError
            If the returned genome differs from the requested assembly.
        """
        self._validate_request(request)

        input_genome = self._get_genebe_genome(assembly)

        payload = self._request_conversion(
            representation=request.variant,
            representation_type=request.representation_type,
            input_genome=input_genome,
        )

        raw_candidates = self._extract_variants(
            payload=payload,
            expected_genome="hg38",
        )
        unique_candidates, duplicates_removed = (
            self._deduplicate_candidates(raw_candidates)
        )

        if not unique_candidates:
            raise GeneBeNoCandidatesError(
                "GeneBe could not convert diagnostic variant "
                f"{request.variant!r} on {assembly}"
            )

        warnings = self._build_conversion_warnings(
            representation_type=request.representation_type,
            candidate_count=len(unique_candidates),
            duplicates_removed=duplicates_removed,
        )

        return [
            VariantCandidate(
                candidate_id=f"candidate_{index}",
                chromosome=raw_candidate["chr"],
                position=raw_candidate["pos"],
                reference=raw_candidate["ref"],
                alternate=raw_candidate["alt"],
                assembly="GRCh38",
                conversion_warnings=warnings,
            )
            for index, raw_candidate in enumerate(
                unique_candidates,
                start=1,
            )
        ]

    @staticmethod
    def _validate_request(
            request: VariantConfirmationRequest,
    ) -> None:
        """
        Validate the converter input.
        """
        if not isinstance(request, VariantConfirmationRequest):
            raise TypeError(
                "request must be a VariantConfirmationRequest"
            )

        if request.representation_type is None:
            raise ValueError(
                "Variant representation has not been parsed before "
                "GeneBe conversion"
            )

        if (
                request.representation_type
                not in VariantConfirmationRequest.SUPPORTED_REPRESENTATION_TYPES
        ):
            raise ValueError(
                "Unsupported variant representation type: "
                f"{request.representation_type}"
            )

    def _get_genebe_genome(self, assembly: str) -> str:
        """
        Translate an SFtool assembly name into the GeneBe API value.
        """
        if assembly not in self.ASSEMBLY_TO_GENEBE:
            raise ValueError(
                f"Unsupported reference genome: {assembly}"
            )

        return self.ASSEMBLY_TO_GENEBE[assembly]

    def _request_conversion(
            self,
            representation: str,
            representation_type: str,
            input_genome: str,
    ) -> Any:
        """
        Submit one representation to the GeneBe public conversion endpoint.
        """
        params = {}

        if (
                representation_type == "genomic"
                and input_genome == "hg19"
        ):
            params["inputGenome"] = "hg19"

        try:
            response = self.session.post(
                self.API_URL,
                params=params,
                json=[representation],
                headers={
                    "Accept": "application/json",
                    "Content-Type": "application/json",
                },
                timeout=self.timeout,
            )

            response.raise_for_status()

        except requests.Timeout as exc:
            raise GeneBeRequestError(
                "GeneBe conversion request timed out for "
                f"{representation!r}"
            ) from exc

        except requests.ConnectionError as exc:
            raise GeneBeRequestError(
                "Could not connect to GeneBe while converting "
                f"{representation!r}"
            ) from exc

        except requests.HTTPError as exc:
            status_code = self._get_status_code(response=exc.response)

            message = (
                "GeneBe conversion returned an HTTP error"
            )

            if status_code is not None:
                message += f" {status_code}"

            message += f" for {representation!r}"

            response_body = self._get_safe_response_body(
                response=exc.response
            )

            if response_body:
                message += f": {response_body}"

            raise GeneBeRequestError(message) from exc

        except requests.RequestException as exc:
            raise GeneBeRequestError(
                "GeneBe conversion request failed for "
                f"{representation!r}: {exc}"
            ) from exc

        try:
            return response.json()
        except (ValueError, requests.JSONDecodeError) as exc:
            raise GeneBeResponseError(
                "GeneBe returned a non-JSON response for "
                f"{representation!r}"
            ) from exc

    def _extract_variants(
            self,
            payload: Any,
            expected_genome: str,
    ) -> List[Dict[str, Any]]:
        """
        Validate and extract candidate variants from the GeneBe response.
        """
        if not isinstance(payload, list):
            raise GeneBeResponseError(
                "GeneBe response must be a list"
            )

        if len(payload) != 1:
            raise GeneBeResponseError(
                "GeneBe response must contain exactly one conversion result, "
                f"got {len(payload)}"
            )

        conversion_result = payload[0]

        if not isinstance(conversion_result, dict):
            raise GeneBeResponseError(
                "GeneBe conversion result must be a dictionary"
            )

        if "variants" not in conversion_result:
            message = self._extract_genebe_error_message(
                conversion_result
            )

            if message:
                raise GeneBeNoCandidatesError(
                    f"GeneBe conversion failed: {message}"
                )

            raise GeneBeResponseError(
                "GeneBe conversion result does not contain 'variants'"
            )

        variants = conversion_result["variants"]

        if not isinstance(variants, list):
            raise GeneBeResponseError(
                "GeneBe 'variants' field must be a list"
            )

        if not variants:
            raise GeneBeNoCandidatesError(
                "GeneBe returned no genomic candidates"
            )

        return [
            self._validate_raw_variant(
                raw_variant=raw_variant,
                expected_genome=expected_genome,
                candidate_index=index,
            )
            for index, raw_variant in enumerate(
                variants,
                start=1,
            )
        ]

    def _validate_raw_variant(
            self,
            raw_variant: Any,
            expected_genome: str,
            candidate_index: int,
    ) -> Dict[str, Any]:
        """
        Validate one genomic candidate returned by GeneBe.
        """
        if not isinstance(raw_variant, dict):
            raise GeneBeResponseError(
                "GeneBe candidate "
                f"{candidate_index} must be a dictionary"
            )

        required_fields = {
            "genome",
            "chr",
            "pos",
            "ref",
            "alt",
        }

        missing_fields = required_fields - set(raw_variant.keys())

        if missing_fields:
            missing = ", ".join(sorted(missing_fields))

            raise GeneBeResponseError(
                "GeneBe candidate "
                f"{candidate_index} is missing required fields: {missing}"
            )

        genome = self._validate_non_empty_string(
            raw_variant["genome"],
            field_name="genome",
            candidate_index=candidate_index,
        ).lower()

        if genome not in self.GENEBE_TO_ASSEMBLY:
            raise GeneBeResponseError(
                "GeneBe candidate "
                f"{candidate_index} returned unsupported genome "
                f"{genome!r}"
            )

        if genome != expected_genome:
            requested_assembly = self.GENEBE_TO_ASSEMBLY[
                expected_genome
            ]

            returned_assembly = self.GENEBE_TO_ASSEMBLY[
                genome
            ]

            raise GeneBeAssemblyMismatchError(
                "GeneBe returned genome "
                f"{genome} ({returned_assembly}) for a "
                f"{requested_assembly} conversion request"
            )

        chromosome = self._validate_non_empty_string(
            raw_variant["chr"],
            field_name="chr",
            candidate_index=candidate_index,
        )

        position = raw_variant["pos"]

        if (
                isinstance(position, bool)
                or not isinstance(position, int)
                or position < 1
        ):
            raise GeneBeResponseError(
                "GeneBe candidate "
                f"{candidate_index} field 'pos' must be a positive integer"
            )

        reference = self._validate_non_empty_string(
            raw_variant["ref"],
            field_name="ref",
            candidate_index=candidate_index,
        ).upper()

        alternate = self._validate_non_empty_string(
            raw_variant["alt"],
            field_name="alt",
            candidate_index=candidate_index,
        ).upper()

        if reference == ".":
            raise GeneBeResponseError(
                "GeneBe candidate "
                f"{candidate_index} field 'ref' contains an invalid allele"
            )

        if alternate == ".":
            raise GeneBeResponseError(
                "GeneBe candidate "
                f"{candidate_index} field 'alt' contains an invalid allele"
            )

        return {
            "genome": genome,
            "chr": chromosome,
            "pos": position,
            "ref": reference,
            "alt": alternate,
        }

    @staticmethod
    def _validate_non_empty_string(
            value: Any,
            field_name: str,
            candidate_index: int,
    ) -> str:
        """
        Validate a required string field from the GeneBe response.
        """
        if not isinstance(value, str):
            raise GeneBeResponseError(
                "GeneBe candidate "
                f"{candidate_index} field {field_name!r} must be a string"
            )

        value = value.strip()

        if not value:
            raise GeneBeResponseError(
                "GeneBe candidate "
                f"{candidate_index} field {field_name!r} "
                "must be a non-empty string"
            )

        return value

    @staticmethod
    def _deduplicate_candidates(
            raw_candidates: List[Dict[str, Any]],
    ) -> Tuple[List[Dict[str, Any]], bool]:
        """
        Remove duplicate genomic candidates while preserving result order.
        """
        unique_candidates: List[Dict[str, Any]] = []
        seen = set()
        duplicates_removed = False

        for candidate in raw_candidates:
            key = (
                candidate["genome"],
                candidate["chr"],
                candidate["pos"],
                candidate["ref"],
                candidate["alt"],
            )

            if key in seen:
                duplicates_removed = True
                continue

            seen.add(key)
            unique_candidates.append(candidate)

        return unique_candidates, duplicates_removed

    @staticmethod
    def _build_conversion_warnings(
            representation_type: str,
            candidate_count: int,
            duplicates_removed: bool,
    ) -> List[str]:
        """
        Build warnings associated with conversion ambiguity or deduplication.
        """
        warnings: List[str] = []

        if candidate_count > 1:
            if representation_type == "hgvsp":
                warnings.append(
                    "Protein representation resolved to multiple "
                    "genomic candidates"
                )
            else:
                warnings.append(
                    "Variant representation resolved to multiple "
                    "genomic candidates"
                )

        if duplicates_removed:
            warnings.append(
                "GeneBe returned duplicate genomic candidates; "
                "duplicates were removed"
            )

        return warnings

    @staticmethod
    def _extract_genebe_error_message(
            conversion_result: Dict[str, Any],
    ) -> Optional[str]:
        """
        Extract an error message from a GeneBe conversion result when present.

        The public service response format may evolve, so several common field
        names are accepted.
        """
        for field_name in (
                "error",
                "message",
                "errorMessage",
                "error_message",
        ):
            value = conversion_result.get(field_name)

            if isinstance(value, str) and value.strip():
                return value.strip()

        return None

    @staticmethod
    def _get_status_code(
            response: Optional[requests.Response],
    ) -> Optional[int]:
        """
        Safely obtain the HTTP status code from a response.
        """
        if response is None:
            return None

        return getattr(response, "status_code", None)

    @staticmethod
    def _get_safe_response_body(
            response: Optional[requests.Response],
            maximum_length: int = 300,
    ) -> Optional[str]:
        """
        Return a truncated single-line response body for error reporting.
        """
        if response is None:
            return None

        body = getattr(response, "text", None)

        if not isinstance(body, str):
            return None

        body = " ".join(body.split())

        if not body:
            return None

        if len(body) > maximum_length:
            return body[:maximum_length] + "..."

        return body
