from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

from sftool.variant_confirmation.models import VariantCandidate


class VariantLiftoverError(RuntimeError):
    """Base exception for blocking variant liftover failures."""


class GeneBeLiftoverRequestError(VariantLiftoverError):
    """Raised when the GeneBe liftover endpoint cannot be reached."""


class GeneBeLiftoverResponseError(VariantLiftoverError):
    """Raised when GeneBe returns malformed liftover data."""


class GeneBeNoLiftoverCandidatesError(VariantLiftoverError):
    """Raised when a variant cannot be lifted to the target assembly."""


class GeneBeVariantLiftover:
    """
    Lift genomic VariantCandidate objects between genome assemblies using
    the public GeneBe liftover API.
    """

    API_URL = "https://api.genebe.net/cloud/api-public/v1/liftover"

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
        if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
            raise TypeError("timeout must be a positive number")

        if timeout <= 0:
            raise ValueError("timeout must be greater than zero")

        self.session = session or requests.Session()
        self.timeout = float(timeout)

    def lift(
            self,
            candidates: List[VariantCandidate],
            target_assembly: str,
    ) -> List[VariantCandidate]:
        """
        Lift all candidates to the requested assembly.

        Candidates already in the target assembly are returned unchanged.
        """
        if target_assembly not in self.ASSEMBLY_TO_GENEBE:
            raise ValueError(
                f"Unsupported target assembly: {target_assembly}"
            )

        lifted_candidates = []

        for candidate in candidates:
            if candidate.assembly == target_assembly:
                lifted_candidates.append(candidate)
                continue

            lifted_candidates.extend(
                self._lift_candidate(
                    candidate=candidate,
                    target_assembly=target_assembly,
                )
            )

        return lifted_candidates

    def _lift_candidate(
            self,
            candidate: VariantCandidate,
            target_assembly: str,
    ) -> List[VariantCandidate]:
        source_genome = self._get_genebe_genome(candidate.assembly)
        target_genome = self._get_genebe_genome(target_assembly)

        payload = self._request_liftover(
            candidate=candidate,
            source_genome=source_genome,
            target_genome=target_genome,
        )

        variants = self._extract_variants(
            payload=payload,
            expected_source=source_genome,
            expected_destination=target_genome,
        )

        warnings = list(candidate.conversion_warnings)

        if len(variants) > 1:
            warnings.append(
                "Variant liftover resolved to multiple genomic candidates"
            )

        return [
            VariantCandidate(
                candidate_id=(
                    candidate.candidate_id
                    if len(variants) == 1
                    else f"{candidate.candidate_id}_lifted_{index}"
                ),
                chromosome=variant["chr"],
                position=variant["pos"],
                reference=variant["ref"],
                alternate=variant["alt"],
                assembly=target_assembly,
                conversion_warnings=list(warnings),
            )
            for index, variant in enumerate(variants, start=1)
        ]

    def _request_liftover(
            self,
            candidate: VariantCandidate,
            source_genome: str,
            target_genome: str,
    ) -> Any:
        query = self._build_query(candidate)

        try:
            response = self.session.get(
                self.API_URL,
                params={
                    "query": query,
                    "from": source_genome,
                    "dest": target_genome,
                },
                headers={
                    "Accept": "*/*",
                },
                timeout=self.timeout,
            )

            response.raise_for_status()

        except requests.Timeout as exc:
            raise GeneBeLiftoverRequestError(
                f"GeneBe liftover timed out for {query!r}"
            ) from exc

        except requests.RequestException as exc:
            raise GeneBeLiftoverRequestError(
                f"GeneBe liftover failed for {query!r}: {exc}"
            ) from exc

        try:
            return response.json()
        except (ValueError, requests.JSONDecodeError) as exc:
            raise GeneBeLiftoverResponseError(
                f"GeneBe returned a non-JSON liftover response for {query!r}"
            ) from exc

    @staticmethod
    def _build_query(candidate: VariantCandidate) -> str:
        chromosome = candidate.chromosome

        if chromosome.lower().startswith("chr"):
            chromosome = chromosome[3:]

        return (
            f"{chromosome}-"
            f"{candidate.position}-"
            f"{candidate.reference}-"
            f"{candidate.alternate}"
        )

    def _extract_variants(
            self,
            payload: Any,
            expected_source: str,
            expected_destination: str,
    ) -> List[Dict[str, Any]]:
        if not isinstance(payload, dict):
            raise GeneBeLiftoverResponseError(
                "GeneBe liftover response must be a dictionary"
            )

        source = self._validate_assembly_field(
            payload=payload,
            field_name="from",
        )

        destination = self._validate_assembly_field(
            payload=payload,
            field_name="dest",
        )

        if source != expected_source:
            raise GeneBeLiftoverResponseError(
                f"GeneBe liftover returned source genome {source!r}; "
                f"expected {expected_source!r}"
            )

        if destination != expected_destination:
            raise GeneBeLiftoverResponseError(
                f"GeneBe liftover returned destination genome "
                f"{destination!r}; expected {expected_destination!r}"
            )

        variants = payload.get("variants")

        if not isinstance(variants, list):
            raise GeneBeLiftoverResponseError(
                "GeneBe liftover field 'variants' must be a list"
            )

        if not variants:
            raise GeneBeNoLiftoverCandidatesError(
                "GeneBe returned no lifted genomic candidates"
            )

        return [
            self._validate_variant(variant, index)
            for index, variant in enumerate(variants, start=1)
        ]

    @staticmethod
    def _validate_assembly_field(
            payload: Dict[str, Any],
            field_name: str,
    ) -> str:
        value = payload.get(field_name)

        if not isinstance(value, str) or not value.strip():
            raise GeneBeLiftoverResponseError(
                f"GeneBe liftover field {field_name!r} "
                "must be a non-empty string"
            )

        return value.strip().lower()

    @staticmethod
    def _validate_variant(
            variant: Any,
            candidate_index: int,
    ) -> Dict[str, Any]:
        if not isinstance(variant, dict):
            raise GeneBeLiftoverResponseError(
                f"Lifted candidate {candidate_index} must be a dictionary"
            )

        required_fields = {"chr", "pos", "ref", "alt"}
        missing_fields = required_fields - set(variant)

        if missing_fields:
            missing = ", ".join(sorted(missing_fields))
            raise GeneBeLiftoverResponseError(
                f"Lifted candidate {candidate_index} is missing: {missing}"
            )

        chromosome = variant["chr"]
        position = variant["pos"]
        reference = variant["ref"]
        alternate = variant["alt"]

        if not isinstance(chromosome, str) or not chromosome.strip():
            raise GeneBeLiftoverResponseError(
                f"Lifted candidate {candidate_index} has invalid chromosome"
            )

        if (
                isinstance(position, bool)
                or not isinstance(position, int)
                or position < 1
        ):
            raise GeneBeLiftoverResponseError(
                f"Lifted candidate {candidate_index} has invalid position"
            )

        if not isinstance(reference, str) or not reference:
            raise GeneBeLiftoverResponseError(
                f"Lifted candidate {candidate_index} has invalid REF"
            )

        if not isinstance(alternate, str) or not alternate:
            raise GeneBeLiftoverResponseError(
                f"Lifted candidate {candidate_index} has invalid ALT"
            )

        return {
            "chr": chromosome,
            "pos": position,
            "ref": reference.upper(),
            "alt": alternate.upper(),
        }

    def _get_genebe_genome(self, assembly: str) -> str:
        try:
            return self.ASSEMBLY_TO_GENEBE[assembly]
        except KeyError as exc:
            raise ValueError(
                f"Unsupported source assembly: {assembly}"
            ) from exc