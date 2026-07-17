from __future__ import annotations

import csv
import gzip
import re

from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Set, Tuple

from sftool.utils.clinvar_utils import map_review_status
from sftool.variant_confirmation.models import VariantMatch


VariantKey = Tuple[str, int, str, str]


class ClinVarLookupError(RuntimeError):
    """
    Raised when Variant Confirmation ClinVar lookup cannot be completed.
    """


class VariantConfirmationClinVarLookup:
    """
    Retrieve ClinVar information for detected Variant Confirmation matches.

    The lookup uses exact normalized genomic coordinates:

        CHROM + POS + REF + ALT

    No gene catalogue or ``clinvar_evidence`` threshold is applied. Every
    matching ClinVar entry is returned regardless of review status or clinical
    significance.
    """

    ALLOWED_VARIANT_TYPES = {
        "deletion",
        "duplication",
        "insertion",
        "indel",
        "single nucleotide variant",
        "microsatellite",
        "variation",
    }

    def lookup(
            self,
            matches: Sequence[VariantMatch],
            clinvar_db: str | Path,
            clinvar_submission_db: str | Path,
    ) -> Dict[str, dict]:
        """
        Return ClinVar annotations indexed by Variant Confirmation candidate ID.

        Only matches detected in the patient VCF are queried. Candidates not
        found in the sample are omitted from the result.
        """
        validated_matches = self._validate_matches(matches)
        clinvar_db = self._validate_file(clinvar_db, "ClinVar variant database")
        clinvar_submission_db = self._validate_file(
            clinvar_submission_db,
            "ClinVar submission summary database",
        )

        candidate_ids_by_key = self._index_detected_matches(
            validated_matches
        )

        if not candidate_ids_by_key:
            return {}

        annotations = self._read_variant_annotations(
            clinvar_db=clinvar_db,
            candidate_ids_by_key=candidate_ids_by_key,
        )

        variation_ids = {
            str(annotation["clinvar_id"])
            for annotation in annotations.values()
            if annotation.get("clinvar_id") not in (None, "")
        }

        summaries = self._read_submission_summaries(
            clinvar_submission_db=clinvar_submission_db,
            variation_ids=variation_ids,
        )

        for annotation in annotations.values():
            clinvar_id = str(
                annotation.get("clinvar_id", "")
            )
            annotation["clinical_significance_summary"] = summaries.get(
                clinvar_id,
                "",
            )

        return annotations

    @staticmethod
    def _validate_matches(
            matches: Sequence[VariantMatch],
    ) -> list[VariantMatch]:
        if isinstance(matches, (str, bytes)):
            raise TypeError(
                "matches must be a sequence of VariantMatch objects"
            )

        if not isinstance(matches, Sequence):
            raise TypeError(
                "matches must be a sequence of VariantMatch objects"
            )

        validated_matches = list(matches)

        for index, variant_match in enumerate(
                validated_matches,
                start=1,
        ):
            if not isinstance(variant_match, VariantMatch):
                raise TypeError(
                    f"Match {index} must be a VariantMatch, "
                    f"got {type(variant_match).__name__}"
                )

        return validated_matches

    @staticmethod
    def _validate_file(
            path: str | Path,
            description: str,
    ) -> Path:
        if not isinstance(path, (str, Path)):
            raise TypeError(
                f"{description} path must be a string or pathlib.Path"
            )

        path = Path(path)

        if not path.is_file():
            raise FileNotFoundError(
                f"{description} not found: {path}"
            )

        return path

    def _index_detected_matches(
            self,
            matches: Iterable[VariantMatch],
    ) -> Dict[VariantKey, str]:
        candidate_ids_by_key: Dict[VariantKey, str] = {}

        for variant_match in matches:
            if not variant_match.found:
                continue

            key = self._build_variant_key(
                chromosome=variant_match.chromosome,
                position=variant_match.position,
                reference=variant_match.reference,
                alternate=variant_match.alternate,
            )

            if key in candidate_ids_by_key:
                raise ClinVarLookupError(
                    "Detected Variant Confirmation matches contain duplicated "
                    f"normalized coordinates: {self._format_key(key)}"
                )

            candidate_ids_by_key[key] = variant_match.candidate_id

        return candidate_ids_by_key

    def _read_variant_annotations(
            self,
            clinvar_db: Path,
            candidate_ids_by_key: Mapping[VariantKey, str],
    ) -> Dict[str, dict]:
        annotations: Dict[str, dict] = {}

        try:
            with clinvar_db.open(
                    "r",
                    encoding="utf-8",
                    newline="",
            ) as handle:
                reader = csv.DictReader(
                    handle,
                    delimiter="\t",
                )

                self._validate_variant_database_header(
                    reader.fieldnames
                )

                for row in reader:
                    annotation = self._parse_variant_row(row)

                    if annotation is None:
                        continue

                    key = annotation.pop("_variant_key")
                    candidate_id = candidate_ids_by_key.get(key)

                    if candidate_id is None:
                        continue

                    if candidate_id in annotations:
                        raise ClinVarLookupError(
                            "More than one ClinVar record was found for "
                            f"candidate {candidate_id!r} at "
                            f"{self._format_key(key)}"
                        )

                    annotations[candidate_id] = annotation

        except ClinVarLookupError:
            raise
        except (OSError, csv.Error, ValueError) as exc:
            raise ClinVarLookupError(
                f"Could not parse ClinVar variant database {clinvar_db}: {exc}"
            ) from exc

        return annotations

    def _parse_variant_row(
            self,
            row: Mapping[str, str],
    ) -> Optional[dict]:
        variant_type = (
                row.get("Type") or ""
        ).strip().lower()

        if variant_type not in self.ALLOWED_VARIANT_TYPES:
            return None

        position = (
                row.get("PositionVCF") or ""
        ).strip()

        if not position or position == "-1":
            return None

        try:
            position_value = int(position)
        except ValueError as exc:
            raise ClinVarLookupError(
                f"Invalid ClinVar PositionVCF value: {position!r}"
            ) from exc

        chromosome = (
                row.get("Chromosome") or ""
        ).strip()
        reference = (
                row.get("ReferenceAlleleVCF") or ""
        ).strip()
        alternate = (
                row.get("AlternateAlleleVCF") or ""
        ).strip()

        if not chromosome or not reference or not alternate:
            return None

        key = self._build_variant_key(
            chromosome=chromosome,
            position=position_value,
            reference=reference,
            alternate=alternate,
        )

        review_status = (
                row.get("ReviewStatus") or ""
        ).strip()

        stars = (
            map_review_status(review_status)
            if review_status
            else 0
        )

        phenotype_ids = (
                row.get("PhenotypeIDS") or ""
        ).strip()

        rs_number = (
                row.get("RS# (dbSNP)") or ""
        ).strip()

        return {
            "_variant_key": key,
            "variant_name": (
                    row.get("Name") or ""
            ).strip(),
            "gene": (
                    row.get("GeneSymbol") or ""
            ).strip(),
            "clinical_significance": (
                    row.get("ClinicalSignificance") or ""
            ).strip(),
            "clin_sig_simple": (
                    row.get("ClinSigSimple") or ""
            ).strip(),
            "rs": self._format_rs(rs_number),
            "review_status": (
                f"({stars}) {review_status}"
                if review_status
                else ""
            ),
            "stars": stars,
            "clinvar_id": (
                    row.get("VariationID") or ""
            ).strip(),
            "phenotype_ids": phenotype_ids,
            "phenotype_list": (
                    row.get("PhenotypeList") or ""
            ).strip(),
            "orpha": self._extract_identifiers(
                phenotype_ids,
                r"Orphanet:(\d+)",
            ),
            "omim": self._extract_identifiers(
                phenotype_ids,
                r"OMIM:\s*([^,|;]+)",
            ),
        }

    def _read_submission_summaries(
            self,
            clinvar_submission_db: Path,
            variation_ids: Set[str],
    ) -> Dict[str, str]:
        if not variation_ids:
            return {}

        counts: Dict[str, Dict[str, int]] = {
            variation_id: {}
            for variation_id in variation_ids
        }

        try:
            with gzip.open(
                    clinvar_submission_db,
                    "rt",
                    encoding="utf-8",
                    newline="",
            ) as handle:
                header = self._find_submission_header(handle)
                reader = csv.DictReader(
                    handle,
                    fieldnames=header,
                    delimiter="\t",
                )

                for row in reader:
                    variation_id = (
                            row.get("VariationID") or ""
                    ).strip()

                    if variation_id not in counts:
                        continue

                    contributes = (
                            row.get("ContributesToAggregateClassification") or ""
                    ).strip().lower()

                    if contributes != "yes":
                        continue

                    significance = (
                            row.get("ClinicalSignificance") or ""
                    ).strip()

                    if not significance:
                        continue

                    current_counts = counts[variation_id]
                    current_counts[significance] = (
                            current_counts.get(significance, 0) + 1
                    )

        except ClinVarLookupError:
            raise
        except (OSError, csv.Error) as exc:
            raise ClinVarLookupError(
                "Could not parse ClinVar submission summary database "
                f"{clinvar_submission_db}: {exc}"
            ) from exc

        return {
            variation_id: "; ".join(
                f"{significance} ({count})"
                for significance, count in significance_counts.items()
            )
            for variation_id, significance_counts in counts.items()
            if significance_counts
        }

    @staticmethod
    def _find_submission_header(
            handle,
    ) -> list[str]:
        for line in handle:
            if not line.startswith("#"):
                continue

            fields = line.lstrip("#").rstrip("\n").split("\t")

            if "VariationID" in fields:
                return fields

        raise ClinVarLookupError(
            "ClinVar submission summary header was not found"
        )

    @staticmethod
    def _validate_variant_database_header(
            fieldnames: Optional[Sequence[str]],
    ) -> None:
        required_fields = {
            "Type",
            "Name",
            "GeneSymbol",
            "ClinicalSignificance",
            "ClinSigSimple",
            "RS# (dbSNP)",
            "VariationID",
            "PhenotypeIDS",
            "PhenotypeList",
            "Chromosome",
            "ReviewStatus",
            "PositionVCF",
            "ReferenceAlleleVCF",
            "AlternateAlleleVCF",
        }

        available_fields = set(fieldnames or [])
        missing_fields = required_fields - available_fields

        if missing_fields:
            raise ClinVarLookupError(
                "ClinVar variant database is missing required columns: "
                + ", ".join(sorted(missing_fields))
            )

    @staticmethod
    def _build_variant_key(
            chromosome: str,
            position: int,
            reference: str,
            alternate: str,
    ) -> VariantKey:
        chromosome = str(chromosome).strip()

        if chromosome.lower().startswith("chr"):
            chromosome = chromosome[3:]

        return (
            chromosome,
            int(position),
            str(reference).strip().upper(),
            str(alternate).strip().upper(),
        )

    @staticmethod
    def _format_key(
            key: VariantKey,
    ) -> str:
        chromosome, position, reference, alternate = key

        return (
            f"{chromosome}:{position}:{reference}:{alternate}"
        )

    @staticmethod
    def _format_rs(
            rs_number: str,
    ) -> str:
        if not rs_number or rs_number in {"-1", "."}:
            return ""

        if rs_number.lower().startswith("rs"):
            return rs_number

        return f"rs{rs_number}"

    @staticmethod
    def _extract_identifiers(
            value: str,
            pattern: str,
    ) -> str:
        return ",".join(
            dict.fromkeys(
                match.strip()
                for match in re.findall(pattern, value)
                if match.strip()
            )
        )
