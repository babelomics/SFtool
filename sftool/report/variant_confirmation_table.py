from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set

from sftool.report.models import ReportTable
from sftool.variant_confirmation.models import (
    VariantCandidate,
    VariantConfirmationResult,
    VariantMatch,
)


TAB_NAME = "Variant confirmation"
NOT_AVAILABLE = "NA"


def build_variant_confirmation_table(ctx, sample) -> Optional[ReportTable]:
    """
    Build the per-sample Variant Confirmation report table.

    One row is generated for every candidate/GeneBe annotation combination.
    Candidates without GeneBe annotations still generate one row so that
    negative and unannotated results remain visible in the report.
    """
    result = sample.results.get("variant_confirmation")

    if result is None:
        return None

    if not isinstance(result, VariantConfirmationResult):
        raise TypeError(
            "sample.results['variant_confirmation'] must be a "
            "VariantConfirmationResult"
        )

    gene_hpos = _load_gene_hpos(
        ctx.config.references.gene_to_phenotype_file
    )
    sample_hpos = set(sample.hpo_terms or [])

    rows: List[Dict[str, Any]] = []

    for candidate in result.candidates:
        variant_match = result.get_match(candidate.candidate_id)

        if variant_match is None:
            raise ValueError(
                "Variant Confirmation result is incomplete: no match found "
                f"for candidate {candidate.candidate_id!r}"
            )

        annotations = variant_match.annotations or [None]

        for annotation in annotations:
            rows.append(
                _build_report_row(
                    result=result,
                    candidate=candidate,
                    variant_match=variant_match,
                    annotation=annotation,
                    gene_hpos=gene_hpos,
                    sample_hpos=sample_hpos,
                )
            )

    if not rows:
        return None

    return ReportTable(
        tab_name=TAB_NAME,
        rows=rows,
    )


def _build_report_row(
        result: VariantConfirmationResult,
        candidate: VariantCandidate,
        variant_match: VariantMatch,
        annotation: Optional[dict],
        gene_hpos: Mapping[str, Set[str]],
        sample_hpos: Set[str],
) -> Dict[str, Any]:
    annotation = annotation or {}
    clinvar = variant_match.clinvar or {}
    gene = annotation.get("gene") or ""

    return {
        "Input_variant": result.request.variant,
        "Input_format": result.request.representation_type or "",
        "Candidate_ID": candidate.candidate_id,
        "Candidate_variant": candidate.get_genomic_variant(),
        "Normalized_variant": candidate.get_normalized_variant() or "",
        "Found_in_sample": "Yes" if variant_match.found else "No",
        "Variant": variant_match.get_variant(),
        "Gene": gene,
        "Genotype": variant_match.genotype or "",
        "Zygosity": _classify_zygosity(
            variant_match.genotype
        ),
        "rs": (
                annotation.get("dbsnp")
                or clinvar.get("rs")
                or ""
        ),
        "Transcript": annotation.get("transcript") or "",
        "HGVSC": annotation.get("hgvsc") or "",
        "HGVSP": annotation.get("hgvsp") or "",
        "Consequence": annotation.get("consequence") or "",
        "GeneBe_ACMG_Classification": (
                annotation.get("acmg_classification") or ""
        ),
        "GeneBe_ACMG_criteria": _join_values(
            annotation.get("acmg_criteria")
        ),
        "ClinvarClinicalSignificance": (
                clinvar.get("clinical_significance") or NOT_AVAILABLE
        ),
        "ClinvarSummary": (
                clinvar.get("clinical_significance_summary") or NOT_AVAILABLE
        ),
        "ReviewStatus": (
                clinvar.get("review_status") or NOT_AVAILABLE
        ),
        "ClinvarID": (
                clinvar.get("clinvar_id") or NOT_AVAILABLE
        ),
        "Orpha": (
                clinvar.get("orpha") or NOT_AVAILABLE
        ),
        "OMIM_clinvar": (
                clinvar.get("omim") or NOT_AVAILABLE
        ),
        "Quality": (
            "" if variant_match.quality is None else variant_match.quality
        ),
        "Filter": _join_values(variant_match.filters),
        "related_HPOs_for_sample": _get_related_hpos(
            gene=gene,
            gene_hpos=gene_hpos,
            sample_hpos=sample_hpos,
        ),
        "VCF_Sample_FORMAT": _format_sample_format(
            variant_match.sample_format
        ),
        "Conversion_warnings": _join_values(
            candidate.conversion_warnings
        ),
        "Match_warnings": _join_values(
            variant_match.warnings
        ),
    }


def _load_gene_hpos(
        gene_to_phenotype_file: str | Path,
) -> Dict[str, Set[str]]:
    """
    Load the same gene-to-phenotype resource used by PR/RR reporting.

    Only direct intersections with the sample HPO list are reported, matching
    the current behaviour of ``add_patient_HPOterms``.
    """
    gene_to_phenotype_file = Path(gene_to_phenotype_file)
    gene_hpos: Dict[str, Set[str]] = {}

    with gene_to_phenotype_file.open(
            "r",
            encoding="utf-8",
    ) as handle:
        next(handle, None)

        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                raise ValueError(
                    "Invalid gene-to-phenotype row at "
                    f"{gene_to_phenotype_file}:{line_number}"
                )

            gene_symbol = fields[1].strip()
            hpo_id = fields[2].strip()

            if not gene_symbol or not hpo_id:
                continue

            gene_hpos.setdefault(
                gene_symbol,
                set(),
            ).add(
                hpo_id
            )

    return gene_hpos


def _get_related_hpos(
        gene: str,
        gene_hpos: Mapping[str, Set[str]],
        sample_hpos: Set[str],
) -> str:
    if not gene or not sample_hpos:
        return NOT_AVAILABLE

    related_hpos = sorted(
        gene_hpos.get(gene, set()) & sample_hpos
    )

    return (
        ",".join(related_hpos)
        if related_hpos
        else NOT_AVAILABLE
    )


def _classify_zygosity(
        genotype: Optional[str],
) -> str:
    if not genotype or genotype in {".", "./.", ".|."}:
        return ""

    separator = (
        "/"
        if "/" in genotype
        else "|"
        if "|" in genotype
        else None
    )

    if separator is None:
        # Haploid alternate genotype, typically on sex chromosomes.
        return "HEMI" if genotype not in {"0", "."} else ""

    alleles = genotype.split(separator)

    if len(alleles) != 2 or "." in alleles:
        return ""

    if alleles[0] == alleles[1]:
        return "HOM"

    return "HET"


def _format_sample_format(
        sample_format: Mapping[str, Any],
) -> str:
    """
    Serialize all patient VCF FORMAT fields into one Excel-safe value.

    Example:
        GT=0/1; AD=46,38; DP=84; GQ=99
    """
    if not sample_format:
        return ""

    return "; ".join(
        f"{key}={_format_value(value)}"
        for key, value in sample_format.items()
    )


def _format_value(value: Any) -> str:
    if value is None:
        return ""

    if isinstance(value, (list, tuple)):
        return ",".join(
            str(item)
            for item in value
        )

    return str(value)


def _join_values(
        values: Any,
) -> str:
    if values is None:
        return ""

    if isinstance(values, str):
        return values

    if isinstance(values, Iterable):
        return ",".join(
            str(value)
            for value in values
        )

    return str(values)
