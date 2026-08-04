import os

from sftool.report.models import ReportTable
from sftool.utils.runtime import get_runtime_versions
from sftool.report.variant_confirmation_table import build_variant_confirmation_table


def build_sample_tables(ctx, sample):
    selected_variants = sample.variant_selection

    tables = [
        _build_versions_and_paths_table(ctx, sample),
        build_variant_confirmation_table(ctx, sample),
        _build_pr_rr_snv_indels_table(selected_variants, "PR"),
        _build_pr_rr_snv_indels_table(selected_variants, "RR"),
        _build_rr_str_table(selected_variants),
        _build_rr_smn1_table(selected_variants),
        _build_pgx_table(selected_variants, sample),
    ]

    return [table for table in tables if table is not None]

# ===============================
# Table builders
# ===============================

# Versions and paths
def _build_versions_and_paths_table(ctx, sample):

    category_string = ''
    for category in sample.categories:
        if category == 'PR':
            category_string += 'PR (Personal Risk), '
        elif category == 'RR':
            category_string += 'RR (Reproductive Risk), '
        elif category == 'PGx':
            category_string += 'PGx (Pharmacogenetic Risk), '

    category_string = category_string.rstrip(", ")
    variant_classification_sources_string = ", ".join(ctx.variant_classification_sources)

    versions = get_runtime_versions(ctx.config.paths)

    uses_clinvar = (
            "clinvar"
            in ctx.variant_classification_sources
    )

    rows = [
        {"Field": "SF tool version", "Value": ctx.config.version},
        {"Field": "SF tool general mode", "Value": ", ".join(ctx.modes)},
        {"Field": "Categories", "Value": category_string},
        {"Field": "Reproductive Risk mode", "Value": ctx.RR_mode},
        {"Field": "SF tool variant classification sources", "Value": variant_classification_sources_string},
        {"Field": "Sample ID", "Value": sample.sample_id},
        {"Field": "Sample sex", "Value": sample.sex},
        {"Field": "Sample role", "Value": sample.role},
        {"Field": "HPO list", "Value": ",".join(sample.hpo_terms)},
        {"Field": "Input VCF file", "Value": str(sample.vcf)},
        {"Field": "SMAca file", "Value": _optional_input_path(sample.smaca_path)},
        {"Field": "STRipy file", "Value": _optional_input_path(sample.stripy_path)},
        {"Field": "Personal Risk catalogue file", "Value": _catalog_provenance(ctx, sample, "PR")},
        {"Field": "Reproductive Risk catalogue file", "Value": _catalog_provenance(ctx, sample, "RR")},
        {"Field": "Base output dir", "Value": ctx.base_output_dir},
        {"Field": "Run dir", "Value": ctx.run_dir},
        {"Field": "Temporal dir", "Value": ctx.tmp_dir},
        {"Field": "Human assembly", "Value": "hg19" if ctx.assembly == "GRCh37" else "hg38" },
        {
            "Field": "Reference genome source",
            "Value": (
                ctx.resources.reference_genome_source
            ),
        },
        {
            "Field": "Reference genome",
            "Value": (
                ctx.resources.resource_relative_path(
                    ctx.resources.reference_genome
                )
            ),
        },
        {
            "Field": "ClinVar version",
            "Value": (
                ctx.resources.clinvar_version
                if uses_clinvar
                else "Not used"
            ),
        },
        {"Field": "Clinvar evidence level", "Value": str(ctx.clinvar_evidence) if uses_clinvar else "Not used"},
        {
            "Field": "ClinVar PR resource",
            "Value": _clinvar_provenance(
                ctx,
                sample,
                "PR",
                uses_clinvar,
            ),
        },
        {
            "Field": "ClinVar RR resource",
            "Value": _clinvar_provenance(
                ctx,
                sample,
                "RR",
                uses_clinvar,
            ),
        },
        {
            "Field": "GeneBe version",
            "Value": _genebe_version(
                sample,
                versions,
            ),
        },
        {"Field": "GeneBe path", "Value": ctx.config.paths.genebe},
        {"Field": "bcftools version", "Value": versions["bcftools"].version},
        {"Field": "pharmCAT version", "Value": versions["pharmcat"].version if 'PGx' in sample.categories else "Not used"},
        {
            "Field": "PharmCAT resource version",
            "Value": (
                ctx.resources.pharmcat_version
                if "PGx" in sample.categories
                else "Not used"
            ),
        },
        {
            "Field": "HPO genes-to-phenotype version",
            "Value": ctx.resources.hpo_version,
        },
        {
            "Field": "HPO genes-to-phenotype resource",
            "Value": (
                ctx.resources.resource_relative_path(
                    ctx.resources.hpo_file
                )
            ),
        },
        {
            "Field": "Resource manifest",
            "Value": str(ctx.resources.manifest_path),
        },
        {
            "Field": "Resource manifest schema",
            "Value": str(ctx.resources.manifest["schema_version"])
        },
        {
            "Field": "Resource version policy",
            "Value": (
                ctx.resources.manifest.get(
                    "resource_version_policy",
                    "Not specified",
                )
            ),
        },
    ]

    return ReportTable("Versions and paths", rows)

def _catalog_provenance(
        ctx,
        sample,
        category,
):
    if category not in sample.categories:
        return "Not used"

    path = ctx.resources.catalog_json(category)
    version = ctx.resources.catalog_version(
        category
    )

    relative_path = (
        ctx.resources.resource_relative_path(path)
    )

    if version:
        return f"{version} ({relative_path})"

    return relative_path


def _clinvar_provenance(
        ctx,
        sample,
        category,
        uses_clinvar,
):
    if (
            not uses_clinvar
            or category not in sample.categories
    ):
        return "Not used"

    return ctx.resources.resource_relative_path(
        ctx.resources.clinvar_file(category)
    )


def _optional_input_path(value):
    if not value:
        return "Not provided"

    return str(value)


def _genebe_version(sample, versions):
    if not {"PR", "RR"} & set(sample.categories):
        return "Not used"

    return versions["genebe"].version
# PR and RR SNVs & Indels
def _build_pr_rr_snv_indels_table(variant_selection, category):

    pr_rr_data = variant_selection.get(category)
    snv_indels_data = pr_rr_data.get("snv_indels_genebe_clinvar")
    if not snv_indels_data:
        return None

    rows = []

    for variant_id, gene_entries in snv_indels_data.items():
        # gene_entries could be a list for a variant overlapping more than a gene
        # Normalize to list

        if isinstance(gene_entries, dict):
            gene_entries = [gene_entries]
        for entry in gene_entries:
            row = {
                "Variant": variant_id,
                **entry
            }
            rows.append(row)

    return ReportTable(
        tab_name=category + " results",
        rows=rows
    )


# RR STRs
def _build_rr_str_table(variant_selection):
    rr_data = variant_selection["RR"]
    str_data = rr_data.get("STRs")

    if not str_data:
        return None

    rows = list(str_data.values())
    return ReportTable("RR-STRs", rows)

    ## RR SMA
def _build_rr_smn1_table(variant_selection):
    rr_data = variant_selection["RR"]
    smn1_data = rr_data.get("SMN1_copy")

    if not smn1_data:
        return None

    rows = [smn1_data]
    return ReportTable("RR_SMN1-copy", rows)

# PGx
def _build_pgx_table(variant_selection, sample):
    pgx_data = variant_selection["PGx"]
    pharmCAT_data = pgx_data.get("pharmCAT_variants")

    if not pharmCAT_data:
        return None

    rows = []

    for gene, phenotype_entries in pharmCAT_data.items():

        if isinstance(phenotype_entries, dict):
            phenotype_entries = [phenotype_entries]

        for entry in phenotype_entries:
            genotype = entry.get("genotype", "")
            phenotype = entry.get("phenotype", "")
            rows.append({
                "Gene": gene,
                "Genotype": genotype,
                "Phenotype": phenotype,
            })

    return ReportTable(
        tab_name = "PGx",
        rows=rows,
        metadata={
            "footer": f"PharmCAT full results are available at: {sample.reports['PGx']}",
            "footer_merge": True
        })


