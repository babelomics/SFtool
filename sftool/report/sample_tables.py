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
        {"Field": "SMAca file", "Value": "Not provided" if sample.smaca_path == '' else sample.smaca_path},
        {"Field": "STRipy file", "Value": "Not provided" if sample.stripy_path == '' else sample.stripy_path},
        {"Field": "Personal Risk catalogue file", "Value": ctx.config.catalogs.personal_risk_geneset if 'PR' in sample.categories else "Not used"},
        {"Field": "Reproductive Risk catalogue file", "Value": ctx.config.catalogs.reproductive_risk_geneset if 'RR' in sample.categories else "Not used"},
        {"Field": "Base output dir", "Value": ctx.base_output_dir},
        {"Field": "Run dir", "Value": ctx.run_dir},
        {"Field": "Temporal dir", "Value": ctx.tmp_dir},
        {"Field": "Human assembly", "Value": "hg19" if ctx.assembly == "GRCh37" else "hg38" },
        {"Field": "Reference genome path", "Value": ctx.config.references.genomes["GRCh37"] if ctx.assembly == "GRCh37" else ctx.config.references.genomes["GRCh38"]},
        {"Field": "Clinvar version", "Value": ctx.config.clinvar.version if 'clinvar' in ctx.variant_classification_sources else "Not used"},
        {"Field": "Clinvar path", "Value": ctx.config.clinvar.db_path if 'clinvar' in ctx.variant_classification_sources else "Not used"},
        {"Field": "Clinvar evidence level", "Value": str(ctx.clinvar_evidence) if 'clinvar' in ctx.variant_classification_sources else "Not used"},
        {"Field": "GeneBe version", "Value": "Not used" if ("PR" not in sample.categories and "rr" not in sample.categories) else versions["genebe"].version},
        {"Field": "GeneBe path", "Value": ctx.config.paths.genebe},
        {"Field": "bcftools version", "Value": versions["bcftools"].version},
        {"Field": "pharmCAT version", "Value": versions["pharmcat"].version if 'PGx' in sample.categories else "Not used"},
        {"Field": "HPO genes to phenotype version", "Value": os.path.splitext(os.path.basename(ctx.config.references.gene_to_phenotype_file))[0].split("_")[-1]},
        {"Field": "HPO genes to phenotype path", "Value": ctx.config.references.gene_to_phenotype_file}
    ]

    return ReportTable("Versions and paths", rows)

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


