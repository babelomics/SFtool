from sftool.report.models import ReportTable
from sftool.report.couple_rules import RRCoupleRule

def build_couple_tables(ctx, sample_a, sample_b):
    rr_mode = ctx.RR_mode

    tables = []

    if rr_mode == "screening":
        tables.append(
            _build_screening_rr_couple_description(sample_a, sample_b, rr_mode)
        )
    else:
        tables.append(
            _build_advanced_rr_couple_description(sample_a, sample_b, rr_mode)
        )

    rr_catalog_path = ctx.resources.catalog_json("RR")

    table_rows = RRCoupleRule(rr_catalog_path).build_tables(
        sample_a,
        sample_b,
        rr_mode
    )

    for tab_name, rows in table_rows.items():
        if rows:
            tables.append(
                ReportTable(
                    tab_name=tab_name,
                    rows=rows
                )
            )

    return [table for table in tables if table is not None]


# ===============================
# Screening Couple report description (RR only)
# ===============================

def _build_screening_rr_couple_description(sample_1, sample_2, rr_mode):

    description = (
            "This Excel report summarizes reproductive risk findings identified in the analyzed couple (screening mode). Results are organized into " +
            "1) SNVs/Indels in autosomal and X chromosomes tab. Contains: Variants in HET in both parents (same variant or compound heterogyzosity) or " +
            " variants in HET in a single parent for those genes with AR and AD inheritance mode (GJB2, CHRNE, ABCC8, AIRE and ALPL). " +
            "Variants in X chromosomes are only reported for females. " +
            "2) SNVs/Indels and STRs in FXN gene. Contains variants in HET in both parents. " +
            "3) SMN1-copy. Results from SMAca software are shown for both parents (1-copy carrier / silent carrier). " +
            "4) STRs. Variants in HET are shown only for females"
    )


    rows = [
        {
            "Field": "Sample ID",
            "Sample 1": sample_1.sample_id,
            "Sample 2": sample_2.sample_id
        },
        {
            "Field": "Sample sex",
            "Sample 1": sample_1.sex,
            "Sample 2": sample_2.sex
        },
        {
            "Field": "Sample role",
            "Sample 1": sample_1.role,
            "Sample 2": sample_2.role
        },
        {
            "Field": "HPO list",
            "Sample 1": ",".join(sample_1.hpo_terms),
            "Sample 2": ",".join(sample_2.hpo_terms)
        },
        {
            "Field": "Input VCF file",
            "Sample 1": str(sample_1.vcf),
            "Sample 2": str(sample_2.vcf)
        },
        {
            "Field": "SMAca file",
            "Sample 1": "Not provided" if sample_1.smaca_path == '' else sample_1.smaca_path,
            "Sample 2": "Not provided" if sample_2.smaca_path == '' else sample_2.smaca_path
        },
        {
            "Field": "STRipy file",
            "Sample 1": "Not provided" if sample_1.stripy_path == '' else sample_1.stripy_path,
            "Sample 2": "Not provided" if sample_2.stripy_path == '' else sample_2.stripy_path
        }
    ]

    return ReportTable(
        tab_name="Report information",
        rows=rows,
        metadata={
            "rr_mode": rr_mode,
            "description": description,
            "description_merge": True
        }
    )

# ===============================
# Advanced Couple report description (RR only)
# ===============================

def _build_advanced_rr_couple_description(sample_1, sample_2, rr_mode):

    description = (
            "This Excel report summarizes reproductive risk findings identified in the analyzed couple (advanced mode). Results are organized into " +
            "1) SNVs/Indels in autosomal and X chromosomes tab. Contains: HET/HOM variants in both parents for the same gene or " +
            " HET/HOM variants in a single parent for those genes with AR and AD inheritance mode (GJB2, CHRNE, ABCC8, AIRE and ALPL). " +
            "For X chromosome, HOM variants in male and HET/HOM variants in female for the same gene are reported." +
            "2) SNVs/Indels and STRs in FXN gene. Contains HET/HOM SNV/Indels/STRs in both parents. " +
            "3) STRs. HOM variants in male and HET/HOM variants in female for the same gene are reported. " +
            "No results are shown for SMN1-copy alterations since SMAca is a carrier screening tool not suitable for SMA diagnosis."
    )


    rows = [
        {
            "Field": "Sample ID",
            "Sample 1": sample_1.sample_id,
            "Sample 2": sample_2.sample_id
        },
        {
            "Field": "Sample sex",
            "Sample 1": sample_1.sex,
            "Sample 2": sample_2.sex
        },
        {
            "Field": "Sample role",
            "Sample 1": sample_1.role,
            "Sample 2": sample_2.role
        },
        {
            "Field": "HPO list",
            "Sample 1": ",".join(sample_1.hpo_terms),
            "Sample 2": ",".join(sample_2.hpo_terms)
        },
        {
            "Field": "Input VCF file",
            "Sample 1": str(sample_1.vcf),
            "Sample 2": str(sample_2.vcf)
        },
        {
            "Field": "STRipy file",
            "Sample 1": "Not provided" if sample_1.stripy_path == '' else sample_1.stripy_path,
            "Sample 2": "Not provided" if sample_2.stripy_path == '' else sample_2.stripy_path
        }

    ]

    return ReportTable(
        tab_name="Report information",
        rows=rows,
        metadata={
            "rr_mode": rr_mode,
            "description": description,
            "description_merge": True
        }
    )
