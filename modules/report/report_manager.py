# modules/report/report_manager.py

from modules.report.sample_report import SampleReport, ReportTable
from modules.report.couple_report import CoupleReport
from modules.writers.excel_writer import ExcelWriter


class ReportManager:

    def __init__(self, ctx):
        self.ctx = ctx
        self.writer = ExcelWriter(ctx) # TO BE DONE: better ExcelWriter(output_dir=ctx.config.output_dir)

    # ===============================
    # Individual reports
    # ===============================

    def build_sample_report(self, sample) -> SampleReport:
        report = SampleReport(sample.sample_id)

        selected_variants = sample.variant_selection

        # PR
        report.add_table(self._build_pr_rr_snv_indels_table(selected_variants, "PR"))

        # RR
        report.add_table(self._build_pr_rr_snv_indels_table(selected_variants, "RR"))

        # RR-STR
        report.add_table(self._build_rr_str_table(selected_variants))

        # RR-SMN1-copy
        report.add_table(self._build_rr_smn1_table(selected_variants))

        # PGx
        report.add_table(self._build_pgx_table(selected_variants))

        return report

    # def write_sample_report(self, sample_report: SampleReport):
    #     self.writer.write_sample_report(sample_report)
    #
    #
    # def write_couple_report(self, couple_report: CoupleReport):
    #     self.writer.write_couple_report(couple_report)

    # ===============================
    # Table builders
    # ===============================

    # PR and RR SNVs & Indels
    def _build_pr_rr_snv_indels_table(self, variant_selection, category):

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
    def _build_rr_str_table(self, variant_selection):
        rr_data = variant_selection["RR"]
        str_data = rr_data.get("STRs")

        if not str_data:
            return None

        rows = list(str_data.values())
        return ReportTable("RR-STRs", rows)

    ## RR SMA
    def _build_rr_smn1_table(self, variant_selection):
        rr_data = variant_selection["RR"]
        smn1_data = rr_data.get("SMN1_copy")

        if not smn1_data:
            return None

        rows = list(smn1_data.values())
        return ReportTable("RR_SMN1-copy", rows)

    # PGx
    def _build_pgx_table(self, variant_selection):
        pgx_data = variant_selection["PGx"]
        pharmCAT_data = pgx_data.get("pharmCAT_variants")

        if not pharmCAT_data:
            return None

        rows = list(pharmCAT_data.values())
        return ReportTable("PGx", rows)

    # # ===============================
    # # Couple reports (RR only)
    # # ===============================
    #
    # def build_couple_report(self, sample_a, sample_b) -> CoupleReport:
    #     rr_mode = self.ctx.config.rr_mode
    #
    #     couple_report = CoupleReport(
    #         sample_a.sample_id,
    #         sample_b.sample_id,
    #         rr_mode
    #     )
    #
    #     rr_a = sample_a.variant_selection.get("RR", {})
    #     rr_b = sample_b.variant_selection.get("RR", {})
    #
    #     if rr_mode == "screening":
    #         table = self._build_screening_rr_couple_table(rr_a, rr_b)
    #     else:
    #         table = self._build_advanced_rr_couple_table(rr_a, rr_b)
    #
    #     couple_report.add_table(table)
    #
    #     return couple_report
    #
    # def _build_screening_rr_couple_table(self, rr_a, rr_b) -> ReportTable:
    #     rows = []
    #     return ReportTable("RR-Couple-Screening", rows)
    #
    # def _build_advanced_rr_couple_table(self, rr_a, rr_b) -> ReportTable:
    #     rows = []
    #     return ReportTable("RR-Couple-Advanced", rows)
