# modules/report/report_manager.py

from modules.report.sample_report import SampleReport, ReportTable
from modules.report.couple_report import CoupleReport
from modules.writers.excel_writer import ExcelWriter


class ReportManager:

    def __init__(self, ctx):
        self.ctx = ctx
        self.writer = ExcelWriter(ctx)

    # ===============================
    # Individual reports
    # ===============================

    def build_sample_report(self, sample) -> SampleReport:
        report = SampleReport(sample.sample_id)

        selected = sample.variant_selection

        # PR
        if "PR" in selected:
            table = self._build_pr_table(selected["PR"])
            report.add_table(table)

        # RR
        if "RR" in selected:
            table = self._build_rr_table(selected["RR"])
            report.add_table(table)

        # RR-STR
        if "RR_STR" in selected:
            table = self._build_rr_str_table(selected["RR_STR"])
            report.add_table(table)

        # RR-SMN1
        if "RR_SMN1" in selected:
            table = self._build_rr_smn1_table(selected["RR_SMN1"])
            report.add_table(table)

        # PGx
        if "PGx" in selected:
            table = self._build_pgx_table(selected["PGx"])
            report.add_table(table)

        return report

    def write_sample_report(self, sample_report: SampleReport):
        self.writer.write_sample_report(sample_report)

    # ===============================
    # Couple reports (RR only)
    # ===============================

    def build_couple_report(self, sample_a, sample_b) -> CoupleReport:
        rr_mode = self.ctx.config.rr_mode

        couple_report = CoupleReport(
            sample_a.sample_id,
            sample_b.sample_id,
            rr_mode
        )

        rr_a = sample_a.variant_selection.get("RR", {})
        rr_b = sample_b.variant_selection.get("RR", {})

        if rr_mode == "screening":
            table = self._build_screening_rr_couple_table(rr_a, rr_b)
        else:
            table = self._build_advanced_rr_couple_table(rr_a, rr_b)

        couple_report.add_table(table)

        return couple_report

    def write_couple_report(self, couple_report: CoupleReport):
        self.writer.write_couple_report(couple_report)

    # ===============================
    # Table builders (pure logic)
    # ===============================

    def _build_pr_table(self, pr_data) -> ReportTable:
        rows = []  # TODO: normalize PR rows
        return ReportTable("PR", rows)

    def _build_rr_table(self, rr_data) -> ReportTable:
        rows = []
        return ReportTable("RR", rows)

    def _build_rr_str_table(self, rr_str_data) -> ReportTable:
        rows = []
        return ReportTable("RR-STR", rows)

    def _build_rr_smn1_table(self, rr_smn1_data) -> ReportTable:
        rows = []
        return ReportTable("RR-SMN1-copy", rows)

    def _build_pgx_table(self, pgx_data) -> ReportTable:
        rows = []
        return ReportTable("PGx", rows)

    def _build_screening_rr_couple_table(self, rr_a, rr_b) -> ReportTable:
        rows = []
        return ReportTable("RR-Couple-Screening", rows)

    def _build_advanced_rr_couple_table(self, rr_a, rr_b) -> ReportTable:
        rows = []
        return ReportTable("RR-Couple-Advanced", rows)
