# modules/report/report_manager.py

from sftool.report.models import SampleReport, CoupleReport
from sftool.report.sample_tables import build_sample_tables
from sftool.report.couple_tables import build_couple_tables
from sftool.writers.excel_writer import ExcelWriter


class ReportManager:

    def __init__(self, ctx):
        self.ctx = ctx
        self.writer = ExcelWriter(ctx)

    # ===============================
    # Individual reports
    # ===============================

    def build_sample_report(self, sample) -> SampleReport:
        report = SampleReport(sample.sample_id)

        for table in build_sample_tables(self.ctx, sample):
            report.add_table(table)

        return report

    def write_sample_report(self, sample_report: SampleReport):
        self.writer.write_sample_report(sample_report)

    def build_couple_report(self, sample_a, sample_b) -> CoupleReport:
        report = CoupleReport(
            sample_a.sample_id,
            sample_b.sample_id,
            self.ctx.RR_mode
        )

        for table in build_couple_tables(self.ctx, sample_a, sample_b):
            report.add_table(table)

        return report

    def write_couple_report(self, couple_report: CoupleReport):
         self.writer.write_couple_report(couple_report)

