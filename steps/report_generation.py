# steps/report_generation.py

from modules.report.report_manager import ReportManager


def run(ctx):
    """
    Pipeline step: report_generation

    Responsibilities:
    - Build one SampleReport per sample
    - Write one Excel per sample
    - If two samples are present, build and write one RR-only couple report
    """

    manager = ReportManager(ctx)

    # Build & write individual reports
    for sample in ctx.samples:
        sample_report = manager.build_sample_report(sample)
        manager.write_sample_report(sample_report)

    # Couple report (RR only)
    if len(ctx.samples) == 2:
        couple_report = manager.build_couple_report(
            ctx.samples[0],
            ctx.samples[1]
        )
        manager.write_couple_report(couple_report)
