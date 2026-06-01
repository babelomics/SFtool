# modules/writers/excel_writer.py

import pandas as pd
from pathlib import Path


class ExcelWriter:

    def __init__(self, ctx):
        self.ctx = ctx

    # =====================================
    # Sample report writing
    # =====================================

    def write_sample_report(self, sample_report):
        output_path = self._get_sample_output_path(sample_report.sample_id)

        with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
            workbook = writer.book
            bold = workbook.add_format({"bold": True})

            for table in sample_report.tables:
                if table.tab_name == "Versions and paths":
                    df = pd.DataFrame(table.rows)
                    df.to_excel(
                        writer,
                        sheet_name=table.tab_name,
                        index=False,
                        header=False
                    )
                    worksheet = writer.sheets[table.tab_name]

                    # Make first column bold
                    worksheet.set_column(0, 0, 40, bold)

                    # Optional widths
                    worksheet.set_column(1, 1, 150)

                else:
                    df = pd.DataFrame(table.rows)
                    df.to_excel(writer, sheet_name=table.tab_name, index=False)

    # =====================================
    # Couple report writing
    # =====================================

    def write_couple_report(self, couple_report):
        output_path = self._get_couple_output_path(
            couple_report.sample_a_id,
            couple_report.sample_b_id
        )

        with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
            for table in couple_report.tables:
                df = pd.DataFrame(table.rows)
                df.to_excel(writer, sheet_name=table.name, index=False)

    # =====================================
    # Path resolution
    # =====================================

    def _get_sample_output_path(self, sample_id):
        outdir = Path(self.ctx.run_dir)
        return outdir / f"{sample_id}_SFtool_report.xlsx"

    def _get_couple_output_path(self, sample_a_id, sample_b_id):
        outdir = Path(self.ctx.run_dir)
        return outdir / f"{sample_a_id}_{sample_b_id}_RR_couple_report.xlsx"
