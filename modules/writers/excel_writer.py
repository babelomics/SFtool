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

    # def write_couple_report(self, couple_report):
    #     output_path = self._get_couple_output_path(
    #         couple_report.sample_a_id,
    #         couple_report.sample_b_id
    #     )
    #
    #     with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
    #         for table in couple_report.tables:
    #             df = pd.DataFrame(table.rows)
    #             df.to_excel(writer, sheet_name=table.tab_name, index=False)

    def write_couple_report(self, couple_report):
        output_path = self._get_couple_output_path(
            couple_report.sample_a_id,
            couple_report.sample_b_id
        )

        with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:

            workbook = writer.book

            description_format = workbook.add_format({
                "bold": False,
                "text_wrap": True,
                "valign": "top"
            })

            wrap_format = workbook.add_format({
                "text_wrap": True,
                "valign": "top"
            })

            for table in couple_report.tables:

                df = pd.DataFrame(table.rows)

                df.to_excel(
                    writer,
                    sheet_name=table.tab_name,
                    index=False
                )

                worksheet = writer.sheets[table.tab_name]

                if table.tab_name == "Report information":
                    worksheet.set_column(0, 0, 25, wrap_format)  # Field
                    worksheet.set_column(1, 2, 80, wrap_format)  # Sample 1 and Sample 2
                else:
                    worksheet.set_column(0, len(df.columns) - 1, 25, wrap_format)
                    for row_num in range(1, len(df) + 1):
                        worksheet.set_row(row_num, 25)

                description = table.metadata.get("description")

                if description:

                    # Leave 2 blank rows after dataframe
                    description_row = len(df) + 4

                    worksheet.merge_range(
                        description_row,
                        0,
                        description_row + 2,
                        len(df.columns) - 1,
                        description,
                        description_format
                    )

    # =====================================
    # Path resolution
    # =====================================

    def _get_sample_output_path(self, sample_id):
        outdir = Path(self.ctx.run_dir)
        return outdir / f"{sample_id}_SFtool_report.xlsx"

    def _get_couple_output_path(self, sample_a_id, sample_b_id):
        outdir = Path(self.ctx.run_dir)
        rr_mode = self.ctx.RR_mode
        if rr_mode == "screening":
            return outdir / f"{sample_a_id}_{sample_b_id}_RR_couple_screening_report.xlsx"
        else:
            return outdir / f"{sample_a_id}_{sample_b_id}_RR_couple_advanced_report.xlsx"
