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

            wrap_format = workbook.add_format({
                "text_wrap": True,
                "valign": "top"
            })

            for table in sample_report.tables:
                df = pd.DataFrame(table.rows)
                if table.tab_name == "Versions and paths":
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
                    worksheet.set_column(1, 1, 150, wrap_format)

                elif table.tab_name == "PGx":

                    df.to_excel(
                        writer,
                        sheet_name=table.tab_name,
                        index=False
                    )

                    worksheet = writer.sheets[table.tab_name]

                    worksheet.set_column(0, 0, 20, wrap_format)   # Gene
                    worksheet.set_column(1, 1, 40, wrap_format)   # Genotype
                    worksheet.set_column(2, 2, 40, wrap_format)   # Phenotype
                    worksheet.set_column(3, 3, 20, wrap_format)   # Source


                    footer = table.metadata.get("footer")

                    if footer:
                        footer_row = len(df) + 2    # one empty row after the table
                        worksheet.merge_range(
                            footer_row,
                            0,
                            footer_row,
                            3,
                            footer,
                            wrap_format
                        )

                else:
                    df.to_excel(writer, sheet_name=table.tab_name, index=False)
                    worksheet = writer.sheets[table.tab_name]
                    self._autosize_columns(
                        worksheet,
                        df,
                        wrap_format,
                        min_width=15,
                        max_width=60
                    )

    # =====================================
    # Couple report writing
    # =====================================


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

            rr_mode_format = workbook.add_format({
                "bold": True,
                "text_wrap": True,
                "valign": "top"
            })

            for table in couple_report.tables:

                df = pd.DataFrame(table.rows)

                if table.tab_name == "SMN1-copy" and "Gene" not in df.columns:
                    df["Gene"] = "SMN1"

                # Put Gene as second column when present
                if "Gene" in df.columns:
                    columns = list(df.columns)
                    columns.remove("Gene")
                    columns.insert(1, "Gene")
                    df = df[columns]

                rr_mode = table.metadata.get("rr_mode")

                if table.tab_name == "Report information" and rr_mode:
                    startrow = 2
                else:
                    startrow = 0

                df.to_excel(
                    writer,
                    sheet_name=table.tab_name,
                    index=False,
                    startrow=startrow
                )

                worksheet = writer.sheets[table.tab_name]

                if table.tab_name == "Report information" and rr_mode:
                    worksheet.merge_range(0, 0, 0, len(df.columns) - 1, f"Reproductive Risk mode: {rr_mode}", rr_mode_format)

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


    def _autosize_columns(self, worksheet, df, cell_format,
                          min_width=15, max_width=60):
        for i, column in enumerate(df.columns):

            max_len = max(
                len(str(column)),
                *(len(str(v)) for v in df[column].fillna(""))
            )

            width = min(max(max_len + 2, min_width), max_width)

            worksheet.set_column(i, i, width, cell_format)