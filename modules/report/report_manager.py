# modules/report/report_manager.py

from modules.report.sample_report import SampleReport, ReportTable
from modules.report.couple_report import CoupleReport
from modules.report.couple_rules import (ScreeningRRCoupleRule, AdvancedRRCoupleRule)
from modules.writers.excel_writer import ExcelWriter
import os
import subprocess
import re


class ReportManager:

    def __init__(self, ctx):
        self.ctx = ctx
        self.writer = ExcelWriter(ctx)

    # ===============================
    # Individual reports
    # ===============================

    def build_sample_report(self, sample) -> SampleReport:
        report = SampleReport(sample.sample_id)


        # Versions and paths
        report.add_table(self._build_versions_and_paths_table(sample))

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

    def write_sample_report(self, sample_report: SampleReport):
        self.writer.write_sample_report(sample_report)
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

        rows = [smn1_data]
        return ReportTable("RR_SMN1-copy", rows)

    # PGx
    def _build_pgx_table(self, variant_selection):
        pgx_data = variant_selection["PGx"]
        pharmCAT_data = pgx_data.get("pharmCAT_variants")

        if not pharmCAT_data:
            return None

        rows = list(pharmCAT_data.values())
        return ReportTable("PGx", rows)

    # Versions and paths
    def _build_versions_and_paths_table(self, sample):

        try:
            cmd = [self.ctx.config.paths.java,
                   "-jar",
                   self.ctx.config.paths.genebe,
                   "version"
                   ]

            # Run command and get output
            genebe_process = subprocess.Popen(cmd, stdout = subprocess.PIPE)
            genebe_out, genebe_err = genebe_process.communicate()

        except subprocess.CalledProcessError as e:
            print(f"Error running GeneBe: {e.output}")


        # Bcftools version
        try:
            cmd = [self.ctx.config.paths.bcftools, "--version"]

            bcftools_process = subprocess.Popen(cmd, stdout= subprocess.PIPE)
            bcftools_out, bcftools_err = bcftools_process.communicate()

        except subprocess.CalledProcessError as e:
            print(f"Error running bcftools: {e.output}")

        # PharmCAT version
        try:
            pharmCAT_command = [self.ctx.config.paths.java, "-jar", self.ctx.config.paths.pharmCAT , "-version"]
            pharmCAT_process = subprocess.Popen(pharmCAT_command, stdout = subprocess.PIPE)
            pharmCAT_output, pharmCAT_err = pharmCAT_process.communicate()

        except subprocess.CalledProcessError as e:
            print(f"Error running pharmCAT: {e.output}")

        category_string = ''
        for category in sample.categories:
            if category == 'PR':
                category_string += 'PR (Personal Risk), '
            elif category == 'RR':
                category_string += 'RR (Reproductive Risk), '
            elif category == 'PGx':
                category_string += 'PGx (Pharmacogenetic Risk), '

        category_string = category_string.rstrip(", ")


        rows = [
            {"Field": "SF tool version", "Value": self.ctx.config.version},
            {"Field": "SF tool general mode", "Value": self.ctx.mode},
            {"Field": "Categories", "Value": category_string},
            {"Field": "Reproductive Risk mode", "Value": self.ctx.RR_mode},
            {"Field": "SF tool pathogenicity profile", "Value": self.ctx.profile},
            {"Field": "Sample ID", "Value": sample.sample_id},
            {"Field": "Sample sex", "Value": sample.sex},
            {"Field": "Sample role", "Value": sample.role},
            {"Field": "HPO list", "Value": ",".join(sample.hpo_terms)},
            {"Field": "Input VCF file", "Value": str(sample.vcf)},
            {"Field": "SMAca file", "Value": "Not provided" if sample.smaca_path == '' else sample.smaca_path},
            {"Field": "STRipy file", "Value": "Not provided" if sample.stripy_path == '' else sample.stripy_path},
            {"Field": "Personal Risk catalogue file", "Value": self.ctx.config.catalogs.personal_risk_geneset if 'PR' in sample.categories else "Not used"},
            {"Field": "Reproductive Risk catalogue file", "Value": self.ctx.config.catalogs.reproductive_risk_geneset if 'RR' in sample.categories else "Not used"},
            {"Field": "Base output dir", "Value": self.ctx.base_output_dir},
            {"Field": "Run dir", "Value": self.ctx.run_dir},
            {"Field": "Temporal dir", "Value": self.ctx.tmp_dir},
            {"Field": "Human assembly", "Value": "hg19" if self.ctx.assembly == "GRCh37" else "hg38" },
            {"Field": "Reference genome path", "Value": self.ctx.config.references.genomes["GRCh37"] if self.ctx.assembly == "GRCh37" else self.ctx.config.reference.genomes["GRCh38"]},
            {"Field": "Clinvar version", "Value": self.ctx.config.clinvar.version if self.ctx.profile == "advanced" else "Not used"},
            {"Field": "Clinvar path", "Value": self.ctx.config.clinvar.db_path if self.ctx.profile == "advanced" else "Not used"},
            {"Field": "Clinvar evidence level", "Value": str(self.ctx.clinvar_evidence) if self.ctx.profile == "advanced" else "Not used"},
            {"Field": "GeneBe version", "Value": "Not used" if ("PR" not in sample.categories and "rr" not in sample.categories) else re.search(r'version:\s*(.*?)\s*::', str(genebe_out)).group(1)},
            {"Field": "GeneBe path", "Value": self.ctx.config.paths.genebe},
            {"Field": "bcftools version", "Value": str(bcftools_out).split(" ")[1].split("\\n")[0]},
            {"Field": "pharmCAT version", "Value": pharmCAT_output.decode().strip() if 'PGx' in sample.categories else "Not used"},
            {"Field": "HPO genes to phenotype version", "Value": os.path.splitext(os.path.basename(self.ctx.config.references.gene_to_phenotype_file))[0].split("_")[-1]},
            {"Field": "HPO genes to phenotype path", "Value": self.ctx.config.references.gene_to_phenotype_file}
        ]

        return ReportTable("Versions and paths", rows)


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

        if rr_mode == "screening":
            table = self._build_screening_rr_couple_table(sample_a, sample_b)
        else:
            table = self._build_advanced_rr_couple_table(sample_a, sample_b)

        couple_report.add_table(table)

        return couple_report

    def _build_screening_rr_couple_table(self, sample_a, sample_b) -> ReportTable:
        rows = ScreeningRRCoupleRule().build_rows(sample_a, sample_b)

        if not rows:
            return None

        return ReportTable("RR Couple Screening", rows)

    def _build_advanced_rr_couple_table(self, sample_a, sample_b) -> ReportTable:
        rows = AdvancedRRCoupleRule().build_rows(sample_a, sample_b)

        if not rows:
            return None
        return ReportTable("RR Couple Advanced", rows)
