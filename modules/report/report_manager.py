# modules/report/report_manager.py

from modules.report.sample_report import SampleReport, ReportTable
from modules.report.couple_report import CoupleReport
from modules.report.couple_rules import RRCoupleRule
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
    def write_couple_report(self, couple_report: CoupleReport):
         self.writer.write_couple_report(couple_report)

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
    # Screening Couple report description (RR only)
    # ===============================

    def _build_screening_rr_couple_description(self, sample_1, sample_2):

        description = (
            "This excel report summarizes reproductive risk findings identified in the analyzed couple (screening mode). Results are organized into " +
            "1) SNVs/Indels in autosomal and X chromosomes tab. Contains: Variants in HET in both parents (same vriant or compound heterogizosity) or " +
            " variants in HET in a single parent for those genes with AR and AD inheritance mode (GJB2, CHRNE, ABCC8, AIRE and ALPL). " +
            "Variants in X chromosomes are only reported for females. " +
            "2) SNVs/Indels and STRs in FXN gene. Contains variants in HET in both parents. " +
            "3) SMN1-copy. Results from SMAca software are showed for both parents (1-copy carrier / silent carrier). " +
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
                "description": description,
                "description_merge": True
            }
        )

    # ===============================
    # Screening Couple report description (RR only)
    # ===============================

    def _build_advanced_rr_couple_description(self, sample_1, sample_2):

        description = (
                "This excel report summarizes reproductive risk findings identified in the analyzed couple (advanced mode). Results are organized into " +
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
                "description": description,
                "description_merge": True
            }
        )


    # ===============================
    # Couple reports (RR only)
    # ===============================

    def build_couple_report(self, sample_a, sample_b) -> CoupleReport:
        rr_mode = self.ctx.RR_mode

        couple_report = CoupleReport(
            sample_a.sample_id,
            sample_b.sample_id,
            rr_mode
        )

        tables = self._build_rr_couple_tables(
            sample_a,
            sample_b,
            rr_mode
        )

        if rr_mode == "screening":
            # Table corresponding to description of screening rr
            couple_report.add_table(self._build_screening_rr_couple_description(sample_a, sample_b))
        else:
            couple_report.add_table(self._build_advanced_rr_couple_description(sample_a, sample_b))

        for table in tables:
            couple_report.add_table(table)

        return couple_report


    def _build_rr_couple_tables(self, sample_a, sample_b, rr_mode) -> ReportTable:

        table_rows = RRCoupleRule(self.ctx.outputs['catalogs']['json_files']['RR']).build_tables(sample_a, sample_b, rr_mode)
        tables = []

        for tab_name, rows in table_rows.items():
            if rows:
                tables.append(
                    ReportTable(
                        tab_name=tab_name,
                        rows=rows
                    )
                )

        return tables


