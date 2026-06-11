import json

from modules.report.utils import ReportUtilsMixin
class ScreeningRRCoupleRule(ReportUtilsMixin):

    X_LINKED_INHERITANCE = {"XL", "XLR", "XLD"}

    def __init__(self, rr_catalog_path):
        self.rr_catalog = self._load_rr_catalog(rr_catalog_path)
        self.rr_catalog_by_gene = self._build_catalog_index(self.rr_catalog)

    def build_tables(self, sample_a, sample_b):

        return {
            "SNV&Indel (1-22 and X)": self._build_snv_indels_rows(sample_a, sample_b),
            "SNV&Indel and STRs (FXN gene)": self._build_fxn_rows(sample_a, sample_b),
            "SMN1-copy": self._build_smn1_rows(sample_a, sample_b),
            "STRs (X-linked genes)": self._build_str_xlinked_rows(sample_a, sample_b)
        }


    ##################################
    # 1. SNVs and Indels (Autosomal and X-linked)
    ##################################
    def _build_snv_indels_rows(self, sample_a, sample_b):
        rows = []

        # Get all genes from variants in both samples
        genes = sorted(
            self._get_genes_from_sample(sample_a)
            |
            self._get_genes_from_sample(sample_b)
        )

        for gene in genes:

            if gene not in {"FXN", "SMN1"} and self._catalog_has_inheritance(gene, "AR") and not self._is_xlinked_gene(gene):

                # Get variants in both samples for current gene
                gene_variants_a = [
                    variant
                    for variant in self._get_gene_snv_indels_entries(sample_a, gene)
                    if self._is_het(variant)
                ]

                gene_variants_b = [
                    variant
                    for variant in self._get_gene_snv_indels_entries(sample_b, gene)
                    if self._is_het(variant)
                ]

                if gene_variants_a or gene_variants_b:
                    is_ar_ad_gene = self._catalog_has_inheritance(gene, "AD")

                    # --------------------------------------------------
                    # AR + AD genes
                    # Report HET variants even if present in only one partner
                    # --------------------------------------------------
                    if is_ar_ad_gene:
                        # Sample A
                        for var_a in gene_variants_a:
                            rows.append({
                                "Case study": "Autosomal AR+AD gene",
                                "Sample": sample_a.sample_id,
                                "Sample role": sample_a.role,
                                "Sample sex": sample_a.sex,
                                "Sample partner": sample_b.sample_id,
                                "Partner role": sample_b.role,
                                "Partner sex": sample_b.sex,
                                "Partner has variant in the same gene": (
                                    "Yes" if gene_variants_b else "No"
                                ),
                                **self._flatten_entry(var_a),
                                "Warning": self._build_ad_warning(gene, var_a)
                            })
                        # Sample B
                        for var_b in gene_variants_b:
                            rows.append({
                                "Case study": "Autosomal AR+AD gene",
                                "Sample": sample_b.sample_id,
                                "Sample role": sample_b.role,
                                "Sample sex": sample_b.sex,
                                "Sample partner": sample_a.sample_id,
                                "Partner role": sample_a.role,
                                "Partner sex": sample_a.sex,
                                "Partner has variant in the same gene": (
                                    "Yes" if gene_variants_a else "No"
                                ),
                                **self._flatten_entry(var_b),
                                "Warning": self._build_ad_warning(gene, var_b)
                            })
                    # --------------------------------------------------
                    # AR-only genes
                    # Report only if both partners have HET variants
                    # in the same gene
                    # --------------------------------------------------
                    else:
                        if gene_variants_a and gene_variants_b:
                            # Sample A
                            for var_a in gene_variants_a:
                                rows.append({
                                    "Case study": "Autosomal AR gene",
                                    "Sample": sample_a.sample_id,
                                    "Sample role": sample_a.role,
                                    "Sample sex": sample_a.sex,
                                    "Sample partner": sample_b.sample_id,
                                    "Partner role": sample_b.role,
                                    "Partner sex": sample_b.sex,
                                    "Partner has variant in the same gene": "Yes",
                                    **self._flatten_entry(var_a),
                                    "Warning": ""
                                })
                            # Sample B
                            for var_b in gene_variants_b:
                                rows.append({
                                    "Case study": "Autosomal AR gene",
                                    "Sample": sample_b.sample_id,
                                    "Sample role": sample_b.role,
                                    "Sample sex": sample_b.sex,
                                    "Sample partner": sample_a.sample_id,
                                    "Partner role": sample_a.role,
                                    "Partner sex": sample_a.sex,
                                    "Partner has variant in the same gene": "Yes",
                                    **self._flatten_entry(var_b),
                                    "Warning": ""
                                })


        # X-linked SNV / Indels
        # Only show variants in X genes for female
        female_sample = self._get_female_sample(sample_a, sample_b)

        if female_sample is None:
            return rows

        variants = self._get_snv_indels(female_sample)

        for variant_id, entries in variants.items():
            entries = self._as_list(entries)

            for entry in entries:
                gene = entry.get("Gene")

                if gene and self._is_xlinked_gene(gene) and self._is_het(entry):
                    var = self._add_variant_id(variant_id, entry)
                    rows.append({
                        "Case study": "X-linked variant in female partner",
                        "Sample": female_sample.sample_id,
                        "Sample role": female_sample.role,
                        "Sample sex": female_sample.sex,
                        "Sample partner": "-",
                        "Partner role": "-",
                        "Partner sex": "-",
                        "Partner has variant in the same gene": "-",
                        **self._flatten_entry(var),
                        "Warning": "X-linked variants are only reported for female partner"
                    })




        return rows

    # --------------------------------------------------
    # 2. FXN (SNVs/Indels and STRs)
    # --------------------------------------------------

    def _build_fxn_rows(self, sample_a, sample_b):
        rows = []

        fxn_snv_a = self._get_gene_snv_indels_entries(sample_a, "FXN")
        fxn_snv_b = self._get_gene_snv_indels_entries(sample_b, "FXN")

        fxn_str_a = self._get_gene_str_entries(sample_a, "FXN")
        fxn_str_b = self._get_gene_str_entries(sample_b, "FXN")

        # HET SNV/indels in both paretns
        for var_a in fxn_snv_a:
            for var_b in fxn_snv_b:
                if self._is_het(var_a) and self._is_het(var_b):

                    if var_a["Variant"] == var_b["Variant"]:
                        rule = "FXN same heterozygous SNV/indel in both parents"
                    else:
                        rule = "FXN compound heterozygous SNV/indel variants in both parents"

                    rows.append({
                        "Case study": rule,
                        "Sample": sample_a.sample_id,
                        "Sample role": sample_a.role,
                        "Sample sex": sample_a.sex,
                        "Sample partner": sample_b.sample_id,
                        "Partner role": sample_b.role,
                        "Partner sex": sample_b.sex,
                        "Partner has variant in the same gene": "Yes",
                        **self._flatten_entry(var_a)
                    })
                    rows.append({
                        "Case study": rule,
                        "Sample": sample_b.sample_id,
                        "Sample role": sample_b.role,
                        "Sample sex": sample_b.sex,
                        "Sample partner": sample_a.sample_id,
                        "Partner role": sample_a.role,
                        "Partner sex": sample_a.sex,
                        "Partner has variant in the same gene": "Yes",
                        **self._flatten_entry(var_b)
                    })

        # HET STR in both members
        for str_a in fxn_str_a:
            for str_b in fxn_str_b:
                if str_a['zigosity'] == 'HET' and str_b['zigosity'] == 'HET':
                    rows.append({
                        "Case study": "FXN pathogenic STR in both parents",
                        "Sample": sample_a.sample_id,
                        "Sample role": sample_a.role,
                        "Sample sex": sample_a.sex,
                        "Sample partner": sample_b.sample_id,
                        "Partner role": sample_b.role,
                        "Partner sex": sample_b.sex,
                        "Partner has variant in the same gene": "Yes",
                        **self._flatten_entry(str_a)
                    })
                    rows.append({
                        "Case study": "FXN pathogenic STR in both parents",
                        "Sample": sample_b.sample_id,
                        "Sample role": sample_b.role,
                        "Sample sex": sample_b.sex,
                        "Sample partner": sample_a.sample_id,
                        "Partner role": sample_a.role,
                        "Partner sex": sample_a.sex,
                        "Partner has variant in the same gene": "Yes",
                        **self._flatten_entry(str_b)
                    })

        # One member with HET SNV/indel and the other member with HET STR
        for var_a in fxn_snv_a:
            for str_b in fxn_str_b:
                if self._is_het(var_a) and str_b['zigosity'] == 'HET':
                    rows.append({
                        "Case study": "FXN HET SNV/indel in one parent and HET STR in the other",
                        "Sample": sample_a.sample_id,
                        "Sample role": sample_a.role,
                        "Sample sex": sample_a.sex,
                        "Sample partner": sample_b.sample_id,
                        "Partner role": sample_b.role,
                        "Partner sex": sample_b.sex,
                        "Partner has variant in the same gene": "Yes",
                        **self._flatten_entry(var_a)
                    })
                    rows.append({
                        "Case study": "FXN HET SNV/indel in one parent and HET STR in the other",
                        "Sample": sample_b.sample_id,
                        "Sample role": sample_b.role,
                        "Sample sex": sample_b.sex,
                        "Sample partner": sample_a.sample_id,
                        "Partner role": sample_a.role,
                        "Partner sex": sample_a.sex,
                        "Partner has variant in the same gene": "Yes",
                        **self._flatten_entry(str_b)
                    })

        # One member with HET STR and the other member with HET
        for str_a in fxn_str_a:
            for var_b in fxn_snv_b:
                if str_a["zigosity"] == 'HET' and self._is_het(var_b):
                    rows.append({
                        "Case study": "HET FXN STR in one parent and HET SNV/indel in the other",
                        "Sample": sample_a.sample_id,
                        "Sample role": sample_a.role,
                        "Sample sex": sample_a.sex,
                        "Sample partner": sample_b.sample_id,
                        "Partner role": sample_b.role,
                        "Partner sex": sample_b.sex,
                        "Partner has variant in the same gene": "Yes",
                        **self._flatten_entry(str_a)
                    })
                    rows.append({
                        "Case study": "HET FXN STR in one parent and HET SNV/indel in the other",
                        "Sample": sample_b.sample_id,
                        "Sample role": sample_b.role,
                        "Sample sex": sample_b.sex,
                        "Sample partner": sample_a.sample_id,
                        "Partner role": sample_a.role,
                        "Partner sex": sample_a.sex,
                        "Partner has variant in the same gene": "Yes",
                        **self._flatten_entry(var_b)
                    })

        return rows

    # --------------------------------------------------
    # 3. SMN1-copy
    # --------------------------------------------------

    def _build_smn1_rows(self, sample_a, sample_b):
        rows = []

        smn1_a = self._get_smn1_data(sample_a)
        smn1_b = self._get_smn1_data(sample_b)

        status_a = self._classify_smn1_carrier(smn1_a)
        status_b = self._classify_smn1_carrier(smn1_b)

        reportable_no_warning = {
            "LIKELY_SMA_CARRIER (1-copy SMN1)", "PUTATIVE_SILENT_SMA_CARRIER"
        }

        reportable_warning = {"Inconclusive"}

        warning = ''

        if status_a and status_b:
            if (status_a in reportable_no_warning and status_b in reportable_warning) or (status_a in reportable_warning and status_b in reportable_no_warning):
                warning = 'Confirmation using an orthogonal technique is highly recommended (e.g. MLPA)'

            rows.append({
                    "Case study": "SMN1-copy",
                    "Sample": sample_a.sample_id,
                    "Sample role": sample_a.role,
                    "Sample sex": sample_a.sex,
                    "Sample partner": sample_b.sample_id,
                    "Partner role": sample_b.role,
                    "Partner sex": sample_b.sex,
                    **self._flatten_entry(smn1_a),
                    "Warning": warning
            })

            rows.append({
                "Case study": "SMN1-copy",
                "Sample": sample_b.sample_id,
                "Sample role": sample_b.role,
                "Sample sex": sample_b.sex,
                "Sample partner": sample_a.sample_id,
                "Partner role": sample_a.role,
                "Partner sex": sample_a.sex,
                **self._flatten_entry(smn1_b),
                "Warning": warning
            })

        return rows

    # --------------------------------------------------
    # 4. STRs in X-linked genes
    # --------------------------------------------------

    def _build_str_xlinked_rows(self, sample_a, sample_b):
        rows = []

        # Only show STRs in X genes for female
        female_sample = self._get_female_sample(sample_a, sample_b)

        if female_sample is None:
            return rows

        # STRs in X-linked genes
        str_data = self._get_str_data(female_sample)

        for item in str_data.items():
            variant_id, entry = item

            gene = entry.get("Gene")

            if gene in {'AFF2', 'DMD', 'ARX', 'FMR1'}:
                rows.append({
                    "Case study": "X-linked STR in female partner",
                    "Sample": female_sample.sample_id,
                    "Sample role": female_sample.role,
                    "Sample sex": female_sample.sex,
                    **self._flatten_entry(item),
                    "Warning": ""
                })



        return rows


class AdvancedRRCoupleRule:

    def build_rows(self, sample_a, sample_b):

        # Placeholder for advanced RR-specific logic
        return -1