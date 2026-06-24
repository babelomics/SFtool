import json

from modules.report.utils import ReportUtilsMixin
class RRCoupleRule(ReportUtilsMixin):

    X_LINKED_INHERITANCE = {"XL", "XLR", "XLD"}

    def __init__(self, rr_catalog_path):
        self.rr_catalog = self._load_rr_catalog(rr_catalog_path)
        self.rr_catalog_by_gene = self._build_catalog_index(self.rr_catalog)

    def build_tables(self, sample_a, sample_b, rr_mode):

        return {
            "SNV&Indel (1-22 and X)": self._build_snv_indels_rows(sample_a, sample_b, rr_mode),
            "SNV&Indel and STRs (FXN gene)": self._build_fxn_rows(sample_a, sample_b, rr_mode),
            "SMN1-copy": self._build_smn1_rows(sample_a, sample_b, rr_mode),
            "STRs (X-linked genes)": self._build_str_xlinked_rows(sample_a, sample_b, rr_mode)
        }


    ##################################
    # 1. SNVs and Indels (Autosomal and X-linked)
    ##################################
    def _build_snv_indels_rows(self, sample_a, sample_b, rr_mode):
        rows = []

        # Get all genes from variants in both samples
        genes = sorted(
            self._get_genes_from_sample(sample_a, "snv/indel")
            |
            self._get_genes_from_sample(sample_b, "snv/indel")
        )

        for gene in genes:
            if gene not in {"FXN", "SMN1"} and self._catalog_has_inheritance(gene, "AR") and not self._is_xlinked_gene(gene):

                # Get variants in both samples for current gene. If screening mode, ensure variants are in HET. Otherwise, HET or HOM are valid
                gene_variants_a = [
                    variant
                    for variant in self._get_gene_snv_indels_entries(sample_a, gene)
                    if (self._is_het(variant) and rr_mode == 'screening') or rr_mode == 'advanced'
                ]

                gene_variants_b = [
                    variant
                    for variant in self._get_gene_snv_indels_entries(sample_b, gene)
                    if (self._is_het(variant) and rr_mode == 'screening') or rr_mode == 'advanced'
                ]

                if gene_variants_a or gene_variants_b:
                    is_ar_ad_gene = self._catalog_has_inheritance(gene, "AD")

                    # --------------------------------------------------
                    # AR + AD genes
                    # Report HET variants in at least one member (screening mode)
                    # or
                    # Report HET or HOM variants in at least one member (advanced mode)
                    # --------------------------------------------------
                    if is_ar_ad_gene:
                        # Sample A
                        for var_a in gene_variants_a:
                            rows.append({
                                "Case study": var_a["inheritance"] + " gene",
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
                                "Case study": var_b["inheritance"] + " gene",
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
                    # Report variants only if both partners have HET variants
                    # in the same gene (screening mode)
                    # or
                    # Report variants only if both partners have HET/HOM variants
                    # in the same gene (advanced mode)
                    # --------------------------------------------------
                    else:
                        if gene_variants_a and gene_variants_b:
                            # Sample A
                            for var_a in gene_variants_a:
                                rows.append({
                                    "Case study": var_a["inheritance"] + " gene",
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
                                    "Case study": var_b["inheritance"] + " gene",
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
        # Only show variants in X genes for female in screening mode

        female_sample = self._get_gender_sample(sample_a, sample_b, "female")

        if female_sample is None:
            return rows

        variants_female = self._get_snv_indels(female_sample)


        if rr_mode == "screening":
            for variant_id, entries in variants_female.items():
                entries = self._as_list(entries)

                for entry in entries:
                    gene = entry.get("Gene")

                    if gene and self._is_xlinked_gene(gene) and self._is_het(entry):
                        var = self._add_variant_id(variant_id, entry)
                        rows.append({
                            "Case study": "X-linked variant in female partner. " + var["inheritance"] + " gene",
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
        else:
            # rr_mode = 'advanced'
            male_sample = self._get_gender_sample(sample_a, sample_b, "male")

            # Get all genes from variants in both samples
            genes = sorted(
                self._get_genes_from_sample(female_sample, "snv/indel")
                |
                self._get_genes_from_sample(male_sample, "snv/indel")
            )

            for gene in genes:
                if self._is_xlinked_gene(gene):
                    # Get variants in both samples for current gene. In advanced mode, HOM variants in male and HET/HOM in female
                    gene_variants_male = [
                        variant
                        for variant in self._get_gene_snv_indels_entries(male_sample, gene)
                    ]

                    gene_variants_female = [
                        variant
                        for variant in self._get_gene_snv_indels_entries(female_sample, gene)
                    ]


                    # Avoid printing variants that have been reported before
                    visited_variants_male = []
                    visited_variants_female = []
                    # Male Sample
                    for var_male in gene_variants_male:
                        for var_female in gene_variants_female:
                            if self._is_hom(var_male) and var_male not in visited_variants_male:
                                rows.append({
                                    "Case study": "X-linked variant. " + var_male["inheritance"] + " gene",
                                    "Sample": male_sample.sample_id,
                                    "Sample role": male_sample.role,
                                    "Sample sex": male_sample.sex,
                                    "Sample partner": female_sample.sample_id,
                                    "Partner role": female_sample.role,
                                    "Partner sex": female_sample.sex,
                                    "Partner has variant in the same gene": "-",
                                    **self._flatten_entry(var_male),
                                    "Warning": ""
                                })
                                visited_variants_male.append(var_male["Variant"])

                            if var_female not in visited_variants_female: # HET or HOM variants in female
                                rows.append({
                                    "Case study": "X-linked variant. " + var_female["inheritance"] + " gene",
                                    "Sample": female_sample.sample_id,
                                    "Sample role": female_sample.role,
                                    "Sample sex": female_sample.sex,
                                    "Sample partner": male_sample.sample_id,
                                    "Partner role": male_sample.role,
                                    "Partner sex": male_sample.sex,
                                    "Partner has variant in the same gene": "-",
                                    **self._flatten_entry(var_female),
                                    "Warning": ""
                                })
                                visited_variants_female.append(var_female["Variant"])


        return rows

    # --------------------------------------------------
    # 2. FXN (SNVs/Indels and STRs)
    # --------------------------------------------------

    def _build_fxn_rows(self, sample_a, sample_b, rr_mode):
        rows = []

        fxn_snv_a = self._get_gene_snv_indels_entries(sample_a, "FXN")
        fxn_snv_b = self._get_gene_snv_indels_entries(sample_b, "FXN")

        fxn_str_a = self._get_gene_str_entries(sample_a, "FXN")
        fxn_str_b = self._get_gene_str_entries(sample_b, "FXN")

        # Avoid printing variants that have been reported before
        visited_variants_a = []
        visited_variants_b = []

        # HET SNV/indels in both paretns
        for var_a in fxn_snv_a:
            for var_b in fxn_snv_b:

                meet_criteria = False

                if self._is_het(var_a) and self._is_het(var_b): # Any mode (screening or advanced)
                    if var_a["Variant"] == var_b["Variant"]:
                        rule = "FXN same heterozygous SNV/indel in both parents. " + var_a["inheritance"] + " gene"
                    else:
                        rule = "FXN compound heterozygous SNV/indel variants in both parents. " + var_a["inheritance"] + " gene"

                    meet_criteria = True

                if ((self._is_hom(var_a) and self._is_het(var_b)) or (self._is_het(var_a) and self._is_hom(var_b)) and rr_mode == 'advanced'):
                    if var_a["Variant"] == var_b["Variant"]:
                        rule = "FXN same SNV/indel in both parents (HOM and HET). " + var_a["inheritance"] + " gene"
                    else:
                        rule = "FXN different SNV/indel variants in both parents (HOM and HET). " + var_a["inheritance"] + " gene"

                    meet_criteria = True

                if self._is_hom(var_a) and self._is_hom(var_b) and rr_mode == 'advanced':
                    if var_a["Variant"] == var_b["Variant"]:
                        rule = "FXN same HOM SNV/indel in both parents. " + var_a["inheritance"] + " gene"
                    else:
                        rule = "FXN different HOM SNV/indel variants in both parents. " + var_a["inheritance"] + " gene"

                    meet_criteria = True

                if meet_criteria:
                    if var_a["Variant"] not in visited_variants_a:
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
                        visited_variants_a.append(var_a["Variant"])
                    if var_b["Variant"] not in visited_variants_b: # Avoid printing variants that have been reported before
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
                        visited_variants_b.append(var_b["Variant"])

        # HET STR in both members
        for str_a in fxn_str_a:
            for str_b in fxn_str_b:
                if (str_a['zigosity'] == 'HET' and str_b['zigosity'] == 'HET' and rr_mode == 'screening') or rr_mode == 'advanced':
                    # In advanced mode, any combination (HET/HET, HET/HOM, HOM/HET or HOM/HOM) is allowed
                    if str_a not in visited_variants_a:
                        rows.append({
                            "Case study": "FXN pathogenic STR in both parents. " + str_a["Inheritance"] + " gene",
                            "Sample": sample_a.sample_id,
                            "Sample role": sample_a.role,
                            "Sample sex": sample_a.sex,
                            "Sample partner": sample_b.sample_id,
                            "Partner role": sample_b.role,
                            "Partner sex": sample_b.sex,
                            "Partner has variant in the same gene": "Yes",
                            **self._flatten_entry(str_a)
                        })
                        visited_variants_a.append(str_a["Variant"])
                    if str_b not in visited_variants_b:
                        rows.append({
                            "Case study": "FXN pathogenic STR in both parents. " + str_b["Inheritance"] + " gene",
                            "Sample": sample_b.sample_id,
                            "Sample role": sample_b.role,
                            "Sample sex": sample_b.sex,
                            "Sample partner": sample_a.sample_id,
                            "Partner role": sample_a.role,
                            "Partner sex": sample_a.sex,
                            "Partner has variant in the same gene": "Yes",
                            **self._flatten_entry(str_b)
                        })
                        visited_variants_b.append(str_b["Variant"])

        # One member with HET SNV/indel and the other member with HET STR
        for var_a in fxn_snv_a:
            for str_b in fxn_str_b:
                meet_criteria = False
                if self._is_het(var_a) and str_b['zigosity'] == 'HET':
                    rule = "FXN HET SNV/indel in one parent and HET STR in the other. " + var_a["inheritance"] + " gene"
                    meet_criteria = True
                if self._is_het(var_a) and str_b['zigosity'] == 'HOM':
                    rule = "FXN HET SNV/indel in one parent and HOM STR in the other. " + var_a["inheritance"] + " gene"
                    meet_criteria = True
                if self._is_hom(var_a) and str_b['zigosity'] == 'HET':
                    rule = "FXN HOM SNV/indel in one parent and HET STR in the other. " + var_a["inheritance"] + " gene"
                    meet_criteria = True
                if self._is_hom(var_a) and str_b['zigosity'] == 'HOM':
                    rule = "FXN HOM SNV/indel in one parent and HOM STR in the other. " + var_a["inheritance"] + " gene"
                    meet_criteria = True

                if meet_criteria:
                    if var_a["Variant"] not in visited_variants_a:
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
                        visited_variants_a.append(var_a["Variant"])

                    if str_b not in visited_variants_b:
                        rows.append({
                            "Case study": rule,
                            "Sample": sample_b.sample_id,
                            "Sample role": sample_b.role,
                            "Sample sex": sample_b.sex,
                            "Sample partner": sample_a.sample_id,
                            "Partner role": sample_a.role,
                            "Partner sex": sample_a.sex,
                            "Partner has variant in the same gene": "Yes",
                            **self._flatten_entry(str_b)
                        })
                        visited_variants_b.append(str_b["Variant"])

        # One member with HET STR and the other member with HET
        for str_a in fxn_str_a:
            for var_b in fxn_snv_b:
                meet_criteria = False
                if str_a["zigosity"] == 'HET' and self._is_het(var_b):
                    rule = "FXN HET STR in one parent and HET SNV/Indel in the other. " + var_b["inheritance"] + " gene"
                    meet_criteria = True
                if str_a["zigosity"] == 'HET' and self._is_hom(var_b):
                    rule = "FXN HET STR in one parent and HOM SNV/Indel in the other. " + var_b["inheritance"] + " gene"
                    meet_criteria = True
                if str_a["zigosity"] == 'HOM' and self._is_het(var_b):
                    rule = "FXN HOM STR in one parent and HET SNV/Indel in the other. " + var_b["inheritance"] + " gene"
                    meet_criteria = True
                if str_a["zigosity"] == 'HOM' and self._is_hom(var_b):
                    rule = "FXN HOM STR in one parent and HOM SNV/Indel in the other. " + var_b["inheritance"] + " gene"
                    meet_criteria = True

                if meet_criteria:
                    if str_a not in visited_variants_a:
                        rows.append({
                            "Case study": rule,
                            "Sample": sample_a.sample_id,
                            "Sample role": sample_a.role,
                            "Sample sex": sample_a.sex,
                            "Sample partner": sample_b.sample_id,
                            "Partner role": sample_b.role,
                            "Partner sex": sample_b.sex,
                            "Partner has variant in the same gene": "Yes",
                            **self._flatten_entry(str_a)
                        })
                        visited_variants_a.append(str_a["Variant"])
                    if var_b["Variant"] not in visited_variants_b:
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
                        visited_variants_b.append(var_b["Variant"])

        return rows

    # --------------------------------------------------
    # 3. SMN1-copy
    # --------------------------------------------------

    def _build_smn1_rows(self, sample_a, sample_b, rr_mode):
        rows = []

        if rr_mode == "screening":

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
                        "Case study": "SMN1-copy. AR gene",
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
                    "Case study": "SMN1-copy. AR gene",
                    "Sample": sample_b.sample_id,
                    "Sample role": sample_b.role,
                    "Sample sex": sample_b.sex,
                    "Sample partner": sample_a.sample_id,
                    "Partner role": sample_a.role,
                    "Partner sex": sample_a.sex,
                    **self._flatten_entry(smn1_b),
                    "Warning": warning
                })
        else: # Advanced mode
            rows.append("SMAca software only detects SMA carriers and its not applied in Advanced Mode.")

        return rows

    # --------------------------------------------------
    # 4. STRs in X-linked genes
    # --------------------------------------------------

    def _build_str_xlinked_rows(self, sample_a, sample_b, rr_mode):
        rows = []

        # Only show STRs in X genes for female
        female_sample = self._get_gender_sample(sample_a, sample_b, "female")

        if female_sample is None:
            return rows

        # STRs in X-linked genes
        str_data_female = self._get_str_data(female_sample)

        if rr_mode == "screening": # In screening mode, only STRs in X-linked genes for female are shown
            for item in str_data_female.items():
                variant_id, entry = item

                gene = entry.get("Gene")

                if gene in {'AFF2', 'DMD', 'ARX', 'FMR1'}:
                    rows.append({
                        "Case study": "X-linked STR in female partner. " + entry["Inheritance"] + " gene",
                        "Sample": female_sample.sample_id,
                        "Sample role": female_sample.role,
                        "Sample sex": female_sample.sex,
                        **self._flatten_entry(item),
                        "Warning": ""
                    })

        else: # Advanced mode: STRs in X-linked genes for males and females
            male_sample = self._get_gender_sample(sample_a, sample_b, "male")

            # STRs in X-linked genes
            #str_data_male = self._get_str_data(male_sample)

            # Get all genes from variants in both samples
            genes = sorted(
                self._get_genes_from_sample(female_sample, "str")
                |
                self._get_genes_from_sample(male_sample, "str")
            )

            for gene in genes:
                if gene in {'AFF2', 'DMD', 'ARX', 'FMR1'}:
                    # Avoid printing variants that have been reported before
                    visited_variants_male = []
                    visited_variants_female = []
                    # Male Sample
                    str_male = self._get_gene_str_entries(male_sample, gene)
                    str_female = self._get_gene_str_entries(female_sample, gene)

                    for entry_male in str_male:
                        for entry_female in str_female:
                            if entry_male["zigosity"] == 'HOM' and (entry_female["zigosity"] == 'HET' or entry_female["zigosity"] == 'HOM'):
                                if entry_male["Variant"] not in visited_variants_male:
                                    rows.append({
                                        "Case study": "X-linked STR variant. " + entry_male["Inheritance"] + " gene",
                                        "Sample": male_sample.sample_id,
                                        "Sample role": male_sample.role,
                                        "Sample sex": male_sample.sex,
                                        "Sample partner": female_sample.sample_id,
                                        "Partner role": female_sample.role,
                                        "Partner sex": female_sample.sex,
                                        "Partner has variant in the same gene": "-",
                                        **self._flatten_entry(item_male),
                                        "Warning": ""
                                    })
                                    visited_variants_male.append(entry_male["Variant"])

                                if entry_female["variant"] not in visited_variants_female:
                                    rows.append({
                                        "Case study": "X-linked STR variant. " + entry_female["Inheritance"] + " gene",
                                        "Sample": female_sample.sample_id,
                                        "Sample role": female_sample.role,
                                        "Sample sex": female_sample.sex,
                                        "Sample partner": male_sample.sample_id,
                                        "Partner role": male_sample.role,
                                        "Partner sex": male_sample.sex,
                                        "Partner has variant in the same gene": "-",
                                        **self._flatten_entry(item_female),
                                        "Warning": ""
                                    })
                                    visited_variants_female.append(entry_female["Variant"])



        return rows
