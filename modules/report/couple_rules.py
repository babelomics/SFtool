
class ScreeningRRCoupleRule:
    def build_rows(self, sample_a, sample_b):
        rr_a = sample_a.variant_selection.get("RR", {})
        rr_b = sample_b.variant_selection.get("RR", {})

        # SNVs and Indels

        variants_a = rr_a.get("snv_indels_genebe_clinvar", {})
        variants_b = rr_b.get("snv_indels_genebe_clinvar", {})

        genes_a = self._index_by_gene(sample_a.sample_id, variants_a)
        genes_b = self._index_by_gene(sample_b.sample_id, variants_b)

        shared_genes = sorted(set(genes_a) & set(genes_b))


        # STRs



        # SMA


        rows = []

        for gene in shared_genes:

            for entry_a in genes_a[gene]:
                for entry_b in genes_b[gene]:

                    row = {
                        "Gene": gene,
                        "Sample_1": sample_a.sample_id,
                        "Sample_1_variant": entry_a["Variant"],
                        "Sample_1_genotype": entry_a.get("Genotype"),
                        "Sample_1_hgvsc": entry_a.get("HGVSC"),
                        "Sample_1_hgvsp": entry_a.get("HGVSP"),
                        "Sample_1_clinvar": entry_a.get("ClinvarClinicalSignificance"),
                        "Sample_1_genebe": entry_a.get("GeneBe_ACMG_Classification"),

                        "Sample_2": sample_b.sample_id,
                        "Sample_2_variant": entry_b["Variant"],
                        "Sample_2_genotype": entry_b.get("Genotype"),
                        "Sample_2_hgvsc": entry_b.get("HGVSC"),
                        "Sample_2_hgvsp": entry_b.get("HGVSP"),
                        "Sample_2_clinvar": entry_b.get("ClinvarClinicalSignificance"),
                        "Sample_2_genebe": entry_b.get("GeneBe_ACMG_Classification"),

                        "Inheritance": entry_a.get("inheritance") or entry_b.get("inheritance"),
                        "Phenotype": entry_a.get("Phenotype") or entry_b.get("Phenotype"),
                        "OMIM_disorder": entry_a.get("OMIM_disorder") or entry_b.get("OMIM_disorder"),
                        "Orpha": entry_a.get("Orpha") or entry_b.get("Orpha"),
                        "Couple_risk": "Both partners have reportable variants in the same RR gene"
                    }

                    rows.append(row)

        return rows

    def _index_by_gene(self, sample_id, snv_indels_data):

        index = {}

        for variant_id, gene_entries in snv_indels_data.items():

            if isinstance(gene_entries, dict):
                gene_entries = [gene_entries]

            for entry in gene_entries:

                gene = entry.get("Gene")

                if not gene:
                    continue

                normalized_entry = {
                    "Sample": sample_id,
                    "Variant": variant_id,
                    **entry
                }

                index.setdefault(gene, []).append(normalized_entry)

        return index


class AdvancedRRCoupleRule:

    def build_rows(self, sample_a, sample_b):

        # Placeholder for advanced RR-specific logic
        return -1