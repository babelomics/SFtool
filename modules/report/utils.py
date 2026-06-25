import json

class ReportUtilsMixin:

    def _load_rr_catalog(self, rr_catalog_path):
        with open(rr_catalog_path, "r", encoding="utf-8") as f:
            return json.load(f)

    def _build_catalog_index(self, rr_catalog):
        return {
            gene_entry["gene_symbol"]: gene_entry
            for gene_entry in rr_catalog.get("genes", [])
        }

    def _get_snv_indels(self, sample):
        return (
            sample.variant_selection
            .get("RR", {})
            .get("snv_indels_genebe_clinvar", {})
        )

    def _get_smn1_data(self, sample):
        return (
            sample.variant_selection
            .get("RR", {})
            .get("SMN1_copy")
        )

    def _get_str_data(self, sample):
        return (
            sample.variant_selection
            .get("RR", {})
            .get("STRs", {})
        )

    def _get_gene_snv_indels_entries(self, sample, gene):
        output = []

        for variant_id, entries in self._get_snv_indels(sample).items():
            entries = self._as_list(entries)

            for entry in entries:
                if entry.get("Gene") == gene:
                    output.append(self._add_variant_id(variant_id, entry))

        return output

    def _get_gene_str_entries(self, sample, gene):
        output = []

        for str_id, entry in self._get_str_data(sample).items():
            if entry.get("Gene") == gene:
                output.append(self._add_variant_id(str_id, entry))

        return output

    def _add_variant_id(self, variant_id, entry):
        return {
            "Variant": variant_id,
            **entry
        }

    def _as_list(self, value):
        if isinstance(value, list):
            return value

        if isinstance(value, dict):
            return [value]

        return []

    def _get_catalog_field(self, gene, field):
        return self.rr_catalog_by_gene.get(gene, {}).get(field, "")

    def _get_catalog_inheritance(self, gene):
        return self._get_catalog_field(gene, "inheritance")

    def _catalog_has_inheritance(self, gene, inheritance):
        inheritance_value = self._get_catalog_inheritance(gene)

        tokens = (
            inheritance_value
            .replace(";", ",")
            .replace("/", ",")
            .split(",")
        )

        tokens = {token.strip().upper() for token in tokens if token.strip()}

        return inheritance.upper() in tokens

    def _clinvar_has_inheritance(self, variant_entry, inheritance):
        OMIM_disorders_catalog = variant_entry["OMIM_disorder"] # OMIM terms from corresponding catalog
        inheritance_catalog = variant_entry["inheritance"] # Inheritance mode from corresponding catalog
        OMIM_disorders_clinvar = variant_entry["OMIM_clinvar"] # OMIM terms for the gene in clinvar

        # Get lists of OMIM and their corresponding inheritance mode from the catalog
        catalog_omims = [
            omim.strip()
            for omim in str(OMIM_disorders_catalog).split(";")
        ]

        catalog_inheritances = [
            inh.strip().upper()
            for inh in str(inheritance_catalog).split(";")
        ]

        # Get only OMIM terms from the corresponding catalog whose inheritance is equal to desired value
        target_omims = {
            omim
            for omim, inh in zip(catalog_omims, catalog_inheritances)
            if inh == inheritance.upper()
        }

        # Get matches with CLinvar OMIM terms
        clinvar_omims = {
            omim.strip()
            for omim in str(OMIM_disorders_clinvar).split(",")
            if omim.strip()
        }

        return sorted(target_omims & clinvar_omims)



    def _is_xlinked_gene(self, gene):
        inheritance_value = self._get_catalog_inheritance(gene)

        tokens = (
            inheritance_value
            .replace(";", ",")
            .replace("/", ",")
            .split(",")
        )

        tokens = {token.strip().upper() for token in tokens if token.strip()}

        return bool(tokens & self.X_LINKED_INHERITANCE)

    def _is_het(self, variant):
        genotype = str(variant.get("Genotype", "")).upper()
        return genotype in {"0/1", "1/0", "HET", "HETEROZYGOUS"}

    def _is_hom(self, variant):
        genotype = str(variant.get("Genotype", "")).upper()
        return genotype in {"1/1", "1/1", "HOM", "HOMOZYGOUS"}


    def _is_pathogenic_str(self, str_entry):
        value = " ".join(str(v).lower() for v in str_entry.values())

        return (
                "pathogenic" in value
                or "full_mutation" in value
                or "expanded" in value
        )

    def _classify_smn1_carrier(self, smn1_data):
        call = smn1_data.get("call", "")
        return call

    def _get_gender_sample(self, sample_a, sample_b, current_gender):
        if str(getattr(sample_a, "sex", "")).lower() == current_gender:
            return sample_a

        if str(getattr(sample_b, "sex", "")).lower() == current_gender:
            return sample_b

        return None

    def _get_genes_from_sample(self, sample, type):
        genes = set()

        if type == "snv/indel":
            variant_set = self._get_snv_indels(sample).items()
        elif type == "str":
            variant_set = self._get_str_data(sample).items()


        for variant_id, entries in variant_set:
            entries = self._as_list(entries)

            for entry in entries:
                gene = entry.get("Gene")
                if gene:
                    genes.add(gene)

        return genes

    def _build_ad_warning(self, gene, variant):
        warnings = []

        if self._catalog_has_inheritance(gene, "AD"):
            warnings.append(
                "WARNING: According to RR catalog, this gene is also associated "
                "with autosomal dominant (AD) inheritance."
            )

        matching_omims = self._clinvar_has_inheritance(variant, "AD")

        if matching_omims:
            warnings.append(
                "WARNING: According to ClinVar, this variant has OMIM disorder(s) "
                "consistent with AD inheritance in the RR catalog: "
                + ", ".join(matching_omims)
                + "."
            )

        return " ".join(warnings)

    def _flatten_entry(self, entry: dict) -> dict:
        """
        Convert a variant/STR/SMA entry dictionary into flat Excel columns.
        Nested dicts/lists are converted to JSON strings.
        """
        # Case 1: tuple like ("chrX:147582158-147582203", {...})
        if (
                isinstance(entry, tuple)
                and len(entry) == 2
                and isinstance(entry[1], dict)
        ):
            variant_id, entry_dict = entry

            entry = {
                "Variant": variant_id,
                **entry_dict
            }

        # Case 2: normal dict
        elif isinstance(entry, dict):
            entry = entry

        # Case 3: unsupported object
        else:
            return {"value": entry}

        flattened = {}

        for key, value in entry.items():
            if isinstance(value, (dict, list)):
                flattened[key] = json.dumps(value, ensure_ascii=False)
            else:
                flattened[key] = value

        return flattened