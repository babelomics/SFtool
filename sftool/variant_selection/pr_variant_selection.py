import json
from sftool.variant_selection.utils import check_specific_criteria, combine_variant_and_gene_info, index_by_gene_symbol, classify_genotype, add_patient_HPOterms

# PR module

def pr_variant_selection(snv_indels_pr_collection, pr_json_file, assembly, gene_to_phenotype_file, sample_hpo_terms):


    # Load JSON file for the given category. This file contains the inheritance mode for each gene
    genes_cat = None
    with open(pr_json_file, "r") as genes_cat_file:
        genes_cat = json.load(genes_cat_file)
        gene_cat_index = index_by_gene_symbol(genes_cat['genes'])

    # Create a dictionary with the set of variants to be informed
    snv_indels_selected = {}



    # Move through the set of combined results
    for variant_key, variant_info_list in snv_indels_pr_collection.items():
        for variant_info in variant_info_list: # A variant key might overlap more than a gene
            variant_gene = variant_info["Gene"]
            # Get inheritance mode for the given gene
            if variant_gene in gene_cat_index:
                gene = gene_cat_index[variant_gene]
                inher = gene["inheritance"]
                # Check for specific consequence or genomic variant
                if check_specific_criteria(gene, variant_info, variant_key, assembly):
                    # For Personal Risk module, report variant if autosomic Dominant (AD), SemiDominant (SD) and X-linked (XL) inheritance mode. For Reproductive Risk category, report variant in any case
                    if inher in ['AD', 'SD', 'XL']:
                        # Merge Variant and gene information
                        combined_info = combine_variant_and_gene_info(variant_info, gene)
                        # Add merged info into the final dataset
                        snv_indels_selected[variant_key] = combined_info

                    # For Autosomic Recessive, check genotype and/or other variants in the same gene
                    elif inher == 'AR':
                        # For a variant in HOM, report variant
                        if classify_genotype(variant_info["Genotype"]) == 'HOM':
                            # Merge information for gene and variant
                            combined_info = combine_variant_and_gene_info(variant_info, gene)
                            # Add merged information to the dictionary
                            snv_indels_selected[variant_key] = combined_info

                        # For a variant in HET, only report if there is another variant in the same gene
                        elif classify_genotype(variant_info["Genotype"]) == 'HET':
                            # Look for another variant in the same gene
                            other_variant_in_gene = False
                            for other_variant_key in snv_indels_pr_collection:
                                other_variant_info_list = snv_indels_pr_collection[other_variant_key]
                                for other_variant_info in other_variant_info_list:
                                    if other_variant_info["Gene"] == variant_gene and other_variant_key != variant_key:
                                        other_combined_info = combine_variant_and_gene_info(other_variant_info, gene)
                                        combined_info = combine_variant_and_gene_info(variant_info, gene)

                                        # Add merged information of the two variants of the gene to the dictionary
                                        snv_indels_selected[variant_key] = combined_info
                                        snv_indels_selected[other_variant_key] = other_combined_info

    # Add HPO terms
    snv_indels_selected_with_HPO = {}
    if snv_indels_selected:
        snv_indels_selected_with_HPO = add_patient_HPOterms(snv_indels_selected, 'snv_indels', gene_to_phenotype_file, sample_hpo_terms)
    return snv_indels_selected_with_HPO


    return snv_indels_selected