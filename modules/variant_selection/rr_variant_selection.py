import json
from modules.variant_selection.utils import index_by_gene_symbol, is_autosomal_chromosome, combine_variant_and_gene_info, classify_genotype


def rr_variant_selection(snv_indels_pr_collection, rr_json_file, RR_mode, sample_sex):


    # Load JSON file for the given category. This file contains the inheritance mode for each gene
    genes_cat = None
    with open(rr_json_file, "r") as genes_cat_file:
        genes_cat = json.load(genes_cat_file)
        gene_cat_index = index_by_gene_symbol(genes_cat['genes'])

    # Create a dictionary with the set of variants to be informed
    snv_indels_selected = {}

    # Move through the set of collected SNV and Indels
    for variant_key, variant_info_list in snv_indels_pr_collection.items():
        for variant_info in variant_info_list:
            variant_gene = variant_info["Gene"]
            if variant_gene in gene_cat_index:
                gene = gene_cat_index[variant_gene]
                chr = variant_key.split(':')[0]
                if ((RR_mode == 'screening' and classify_genotype(variant_info["Genotype"]) == 'HET' and \
                     (is_autosomal_chromosome(chr) or ((chr == 'X' or chr == 'chrX') and sample_sex == 'female' and gene != 'FMR1'))) or \
                        RR_mode == 'advanced'):
                    # Merge information for gene and variant
                    combined_info = combine_variant_and_gene_info(variant_info, gene)
                    # Add merged information to the dictionary
                    snv_indels_selected[variant_key] = combined_info



    return snv_indels_selected