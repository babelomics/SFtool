import json
from modules.variant_selection.utils import index_by_gene_symbol, is_autosomal_chromosome, combine_variant_and_gene_info,\
    classify_genotype, classify_pathogenic_STRs


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


def rr_str_selection(STR_collection, RR_mode, sample_sex):
    '''
    Selection of STRs for Reproductive Risk category according to the following rules
    1. RR_mode: screenin
       Females: get those STRs in HET above pathogenic score in FXN (AR), AFF2 (XLR), DMD (XLR) and ARX (XLR) genes. Get those STRs in HET
       or HOM above intermediate score in FMR1 gene (XLD)
       Males: get those STRs in HET above pathogenic score in FXN (AR). STRs in genes AFF2, DMD, ARX and FMR1 (chrX) are not reported in Males
    2. RR_mode: advanced
        Males and Females: get those STRs in HET or HOM above pathogenic score in FXN (AR), AFF2 (XLR), DMD (XLR) and ARX (XLR) genes. Get those STRs in HET or HOM
        above intermediate score in FMR1 gene (XLD)
    :param STR_collection: Set of STRs selected from current sample with at least one allele above pathogenic/intermediate threshold
    :param RR_mode: mode of reproductive risk. Either screening or advanced
    :param sample_sex: sample sex (female or male)
    :return:
    '''

    STR_selected = {}

    AR_GENES = {"FXN", "AFF2", "DMD", "ARX"}
    XLD_GENES = {"FMR1"}

    for STR_key, current_STR in STR_collection.items():
        gene = current_STR["Gene"]
        repeats = current_STR["Genotype"]

        # -----------------------------
        # Screening mode
        # -----------------------------
        if RR_mode == "screening":
            if gene == "FXN" or (gene in AR_GENES and sample_sex == "female"):
                threshold = current_STR["PathogenicThreshold"]
                if classify_pathogenic_STRs(repeats, threshold) == "HET":
                    STR_selected[STR_key] = current_STR

            elif gene in XLD_GENES and sample_sex == "female":
                threshold = current_STR["IntermediateThreshold"]
                if classify_pathogenic_STRs(repeats, threshold) == "HET":
                    STR_selected[STR_key] = current_STR

        # -----------------------------
        # Advanced mode
        # -----------------------------
        elif RR_mode == "advanced":
            if gene in AR_GENES:
                threshold = current_STR["PathogenicThreshold"]

            elif gene in XLD_GENES:
                threshold = current_STR["IntermediateThreshold"]

            else:
                continue

            if classify_pathogenic_STRs(repeats, threshold) in {"HET", "HOM"}:
                STR_selected[STR_key] = current_STR

    return STR_selected


def rr_smn1_copy_selection(SMN1_copy_collection, RR_mode):
    '''
    Function that returns SMA CARRIER, Inconclusive or Silent carrier individuals (only screening mode, since SMAca only detects carriers)
    :param SMN1_copy_collection
    :param RR_mode: either screening or advanced
    :return: a copy of SMN1_copy_collection
    '''

    SMN1_copy_collection = {}
    if SMN1_copy_collection["call"] in ["LIKELY_SMA_CARRIER (1-copy SMN1)", "Inconclusive", "PUTATIVE_SILENT_SMA_CARRIER"] and RR_mode == 'screening':
        return SMN1_copy_collection
