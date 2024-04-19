

from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_cyp2c9_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Assing a diplotype to CYP2C9 according to variants found in this gene and recommendations based on SEFF.
    Essential alleles for SEFF: *2, *3

    Args:
        variants (list): A list with a dictionary of variants present in the set of pharmacogenetic genes
        diplo_pheno_dct (dict): Dictionary with information about diplotype-phenotype relationships
        aggregated_results (list): A list with results of diplotype information accumulated for each gene

    Returns:
        list: List with results of diplotype information accumulateed for each gene, including CUP2C9

    Raises:
        TypeError: Yes if 'variants' isnot a list, 'diplo_pheno_dict' is not a list and 'aggregated_results' is not a list
    """

    # Check variable type
    if not isinstance(variants, list):
        raise TypeError("'variants' argument should be a list of genetic variants.")
    if not isinstance(diplo_pheno_dct, dict):
        raise TypeError("'diplo_pheno_dct' argument should be a dictionary with information about diplotypes and phenotypes.")
    if not isinstance(aggregated_results, list):
        raise TypeError("'results' argument should bea list where diplotype assignment is stored")

    # Gene
    gene = 'CYP2C9'

    # Check and get variants for 'gene'
    found, variants_gene = check_gene_variants(variants, gene)

    # Diplotype *1/*1 (reference) if there are no variants on current gene
    if not found:
        diplotype = '*1/*1'
    else:
        # Check *2
        if 'rs1799853' in variants_gene.keys():
            if variants_gene['rs1799853'] == '1/1':
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2'
                else: # More than a variant
                    if 'rs1057910' in variants_gene.keys():
                        diplotype = '*2/*3'
                    else:
                        diplotype = 'NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        # Check *3
        elif 'rs1057910' in variants_gene.keys():
            if variants_gene['rs1057910'] == '1/1':
                diplotype = '*3/*3'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3'
                else:
                    diplotype = 'NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        else:
            diplotype = 'NA'
            print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')

    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results