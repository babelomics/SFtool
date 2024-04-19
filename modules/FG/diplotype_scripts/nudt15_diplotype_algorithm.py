
from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_nudt15_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
     Assing a diplotype to NUDT15 according to variants found in this gene and recommendations based on SEFF
     Essential alleles for SEFF: *3

     Args:
         variants (list): A list with a dictionary of variants present in the set of pharmacogenetic genes
         diplo_pheno_dct (dict): Dictionary with information about diplotype-phenotype relationships
         aggregated_results (list): A list with results of diplotype information accumulated for each gene

     Returns:
         list: List with results of diplotype information accumulateed for each gene, including CUP2C9


     """
    # Gene
    gene = 'NUDT15'

    # Check and get variants for 'gene'
    found, variants_gene = check_gene_variants(variants, gene)

    # Diplotype *1/*1 (reference) if there are no variants on current gene
    if not found:
        diplotype = '*1/*1'
    else:
        # Check *3
        if 'rs116855232' in variants_gene.keys():
            if len(variants_gene) == 1:
                if variants_gene['rs116855232'] == '1/1':
                    diplotype = '*3/*3'
                else:
                    diplotype = '*1/*3'
        else:
            diplotype='NA'
            print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')


    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results
