
from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_cyp2c19_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Assing a diplotype to CYP2C19 according to variants found in this gene and recommendations based on SEFF
    Essential alleles for SEFF: *2, *3, *17

    Args:
        variants (list): A list with a dictionary of variants present in the set of pharmacogenetic genes
        diplo_pheno_dct (dict): Dictionary with information about diplotype-phenotype relationships
        aggregated_results (list): A list with results of diplotype information accumulated for each gene

    Returns:
        list: List with results of diplotype information accumulateed for each gene, including CUP2C9


    """
    # Gene
    gene = 'CYP2C19'

    # Take into account the following: Variant rs3758581 is considered as reference allele. It matches to CYP2C19 reference sequenceLRG_584 / NG_008384.3.

    # Check and get variants for 'gene'
    found, variants_gene = check_gene_variants(variants, gene)

    # Diplotype *1/*1 (reference) if there are no variants on current gene or there is the variant related to rs3758581
    if not found or (len(variants_gene) == 1 and 'rs3758581' in variants_gene.keys()):
        diplotype = '*1/*1'
    else:
        # Check *2 allele
        if 'rs12769205' in variants_gene.keys() and 'rs4244285' in variants_gene.keys():
            if variants_gene['rs12769205'] == '1/1' and variants_gene['rs4244285'] == '1/1':
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 2 or (len(variants_gene) == 3 and 'rs3758581' in variants_gene.keys()):
                    diplotype = '*1/*2'
                else:
                    if 'rs4986893' in variants_gene.keys():
                        diplotype = '*2/*3'
                    elif 'rs12248560' in variants_gene.keys():
                        diplotype = '*2/*17'
                    else:
                        diplotype = 'NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        # Check *3 allele
        elif 'rs4986893' in variants_gene.keys():
            if variants_gene['rs4986893'] == '1/1':
                diplotype = '*3/*3'
            else:
                if len(variants_gene) == 1 or (len(variants_gene) == 2 and 'rs3758581' in variants_gene.keys()):
                    diplotype = '*1/*3'
                elif 'rs12248560' in variants_gene.keys():
                    diplotype = '*3/*17'
                else:
                    diplotype = 'NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')

        # Check *17 allele
        elif 'rs12248560' in variants_gene.keys():
            if variants_gene['rs12248560'] == '1/1':
                diplotype = '*17/*17'
            elif len(variants_gene) == 1 or (len(variants_gene) == 2 and 'rs3758581' in variants_gene.keys()):
                diplotype = '*1/*17'
            else:
                diplotype = 'NA'
                print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        else:
            diplotype = 'NA'
            print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')

    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results
