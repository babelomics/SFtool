
from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_dpyd_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Assing a diplotype to DPYD according to variants found in this gene and recommendations based on SEFF
    Essential alleles for SEFF: *2A, *13, c.2846A>T, HapB3

    Args:
        variants (list): A list with a dictionary of variants present in the set of pharmacogenetic genes
        diplo_pheno_dct (dict): Dictionary with information about diplotype-phenotype relationships
        aggregated_results (list): A list with results of diplotype information accumulated for each gene

    Returns:
        list: List with results of diplotype information accumulateed for each gene, including CUP2C9

    """
    # Gene
    gene = 'DPYD'

    # Check and get variants for 'gene'
    found, variants_gene = check_gene_variants(variants, gene)

    # Diplotype *1/*1 (reference) if there are no variants on current gene
    if not found:
        diplotype = '*1/*1'
    else:
        # Check *2A:
        if 'rs3918290' in variants_gene.keys():
            if variants_gene['rs3918290'] == '1/1':
                diplotype = '*2A/*2A'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2A'
                elif len(variants_gene) == 2:
                    if 'rs55886062' in variants_gene.keys():
                        diplotype = '*2A/*13'
                    elif 'rs67376798' in variants_gene.keys():
                        diplotype = '*2A/c.2846A>T'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = '*2a/HapB3'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')

        # Check *13:
        elif 'rs55886062' in variants_gene.keys():
            if variants_gene['rs55886062'] == '1/1':
                diplotype = '*13/*13'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*13'
                elif len(variants_gene) == 2:
                    if 'rs67376798' in variants_gene.keys():
                        diplotype = '*13/c.2846A>T'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = '*13/HapB3'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        # Check c.2846A>T:
        elif 'rs67376798' in variants_gene.keys():
            if variants_gene['rs67376798'] == '1/1':
                diplotype = 'c.2846A>T/c.2846A>T'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/c.2846A>T'
                elif len(variants_gene) == 2:
                    diplotype='NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = 'c.2846A>T/HapB3'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        # Check c.2846A>T:
        elif 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
            if variants_gene['rs75017182'] == '1/1' and variants_gene['rs56038477'] == '1/1':
                diplotype = 'HapB3/HapB3'
            else:
                if len(variants_gene) == 2:
                    diplotype = '*1/HapB3'
                elif len(variants_gene) == 3:
                    diplotype='NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        else:
            diplotype='NA'
            print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')

    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results
