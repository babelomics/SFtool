
from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_tpmt_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Assing a diplotype to TPMT according to variants found in this gene and recommendations based on SEFF
    Essential alleles for SEFF: *2A, *3A, *3B, *3C, *4

    Args:
        variants (list): A list with a dictionary of variants present in the set of pharmacogenetic genes
        diplo_pheno_dct (dict): Dictionary with information about diplotype-phenotype relationships
        aggregated_results (list): A list with results of diplotype information accumulated for each gene

    Returns:
        list: List with results of diplotype information accumulateed for each gene, including CUP2C9

    """
    # Gene
    gene = 'TPMT'

    # Check and get variants for 'gene'
    found, variants_gene = check_gene_variants(variants, gene)

    # Diplotype *1/*1 (reference) if there are no variants on current gene
    if not found:
        diplotype = '*1/*1'
    else:
        # Check *2:
        if 'rs1800462' in variants_gene.keys():
            if variants_gene['rs1800462'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2'
                elif len(variants_gene) == 2:
                    if 'rs1800460' in variants_gene.keys():
                        diplotype = '*2/*3B'
                    elif 'rs1142345' in variants_gene.keys():
                        diplotype = '*2/*3C'
                    elif 'rs1800584' in variants_gene.keys():
                        diplotype = '*2/*4'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*2/*3A'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
                else:
                    diplotype='NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')

        elif 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys(): # Check *3A
            if variants_gene['rs1800460'] == '1/1' and variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3A'
            elif variants_gene['rs1800460'] == '0/1' and variants_gene['rs1142345'] == '0/1' and len(variants_gene) == 2:
                diplotype = '*1/*3A or *3B/*3C'
            elif variants_gene['rs1800460'] == '1/1' and variants_gene['rs1142345'] == '0/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3B'
            elif variants_gene['rs1800460'] == '0/1' and variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3C'
            elif 'rs1800584' in variants_gene.keys() and len(variants_gene) == 3:
                diplotype = '*3A/*4'
            else:
                diplotype='NA'
                print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        elif 'rs1800460' in variants_gene.keys(): # Check *3B:
            if variants_gene['rs1800460'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*3B/*3B'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3B'
                elif len(variants_gene) == 2:
                    if 'rs1142345' in variants_gene.keys() and (variants_gene['rs1800460'] == '0/1') and (variants_gene['rs1142345'] == '0/1'):
                        diplotype = '*1/*3A or *3B/*3C'
                    elif 'rs1142345' in variants_gene.keys() and (variants_gene['rs1800460'] == '1/1') and (variants_gene['rs1142345'] == '0/1'):
                        diplotype = '*3A/*3B'
                    elif 'rs1800584' in variants_gene.keys():
                        diplotype = '*3B/*4'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
                else:
                    diplotype='NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        elif 'rs1142345' in variants_gene.keys(): # Check *3C:
            if variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*3C/*3C'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3C'
                elif len(variants_gene) == 2:
                    if 'rs1800584' in variants_gene.keys():
                        diplotype = '*3C/*4'
                    else:
                        diplotype='NA'
                        print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
                else:
                    diplotype='NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
        elif 'rs1800584' in variants_gene.keys(): # Check *4:
            if variants_gene['rs1800584'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*4/*4'
            else:
                if len(variants_gene) == 1 and variants_gene['rs1800584'] == '0/1':
                    diplotype = '*1/*4'
                elif len(variants_gene) == 2:
                    diplotype='NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')
                else:
                    diplotype='NA'
                    print('There are variants in ' + gene + ' that are not considered in the diplotype assignment. Please, review variants manually')

    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results