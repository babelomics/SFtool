
from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_cyp2c19_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Asigna un diplotipo a CYP2C19 basado en las variantes genéticas y almacena los resultados.

    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        aggregated_results (list): Una lista de resultados donde se agregarán los resultados de la asignación.

    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de CYP2C19.

    """
    # Gen
    gene = 'CYP2C19'

    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)

    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False or (len(variants_gene) == 1 and 'rs3758581' in variants_gene.keys()):
        diplotype = '*1/*1'

    else:
        if 'rs12769205' in variants_gene.keys() and 'rs4244285' in variants_gene.keys():
            if variants_gene['rs12769205'] == '1/1'  and variants_gene['rs4244285'] == '1/1': # en realidad habría que comprobar si es len >2, porque puede que esté en fase con otro alelo
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
                        print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        elif 'rs4986893' in variants_gene.keys():
            if variants_gene['rs4986893'] == '1/1':
                diplotype = '*3/*3'
            else:
                if len(variants_gene) == 1 or (len(variants_gene) == 2 and 'rs3758581' in variants_gene.keys()):
                    diplotype = '*1/*3'
                elif 'rs12248560' in variants_gene.keys():
                    diplotype = '*3/*17'
                else:
                    print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # me falta el *17
        elif 'rs12248560' in variants_gene.keys():
            if variants_gene['rs12248560'] == '1/1':
                diplotype = '*17/*17'
            elif len(variants_gene) == 1 or (len(variants_gene) == 2 and 'rs3758581' in variants_gene.keys()):
                diplotype = '*1/*17'
            else:
                print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')


    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results
