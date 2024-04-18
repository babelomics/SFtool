

from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_cyp2c9_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Asigna un diplotipo a CYP2C9 basado en las variantes genéticas y almacena los resultados.

    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        aggregated_results (list): Una lista de resultados donde se agregarán los resultados de la asignación.

    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de CYP2C9.

    Raises:
        TypeError: Si 'variants' no es una lista, si 'diplo_pheno_dct' no es un diccionario o si 'results' no es una lista.
    """
    # Control de errores para los argumentos de entrada
    if not isinstance(variants, list):
        raise TypeError("El argumento 'variants' debe ser una lista de variantes genéticas.")
    if not isinstance(diplo_pheno_dct, dict):
        raise TypeError("El argumento 'diplo_pheno_dct' debe ser un diccionario que contiene información sobre diplotipos y fenotipos genéticos.")
    if not isinstance(aggregated_results, list):
        raise TypeError("El argumento 'results' debe ser una lista para almacenar los resultados de la asignación de diplotipos.")

    # Gen
    gene = 'CYP2C9'

    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)

    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'

    else:
        if 'rs1799853' in variants_gene.keys():
            if variants_gene['rs1799853'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2'
                else:
                    if 'rs1057910' in variants_gene.keys():
                        diplotype = '*2/*3'
                    else:
                        diplotype = 'NA'
                        print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        elif 'rs1057910' in variants_gene.keys():
            if variants_gene['rs1057910'] == '1/1':
                diplotype = '*3/*3'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3'
                else:
                    diplotype = 'NA'
                    print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        else:
            diplotype = 'NA'
            print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results