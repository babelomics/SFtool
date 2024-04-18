
from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_dpyd_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Asigna un diplotipo a DPYD basado en las variantes genéticas y almacena los resultados.

    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        aggregated_results (list): Una lista de resultados donde se agregarán los resultados de la asignación.

    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de DPYD.
    """
    # Gen
    gene = 'DPYD'

    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)

    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'

    else:
        # Chequear *2A:
        if 'rs3918290' in variants_gene.keys():
            if variants_gene['rs3918290'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
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
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = '*2a/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *13:
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
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = '*13/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        # Chequear c.2846A>T:
        elif 'rs67376798' in variants_gene.keys():
            if variants_gene['rs67376798'] == '1/1':
                diplotype = 'c.2846A>T/c.2846A>T'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/c.2846A>T'
                elif len(variants_gene) == 2:
                    print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = 'c.2846A>T/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        # Chequear c.2846A>T:
        elif 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
            if variants_gene['rs75017182'] == '1/1' and variants_gene['rs56038477'] == '1/1':
                diplotype = 'HapB3/HapB3'
            else:
                if len(variants_gene) == 2:
                    diplotype = '*1/HapB3'
                elif len(variants_gene) == 3:
                    print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results
