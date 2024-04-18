
from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_nudt15_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Asigna un diplotipo a NUDT15 basado en las variantes genéticas y almacena los resultados.

    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        aggregated_results (list): Una lista de resultados donde se agregarán los resultados de la asignación.

    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de NUDT15.
    """
    # Gen
    gene = 'NUDT15'

    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)

    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'

    else:
        # Chequear *2A:
        if 'rs116855232' in variants_gene.keys():
            if len(variants_gene) == 1:
                if variants_gene['rs116855232'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                    diplotype = '*3/*3'
                else:
                    diplotype = '*1/*3'
            elif len(variants_gene) == 2:
                if 'rs746071566' in variants_gene.keys():
                    if variants_gene['rs746071566'] == '1/1' and variants_gene['rs116855232'] == '1/1':
                        diplotype = '*2/*2'
                    elif variants_gene['rs746071566'] == '0/1' and variants_gene['rs116855232'] == '1/1': #hbaría que chequear variante si es inserción o del
                        diplotype = '*2/*3'
                    elif variants_gene['rs746071566'] == '0/1' and variants_gene['rs116855232'] == '0/1':
                        diplotype = '*1/*2'
                    else:
                        print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
            else:
                print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        else:
            print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results
