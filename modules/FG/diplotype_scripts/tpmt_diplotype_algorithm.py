
from modules.FG.misc_fg import check_gene_variants
from modules.FG.misc_fg import assign_phenotype_AC

def assign_tpmt_diplotype(variants, diplo_pheno_dct, aggregated_results):
    """
    Asigna un diplotipo a TPMT basado en las variantes genéticas y almacena los resultados.

    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        diplo_pheno_dct (dict): Un diccionario que contiene información sobre diplotipos y fenotipos genéticos.
        aggregated_results (list): Una lista de resultados donde se agregarán los resultados de la asignación.

    Returns:
        list: La lista de resultados actualizada después de agregar el resultado de la asignación del diplotipo de TPMT.
    """
    # Gen
    gene = 'TPMT'

    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)

    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'

    else:
        # Chequear *2:
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
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*2/*3A'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *3B:
        elif 'rs1800460' in variants_gene.keys():
            if variants_gene['rs1800460'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*3B/*3B'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3B'
                elif len(variants_gene) == 2:
                    if 'rs1142345' in variants_gene.keys() and (variants_gene['rs1800460'] == '0/1') and (variants_gene['rs1142345'] == '0/1'):
                        diplotype = '*1/*3A o *3B/*3C'
                    elif 'rs1800584' in variants_gene.keys():
                        diplotype = '*3B/*4'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*3B'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *3C:
        elif 'rs1142345' in variants_gene.keys():
            if variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*3C/*3C'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3C'
                elif len(variants_gene) == 2:
                    if 'rs1800584' in variants_gene.keys():
                        diplotype = '*3C/*4'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*3B'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *4:
        elif 'rs1800584' in variants_gene.keys():
            if variants_gene['rs1800584'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*4/*4'
            else:
                if len(variants_gene) == 1 and variants_gene['rs1800584'] == '0/1':
                    diplotype = '*1/*4'
                elif len(variants_gene) == 2:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*4'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *3A:  IGUAL REORGANIZO TODO TPMT Y ESTO AL INICIO
        elif 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
            if variants_gene['rs1800460'] == '1/1' and variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3A'
            elif variants_gene['rs1800460'] == '1/1' and variants_gene['rs1142345'] == '0/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3B'
            elif variants_gene['rs1800460'] == '0/1' and variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3C'
            elif variants_gene['rs1800460'] == '0/1' and variants_gene['rs1142345'] == '0/1' and len(variants_gene) == 2:
                diplotype = '*1/*3A o *3A/*3C'
            else:
                print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

    # Phenotype and Activity Score assignment for the obtained diplotype
    aggregated_results = assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results)

    return aggregated_results