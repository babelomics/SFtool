
import json
import csv
import vcfpy
import os


def check_gene_variants(variants, gene):
    """
    Verifica la presencia de variantes en un gen específico y recopila información sobre su genotipo.

    Args:
        variants (list): Una lista de diccionarios que contienen información de variantes genéticas.
        gene (str): El nombre del gen objetivo a buscar en las variantes.

    Returns:
        tuple: Una tupla que contiene dos elementos:
            - found (bool): True si se encontró al menos una variante en el gen objetivo, False en caso contrario.
            - variants_gene (dict): Un diccionario que almacena el genotipo de las variantes encontradas en el gen, donde
              la clave es el rsID de la variante y el valor es el genotipo.

    Raises:
        TypeError: Si 'variants' no es una lista o si 'gene' no es una cadena de caracteres.
    """
    # Control de errores para los argumentos de entrada
    if not isinstance(variants, list):
        raise TypeError("El argumento 'variants' debe ser una lista de variantes genéticas.")
    if not isinstance(gene, str):
        raise TypeError("El argumento 'gene' debe ser una cadena de caracteres que representa el gen objetivo.")

    # Flag para verificar si se encontró una variante en el gen objetivo
    found = False
    # Crear diccionario para almacenar genotipo de variantes del gene       #igual mejor guardar la variante también, para el caso de indels con mismo rs y diferente nº de repeticiones
    variants_gene = {}

    # Iterar a través de las variantes
    for variant in variants:
        if variant['Gene'] == gene:
            found = True
            variants_gene[variant['rs']] = variant['GT']

    return found, variants_gene

def annotate_fg_variants(categories_path, norm_vcf_fg, assembly, temp_path):
    """
    Anota variantes genéticas utilizando un archivo JSON de variantes farmacogenéticas.

    Args:
        categories_path (str): Ruta al directorio que contiene archivos relacionados con las categorías.
        norm_vcf_fg (str): Ruta al archivo VCF normalizado e interseccionado.
        assembly (str): Ensamblaje genómico a utilizar.
        temp_path (str): Ruta al directorio que contiene archivos intermedios.

    Returns:
        list: Una lista de diccionarios que contienen los resultados de la anotación.

    Raises:
        FileNotFoundError: Si no se puede encontrar el archivo JSON con variantes farmacogenéticas.
    """

    try:
        # Cargar archivo fg_json con variantes farmacogenéticas
        fg_json_path = f'{categories_path}FG/fg_risk_genes_GRCh{assembly}.json'
        with open(fg_json_path, 'r') as file:
            fg_json = json.load(file)

        annotated_variants = []

        # Read VCF file
        vcf_reader = vcfpy.Reader.from_path(norm_vcf_fg)
        for variant_record in vcf_reader:
            chrom = str(variant_record.CHROM)
            pos = str(variant_record.POS)
            ref = str(variant_record.REF)
            alt = str(variant_record.ALT[0].value)
            sample_name = list(variant_record.call_for_sample.keys())[0]
            variant_key = chrom + ':' + pos + ':' + ref + ':' + alt
            variant_info = fg_json['variants'].get(variant_key, None)

            if variant_info:
                annotated_variant = {
                    "Variant": variant_key,
                    "GT": variant_record.call_for_sample[sample_name].data.get('GT'),
                    "Gene": variant_info["gene_symbol"],
                    "rs": variant_info["rs"]
                }
            else:
                annotated_variant = {
                    "Variant": variant_key,
                    "GT": variant_record.call_for_sample[sample_name].data.get('GT'),
                    "Gene": "Not found",
                    "rs": "Not found",
                }

            annotated_variants.append(annotated_variant)

        # Escribir los resultados en un archivo
        just_filename = os.path.basename(norm_vcf_fg)
        result_file = os.path.join(temp_path, just_filename.split(".norm.FG.vcf.gz")[0] + ".FG.annotated.json")

        with open(result_file, 'w') as fout:
            json.dump(annotated_variants, fout)

        return annotated_variants

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return []


def get_diplotype_phenotype_dictionary(categories_path):
    """
    Obtener un diccionario de diplotipos con fenotipos y puntuaciones de actividad.

    Args:
        categories_path (str): Ruta al directorio categories.

    Returns:
        dict: Un diccionario que asigna diplotipos con su información de fenotipo y puntuación de actividad.
    """
    # Nombre del archivo CSV
    csv_file = f'{categories_path}FG/diplotipo_fenotipo.csv'

    # Definir un diccionario para mapear los diplotipos con su fenotipo y Activity Score
    diplotype_data = {}

    try:
        # Abrir el archivo CSV y cargar los datos en el diccionario
        with open(csv_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gene = row["GENE"]
                diplotype = row["DIPLOTYPE"]
                phenotype = row["Phenotype"]
                activity_score = row["Activity_Score"]
                # Verificar si el gen ya está en el diccionario
                if gene in diplotype_data:
                    # Si el gen ya está en el diccionario, agregar el nuevo diplotipo
                    diplotype_data[gene][diplotype] = {"Phenotype": phenotype, "Activity_Score": activity_score}
                else:
                    # Si el gen no está en el diccionario, crear una entrada para el gen y agregar el diplotipo
                    diplotype_data[gene] = {diplotype: {"Phenotype": phenotype, "Activity_Score": activity_score}}
    except FileNotFoundError:
        raise FileNotFoundError(f"El archivo CSV no se encuentra en la ruta especificada: {csv_file}")

    return(diplotype_data)


def assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results):
    """

    :param diplotype:
    :param gene:
    :param diplo_pheno_dct:
    :param aggregated_results:
    :return:
    """

    if gene in diplo_pheno_dct and diplotype in diplo_pheno_dct[gene]:
        data = diplo_pheno_dct[gene][diplotype]
        phenotype = data['Phenotype']
        activity_score = data['Activity_Score']
    else:
        phenotype = 'NA'
        activity_score = 'NA'

    # Agregar los resultados a la lista de diccionarios
    aggregated_results.append({
        'Gene': gene,
        'Diplotipo': diplotype,
        'Phenotype': phenotype,
        'Activity Score': activity_score
    })

    return aggregated_results
