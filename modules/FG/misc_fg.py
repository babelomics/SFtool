
import json
import csv
import vcfpy
import os
import re
import pandas as pd


def check_gene_variants(variants, gene):
    """
    Check the presence of variants in a given pharma gene. Genotype info is retrieved if variant exists

    Args:
        variants (list): A list of dictionaries that contain information of genetic variants
        gene (str): Gene name to be found

    Returns:
        tuple: tuple with two elements:
            - found (bool): True if at least a variant is found in the given gene. False otherwise
            - variants_gene (dict): Dictionary that stores genotype for the variants found in the gene. Key is rs and genotype is the corresponding value for the key

    Raises:
        TypeError: Error if 'variants' is not a list or if 'gene' is not a string
    """
    # Error control
    if not isinstance(variants, list):
        raise TypeError("Argument 'variants' must be a list of genetic variants.")
    if not isinstance(gene, str):
        raise TypeError("Argument 'gene' must be a string.")

    found = False
    variants_gene = {}

    for variant in variants:
        if variant['Gene Symbol'] == gene:
            found = True
            variants_gene[variant['rs']] = variant['Genotype']

    return found, variants_gene

def annotate_fg_variants(categories_path, norm_vcf_fg, assembly, temp_path):
    """
    Annote a list of variants from a VCF file given a JSON file with information about pharma variants

    Args:
        categories_path (str): Path to directory with files related to pharma category
        norm_vcf_fg (str): Path to VCF file that contains normalized variants in genes related to pharma
        assembly (str): Assembly version
        temp_path (str): Path to temporal directory

    Returns:
        list: The list of annotated variants

    Raises:
        FileNotFoundError: I/O error
    """

    try:
        # Load JSON file with the list of pharma variants
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
                    "Genotype": variant_record.call_for_sample[sample_name].data.get('GT'),
                    "Gene Symbol": variant_info["gene_symbol"],
                    "rs": variant_info["rs"]
                }
            else:
                annotated_variant = {
                    "Variant": variant_key,
                    "Genotype": variant_record.call_for_sample[sample_name].data.get('GT'),
                    "Gene Symbol": "Not found",
                    "rs": "Not found",
                }

            annotated_variants.append(annotated_variant)

        return annotated_variants

    except FileNotFoundError as e:
        print(f"Error when annotating pharma variants: {e}")
        return []


def get_diplotype_phenotype_dictionary(diplotype_phenotype_info_file):
    """
    Build a dictionary with diplotype-phenotype relationships and their corresponding Activity Score


    Args:
        diplotype_phenotype_info_file (str): Path to diplotype-phenotype info file.

    Returns:
        dict: A dictionary in which a give diplotype has phenotype and activity score info
    """

    diplotype_data = {}

    try:
        with open(diplotype_phenotype_info_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gene = row["GENE"]
                diplotype = row["DIPLOTYPE"]
                phenotype = row["Phenotype"]
                activity_score = row["Activity_Score"]
                # Check whether the gene is already in the dictionary
                if gene in diplotype_data:
                    # If gene exists in dictionary, add another diplotype
                    diplotype_data[gene][diplotype] = {"Phenotype": phenotype, "Activity_Score": activity_score}
                else:
                    # If gene does not exist in dictionary, create a new entry and add the corresponding diplotype
                    diplotype_data[gene] = {diplotype: {"Phenotype": phenotype, "Activity_Score": activity_score}}
    except FileNotFoundError:
        raise FileNotFoundError(f"El archivo CSV no se encuentra en la ruta especificada: {diplotype_phenotype_info_file}")

    return diplotype_data


def assign_phenotype_AC(diplotype, gene, diplo_pheno_dct, aggregated_results):
    """
    According to a set of diplotypes in a gene, assign phenotype information and activity scores

    :param diplotype: dictionary with diplotype information for a given gene
    :param gene: Gene
    :param diplo_pheno_dct: dictionary with diplotype-phenotype information
    :param aggregated_results: aggregated results with diplotype-phenotype information from other genes
    :return: arregated (from a set of genes) results with diplotype-phenotype information
    """

    if gene in diplo_pheno_dct and 'or' in diplotype: # There are questionable diplotypes
        all_diplotypes = re.findall(r'\*\d+/\*\d+', diplotype)
        tmp_phenotype = []
        tmp_activity_score = []
        for current_diplotype in all_diplotypes:
            if current_diplotype in diplo_pheno_dct[gene]:
                data = diplo_pheno_dct[gene][current_diplotype]
                tmp_phenotype.append(data['Phenotype'])
                tmp_activity_score.append(data['Activity_Score'])
            else:
                tmp_phenotype.append('NA')
                tmp_activity_score.append('NA')
        phenotype = ', '.join(tmp_phenotype[:-1]) + ' or ' + tmp_phenotype[-1]
        activity_score = ', '.join(tmp_activity_score[:-1]) + ' or ' + tmp_activity_score[-1]
    elif gene in diplo_pheno_dct and diplotype in diplo_pheno_dct[gene]:
        data = diplo_pheno_dct[gene][diplotype]
        phenotype = data['Phenotype']
        activity_score = data['Activity_Score']
    else:
        phenotype = 'NA'
        activity_score = 'NA'

    # Append results to the list dictionary
    aggregated_results.append({
        'Gene Symbol': gene,
        'Diplotype': diplotype,
        'Phenotype': phenotype,
        'Activity Score': activity_score
    })

    return aggregated_results

def write_fg_results_to_tsv(fg_annotated_variants, aggregated_fg_results, output_file_variants, output_file_diplopheno):
    """
    Write pharmacogenietic results to two TSV files

    :param fg_annotated_variants: Dictionary with variants found in VCF file that are in the Pharma list
    :param aggregated_fg_results: Dictionary with variants in CYP2C9, CYP2C19, DPYD, NUDT15 and TPMT genes with Diplotype and Phenotype info
    :param output_file_variants: Output TSV file for fg_annotated_variants
    :param output_file_diplopheno: Output TSV file for aggregated_fg_results
    :return:
    """


    fg_df = pd.DataFrame(fg_annotated_variants)
    fg_df.to_csv(output_file_variants, index=False, sep='\t')

    haplot_df = pd.DataFrame(aggregated_fg_results)
    haplot_df.to_csv(output_file_diplopheno, index=False, sep='\t')
