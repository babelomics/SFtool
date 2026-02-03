# Set of utils functions

from modules.STRipy.STR_collection import parse_threshold

def check_specific_criteria(gene, variant, variant_key, assembly):
    """
    Check specific features for given gene: either specific consequence or specific variant

    :param gene: gene information from current catalogue
    :param variant: variant information under study
    :param variant_key: specific genomic variant under study
    :param assembly: genome assembly version, either 37 or 38
    :return: TRUE if criteria is met, else, FALSE
    """

    meet_criteria = True

    if gene["specific_consequence"] != "":
        variant_consequences = variant["Consequence"].split(",")
        gene_consequences = gene["specific_consequence"].split(",")

        # If required consequence from gene information is not present in the variant consequence, criteria is not meet (variant is discarded)
        if len(set(variant_consequences) & set(gene_consequences)) == 0:
            meet_criteria = False
    elif gene["specific_variant_" + str(assembly)] != "":
        specific_variant = gene["specific_variant_" + str(assembly)].split(',')[0]
        specific_genotype = gene["specific_variant_" + str(assembly)].split(',')[1]

        if specific_variant != variant_key or specific_genotype.lower() != variant["Genotype"]:
            meet_criteria = False

    return meet_criteria

def combine_variant_and_gene_info(variant_info, gene_info):
    """
    Combine information of variant and gene

    Args:
        variant_info (dict): Variant info
        gene_info (dict): Gene info

    Returns:
        dict: Merged information
    """
    combined_info = {
        "Gene": variant_info["Gene"],
        "Genotype": variant_info["Genotype"],
        "rs": variant_info.get("rs", ""),
        "Transcript": variant_info.get("Transcript",""),
        "HGVSC": variant_info.get("HGVSC",""),
        "HGVSP": variant_info.get("HGVSP",""),
        "Consequence": variant_info.get("Consequence", ""),
        "GeneBe_ACMG_Classification": variant_info["GeneBeClassification"],
        "GeneBe_ACMG_criteria": variant_info.get("ACMG_criteria", "-"),
        "ClinvarClinicalSignificance": variant_info.get("ClinvarClinicalSignificance", "-"),
        "ReviewStatus": variant_info.get("ReviewStatus", "-"),
        "ClinvarSummary": variant_info.get("ClinvarSummary", "-"),
        "ClinvarID": variant_info.get("ClinvarID", "-"),
        "Orpha": variant_info.get("Orpha", ""),
        "Phenotype": gene_info["phenotype"],
        "ACMG_version": gene_info.get("ACMG_version", ""),  # Usar get para manejar la falta de 'ACMG_version'
        "OMIM_disorder": gene_info["OMIM_disorder"],
        "inheritance": gene_info["inheritance"],
        "variants_to_report": gene_info.get("variants_to_report", ""),  # Usar get para manejar la falta de 'variants_to_report'
        "related_HPOs_for_sample": 'NA',
        "VCF_Sample_FORMAT": variant_info["VCFSampleFormat"]
    }
    return combined_info

def is_autosomal_chromosome(chrom: str) -> bool:
    """
    Return True if the chromosome corresponds to a human autosome (1–22).

    Accepted inputs:
      - "1" .. "22"
      - "chr1" .. "chr22"
      - case-insensitive ("CHR1", "Chr22", etc.)

    Returns False for sex chromosomes (X, Y), MT/M, or invalid values.
    """
    if chrom is None:
        return False

    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    try:
        chrom_num = int(chrom)
    except ValueError:
        return False

    return 1 <= chrom_num <= 22

def index_by_gene_symbol(entries: list[dict]) -> dict[str, dict]:
    """
    Build a dictionary indexed by gene_symbol.
    """
    return {
        entry["gene_symbol"]: entry
        for entry in entries
        if "gene_symbol" in entry
    }

def classify_genotype(gt: str) -> str | None:
    """
    Classify a genotype string as heterozygous or homozygous.

    Parameters
    ----------
    gt : str
        Genotype string (e.g. '0/1', '1|1', '1/0', '0|0').

    Returns
    -------
    str | None
        'het' for heterozygous genotypes,
        'hom' for homozygous genotypes,
        None if genotype is missing or unsupported.
    """
    if not gt or gt in {".", "./.", ".|."}:
        return None

    # Normalize separator
    sep = "/" if "/" in gt else "|" if "|" in gt else None
    if sep is None:
        return None

    alleles = gt.split(sep)
    if len(alleles) != 2 or "." in alleles:
        return None

    return "HOM" if alleles[0] == alleles[1] else "HET"


def classify_pathogenic_STRs(repeats: str, threshold: str) -> bool:
    """
        Classify repeat values relative to a threshold.

        Parameters
        ----------
        repeats : str
            Repeat values in the form 'value1/value2' (e.g. '28/35').
        threshold : int
            Threshold to classify repeats.

        Returns
        -------
        str | None
            'HOM' if both values are >= threshold,
            'HET' if only one value is >= threshold,
            None if neither value meets the threshold or input is invalid.
        """
    if not repeats or "/" not in repeats:
        return None

    try:
        mode, *values = parse_threshold(threshold)
        v1, v2 = (int(v) for v in repeats.split("/", 1))
    except ValueError:
        return None


    if mode == "gte":
        current_threshold = values[0]
    else: #mode = range
        low, high = values
        current_threshold = low

    if v1 >= current_threshold and v2 >= current_threshold:
        return "HOM"
    if v1 >= current_threshold or v2 >= current_threshold:
        return "HET"

    return None

