
import json
import pandas as pd
import re

def collapse_val(val):
    """
    If val is a dict with 'Min' and 'Max', collapse to "[Min,Max]".
    Otherwise, convert val to string.
    """
    if isinstance(val, dict) and "Min" in val and "Max" in val:
        return f"[{val['Min']},{val['Max']}]"
    return str(val)


def parse_genotype(genotype: str) -> list[int]:
    """
    Convert STRipy genotype string into a list of allele repeat sizes.

    Examples:
        "34"      -> [34]
        "22/45"   -> [22, 45]
    """
    if not genotype:
        return []

    return [int(allele) for allele in genotype.split("/") if allele.isdigit()]

def parse_threshold(threshold: str):
    """
    Parse pathogenic threshold string.

    Supported formats:
        "200"     -> allele > 200
        "[6-25]"  -> 6 <= allele <= 25
    """
    threshold = threshold.strip()

    # Interval: [a-b]
    interval_match = re.fullmatch(r"\[(\d+),(\d+)\]", threshold)
    if interval_match:
        low, high = map(int, interval_match.groups())
        return ("range", low, high)

    # Single integer
    if threshold.isdigit():
        return ("gte", int(threshold))

    raise ValueError(f"Unsupported PathogenicThreshold format: {threshold}")

def is_pathogenic_allele(allele: int, threshold: str) -> bool:
    mode, *values = parse_threshold(threshold)

    if mode == "gte":
        return allele > values[0]

    if mode == "range":
        low, high = values
        return low <= allele

    return False


def is_pathogenic_threshold(genotype: str, threshold: str) -> bool:

    alleles = parse_genotype(genotype)
    return any(is_pathogenic_allele(a, threshold) for a in alleles)




def STR_collection(reproductive_risk_geneset_STR_file, STRipy_output_file):
    '''
    Parse STRipy's JSON file

    :param reproductive_risk_geneset_STR_file:
    :param STRipy_output:
    :return:
    '''


    # Parse JSON file from STRipy
    with open(STRipy_output_file, "r") as fd:
        STRipy_results = json.load(fd)

    # Read reproductive risk geneset related to pathogenic STRs
    rr_STRs_info = pd.read_csv(reproductive_risk_geneset_STR_file)
    STR_genes = list(rr_STRs_info['Gene'])

    rr_STRs_info_indexed = rr_STRs_info.set_index('Gene')

    if 'ARX' in STR_genes: # ARX is split into ARX_1 and ARX_2 in STRipy results
        arx_entry = rr_STRs_info_indexed.loc['ARX']
        rr_STRs_info_indexed = rr_STRs_info_indexed.drop('ARX')

        rr_STRs_info_indexed.loc['ARX_1'] = arx_entry
        rr_STRs_info_indexed.loc['ARX_2'] = arx_entry

        rr_STRs_info = rr_STRs_info_indexed.to_dict(orient='index')
    else:
        rr_STRs_info = rr_STRs_info.to_dict(orient='index')
    STR_genes = list(rr_STRs_info.keys())

    # Get results from STRipy related to the set of genes contained in STR_genes
    matched = []
    for entry in STRipy_results["GenotypingResults"]:
        STRipy_gene = next(iter(entry))
        if STRipy_gene in STR_genes:
            matched.append((STRipy_gene, entry[STRipy_gene]))


    # Access information for the genes of interest
    STRipy_results = {}
    for current_gene, STRipy_info in matched:
        # Get Repeat values of each allele
        repeats = "/".join(str(allele["Repeats"]) for allele in STRipy_info.get("Alleles", "" ))
        # Get result for each allele
        ranges = "/".join(str(allele["Range"]) for allele in STRipy_info.get("Alleles", ""))
        # Get outliers information
        outliers = "/".join(str(allele["IsPopulationOutlier"]) for allele in STRipy_info.get("Alleles", ""))
        # Get coordinates
        coords = STRipy_info["TargetedLocus"]["Coordinates"]
        # Get motif
        motif = STRipy_info["TargetedLocus"]["Motif"]
        # Find matching disease according to rr_STRs_info
        for did, dinfo in STRipy_info["TargetedLocus"]["CorrespondingDisease"].items():
            if dinfo["DiseaseOMIM"] == str(rr_STRs_info.get(current_gene)['OMIM Disorder']):
                intermediate_threshold = collapse_val(dinfo.get("IntermediateRange"))
                normal_threshold = collapse_val(dinfo.get("NormalRange"))
                pathogenic_threshold  = collapse_val(dinfo.get("PathogenicCutoff"))
                omim_disorder = rr_STRs_info.get(current_gene)['OMIM Disorder']
                phenotype = rr_STRs_info.get(current_gene)['Phenotype']
                inheritance = rr_STRs_info.get(current_gene)['Inheritance']


        current_threshold = pathogenic_threshold
        if current_gene == 'FMR1': # For FMR1, intermediate threshold is used
            current_threshold = intermediate_threshold


        if is_pathogenic_threshold(repeats,current_threshold): # Only store those entries above the threshold (pathogenic or intermediate)
            STRipy_results[coords] = {
                "Gene": current_gene,
                "Genotype": repeats,
                "Motif": motif,
                "Outliers": outliers,
                "NormalThreshold": normal_threshold,
                "IntermediateThreshold": intermediate_threshold,
                "PathogenicThreshold": pathogenic_threshold,
                "Result": ranges,
                "Filter": STRipy_info["Filter"],
                "Phenotype": phenotype,
                "OMIMdisorder": omim_disorder,
                "Inheritance": inheritance,
                "related_HPOs_for_sample": 'NA'
            }


    return STRipy_results

