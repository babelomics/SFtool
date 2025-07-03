
import json
import pandas as pd

def collapse_val(val):
    """
    If val is a dict with 'Min' and 'Max', collapse to "[Min,Max]".
    Otherwise, convert val to string.
    """
    if isinstance(val, dict) and "Min" in val and "Max" in val:
        return f"[{val['Min']},{val['Max']}]"
    return str(val)



def run_STRipy_module(reproductive_risk_geneset_STR_file, STRipy_output_file):
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
        STR_genes = list(rr_STRs_info.keys())
    else:
        rr_STRs_info = rr_STRs_info.to_dict(orient='index')


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
            "Inheritance": inheritance
        }


    return STRipy_results

