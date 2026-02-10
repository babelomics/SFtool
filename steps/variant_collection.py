from modules.context import ExecutionContext
from modules.misc.geneBe_utils import parse_genebe_output
from modules.misc.utils import combine_genebe_clinvar_results
from modules.STRipy.STR_collection import STR_collection
from modules.SMAca.SMN1_collection import SMN1_collection
from modules.misc.pharmcat_utils import pharmCAT_collection
import json

def run(ctx: ExecutionContext) -> None:
    """
    Variant collection step to gather variants of interest from different evidences (GeneBe, Clinvar, STRipy and SMAca)
    """

    profile = ctx.profile
    CATEGORY_GENESETS = {
        "PR": ctx.config.catalogs.personal_risk_geneset,
        "RR": ctx.config.catalogs.reproductive_risk_geneset,
    }

    for sample in ctx.samples:
        categories = sample.categories
        for category in categories:
            if category == "PR" or category == "RR":
                category_geneset_file = CATEGORY_GENESETS[category]
                # Collect GENEBE variants (SNV/Indels)
                genebe_file = sample.vcf_outputs["genebe_annotated"][category]
                snv_indels_genebe = parse_genebe_output(genebe_file, profile, category, category_geneset_file)
                # Collect CLINVAR variants (SNV/Indels) and merge with GENEBE variants (only of profile is advanced)
                if profile == "advanced":
                    json_category = category + '_json'
                    clinvar_file = ctx.outputs["clinvar"][json_category]
                    with clinvar_file.open("r", encoding="utf-8") as fh:
                        snv_indels_clinvar = json.load(fh)
                    snv_indels_genebe_clinvar = combine_genebe_clinvar_results(snv_indels_genebe, snv_indels_clinvar)
                    sample.variant_collections[category]["snv_indels_genebe_clinvar"] = snv_indels_genebe_clinvar
                else:
                    sample.variant_collections[category]["snv_indels_genebe_clinvar"] = snv_indels_genebe
                # Collect STR results and SMN1 copy results (only for RR)
                if category == 'RR':
                    STRipy_file = sample.stripy_path
                    if STRipy_file != "None":
                        reproductive_risk_geneset_STR_file = ctx.config.catalogs.reproductive_risk_geneset_STR
                        sample.variant_collections[category]["STRs"] = STR_collection(reproductive_risk_geneset_STR_file, STRipy_file)
                    SMAca_file = sample.smaca_path
                    if SMAca_file != "None":
                        sample.variant_collections[category]["SMN1_copy"] = SMN1_collection(SMAca_file, ctx.config.smaca_thresholds)
            elif category == "PGx":
                sample.variant_collections[category]["pharmCAT_variants"] = pharmCAT_collection(sample.results["PGx"])
