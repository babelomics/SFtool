from modules.context import ExecutionContext
from modules.variant_selection.pr_variant_selection import pr_variant_selection
from modules.variant_selection.rr_variant_selection import rr_variant_selection, rr_str_selection, rr_smn1_copy_selection


def run(ctx: ExecutionContext) -> None:
    for sample in ctx.samples:
        categories = sample.categories
        for category in categories:
            if category == "PR" or category == 'RR':
                snv_indels_collection = sample.variant_collections[category]["snv_indels_genebe_clinvar"]
                json_file = ctx.outputs["catalogs"]["json_files"][category]
                if category == "PR":
                    # SNV / Indels
                    snv_indels_selection = pr_variant_selection(snv_indels_collection, json_file, ctx.assembly)
                elif category == "RR":
                    # SNV / Indels
                    snv_indels_selection = rr_variant_selection(snv_indels_collection, json_file, ctx.RR_mode, sample.sex)
                    # STRs
                    if sample.variant_collections[category]["STRs"]:
                        str_collection = sample.variant_collections[category]["STRs"]
                        str_selection = rr_str_selection(str_collection, ctx.RR_mode, sample.sex)
                        sample.variant_selection[category]["STRs"] = str_selection
                    # SMN1-copy
                    if sample.variant_collections[category]["SMN1_copy"]:
                        SMN1_copy_collection = sample.variant_collections[category]["SMN1_copy"]
                        SMN1_copy_selection = rr_smn1_copy_selection(SMN1_copy_collection, ctx.RR_mode)
                        sample.variant_selection[category]["SMN1_copy"] = SMN1_copy_selection
                else:
                    continue

                sample.variant_selection[category]["snv_indels_genebe_clinvar"] = snv_indels_selection



