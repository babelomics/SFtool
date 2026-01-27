from modules.context import ExecutionContext
from modules.variant_selection.pr_variant_selection import pr_variant_selection
from modules.variant_selection.rr_variant_selection import rr_variant_selection


def run(ctx: ExecutionContext) -> None:
    for sample in ctx.samples:
        categories = sample.categories
        for category in categories:
            if category == "PR" or category == 'RR':
                snv_indels_collection = sample.variant_collections[category]["snv_indels_genebe_clinvar"]
                json_file = ctx.outputs["catalogs"]["json_files"][category]
                if category == "PR":
                    snv_indels_selection = pr_variant_selection(snv_indels_collection, json_file, ctx.assembly)
                elif category == "RR":
                    snv_indels_selection = rr_variant_selection(snv_indels_collection, json_file, ctx.RR_mode, sample.sex)
                sample.variant_selection[category]["snv_indels_genebe_clinvar"] = snv_indels_selection


