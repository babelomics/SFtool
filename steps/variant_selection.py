from modules.context import ExecutionContext
from modules.variant_selection.pr_variant_selection import pr_variant_selection


def run(ctx: ExecutionContext) -> None:
    for sample in ctx.samples:
        categories = sample.categories
        for category in categories:
            if category == "PR":
                snv_indels_pr_collection = sample.variant_collections[category]["snv_indels_genebe_clinvar"]
                json_file = ctx.outputs["catalogs"]["json_files"][category]
                assembly = ctx.assembly
                snv_indels_pr_selection = pr_variant_selection(snv_indels_pr_collection, json_file, assembly)
                sample.variant_selection[category]["snv_indels_genebe_clinvar"] = snv_indels_pr_selection
