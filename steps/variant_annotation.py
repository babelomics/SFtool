
from modules.context import ExecutionContext
from modules.misc.geneBe_utils import run_genebe
from modules.misc.clinvar_utils import run_clinvar


def run(ctx: ExecutionContext) -> None:
    """
    Variant annotation with GeneBe and Clinvar
    """

    genebe_apikey = ctx.config.genebe_credentials.api_key
    genebe_username = ctx.config.genebe_credentials.username
    genebe_path = ctx.config.paths.genebe
    java_path = ctx.config.paths.java
    assembly = ctx.assembly


    for sample in ctx.samples:
        categories = sample.categories
        for category in categories:
            if category == 'PR' or category == 'RR':
                vcf_file = sample.vcf_outputs["intersected"][category]

                # 1. GeneBe annotation
                genebe_output_file = run_genebe(vcf_file, category, assembly, genebe_path, java_path, genebe_apikey, genebe_username)
                sample.vcf_outputs["genebe_annotated"][category] = genebe_output_file

    # # 2. Clinvar annotation (only in advanced mode)
    # if mode == 'advanced':
    #     clinvar_results = run_clinvar(evidence_level, clinvar_db, clinvar_submission, category, category_geneset_file)