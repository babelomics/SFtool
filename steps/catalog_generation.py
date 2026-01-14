# SFtool/steps/catalog_generation.py

from modules.context import ExecutionContext
from modules.catalogs.build_json_bed_files import build_json_bed_files
from pathlib import Path


def run(ctx: ExecutionContext) -> None:
    """
    Build category-specific BED/JSON catalogs if missing
    and register outputs in ctx.outputs.catalogs
    """

    # Get variables needed from the ctx object
    catalogs_cfg = ctx.config.catalogs
    assembly = ctx.assembly

    # Get unique list of categories for all samples
    unique_categories = sorted({
        c for s in ctx.samples for c in s.categories
    })

    # To check whether chr prefix is present in VCF file, get the VCF file from the first sample
    sample = ctx.samples[0]
    vcf_file = str(sample.vcf)

    # ------------------------------------------------------------
    # Personal Risk (PR)
    # ------------------------------------------------------------
    if "PR" in unique_categories:
        personal_risk_geneset_file = catalogs_cfg.personal_risk_geneset
        output_dir = Path(personal_risk_geneset_file).resolve().parent
        bed_file = output_dir / f"PR_risk_genes_{assembly}.bed"
        json_file = output_dir / f"PR_risk_genes.json"
        if not Path(bed_file).exists():
            build_json_bed_files(
                "PR",
                assembly,
                personal_risk_geneset_file,
                bed_file,
                json_file,
                vcf_file,
            )
        # Store BED and JSON files in ctx object
        ctx.outputs["catalogs"]["bed_files"]["PR"] = str(bed_file)
        ctx.outputs["catalogs"]["json_files"]["PR"] = str(json_file)

    # ------------------------------------------------------------
    # Reproductive Risk (RR)
    # ------------------------------------------------------------
    if "RR" in unique_categories:
        reproductive_risk_geneset_file = catalogs_cfg.reproductive_risk_geneset
        output_dir = Path(reproductive_risk_geneset_file).resolve().parent
        bed_file = output_dir / f"RR_risk_genes_{assembly}.bed"
        json_file = output_dir / f"RR_risk_genes.json"
        if not Path(bed_file).exists():
            build_json_bed_files(
                "RR",
                assembly,
                reproductive_risk_geneset_file,
                bed_file,
                json_file,
                vcf_file,
            )
        # Store BED and JSON files in ctx object
        ctx.outputs["catalogs"]["bed_files"]["RR"] = str(bed_file)
        ctx.outputs["catalogs"]["json_files"]["RR"] = str(json_file)

