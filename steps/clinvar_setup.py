# steps/clinvar_setup.py

from modules.context import ExecutionContext
from modules.misc.clinvar_utils import clinvar_manager

def run(ctx: ExecutionContext) -> None:
    """
    Prepare ClinVar database files for the selected assembly
    and register them in ctx.outputs["clinvar"].
    """

    [clinvar_db, clinvar_submission] = clinvar_manager(
        ctx.config.clinvar.db_path,
        ctx.config.clinvar.version,
        ctx.assembly,
    )

    ctx.outputs["clinvar"]["clinvar_db"] = clinvar_db
    ctx.outputs["clinvar"]["clinvar_summary_db"] = clinvar_submission
    ctx.outputs["clinvar"]["clinvar_db_version"] = ctx.config.clinvar.version
    ctx.outputs["clinvar"]["clinvar_db_assembly"] = ctx.assembly
