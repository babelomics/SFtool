import click
import sys
import pickle
from pathlib import Path

from sftool.utils.errors import (
    BootstrapError
)

from sftool.core.bootstrap import bootstrap_execution
from sftool.steps.catalog_selection import run as run_catalog_selection
from sftool.steps.sample_preprocessing import run as run_sample_preprocessing
from sftool.steps.variant_evidence_preparation import run as run_variant_evidence_preparation
from sftool.steps.variant_collection import run as run_variant_collection
from sftool.steps.variant_selection import run as run_variant_selection
from sftool.steps.report_generation import run as run_report_generation
from sftool.steps.variant_confirmation import run as run_variant_confirmation


RUN_HELP = """
Run the complete SFtool workflow.

\b
Example:
  sftool run \\
    --samples samples_info.json \\
    --config config.json \\
    --outdir results \\
    --force
"""

@click.command(
    help=RUN_HELP
)
@click.option(
    "--samples",
    "samples_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the samples_info JSON file (see template in examples/samples_info.schema.json)."
)
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the SFtool configuration JSON file (see template in examples/config.schema.json)."
)
@click.option(
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, writable=True),
    help="Directory where SFtool output will be written."
)
@click.option(
    "--force",
    is_flag=True,
    help="Overwrite or reuse an existing output directory when applicable."
)
@click.option(
    "--debug-dump-ctx",
    is_flag=True,
    hidden=True,
    help="Dump the internal ExecutionContext object for debugging."
)
@click.option(
    "--debug-load-ctx",
    is_flag=True,
    hidden=True,
    help="Load a previously dumped ExecutionContext object for debugging."
)

def run(samples_path, config_path, outdir, force, debug_dump_ctx, debug_load_ctx):
    """
    Run SFtool analysis.
    """

    if not debug_load_ctx:

        # ----------------------------
        # STEP 1: SFtool bootstraping: JSON inputs validation, create execution context object and validate run dependencies
        # ----------------------------

        try:
            ctx = bootstrap_execution(samples_path, config_path, outdir)
        except BootstrapError as e:
            print(f"[ERROR] {e}")
            sys.exit(1)

        # ----------------------------
        # STEP 2: JSON and BED files catalog selection
        # ----------------------------
        run_catalog_selection(ctx)

        # ----------------------------
        # STEP 3: SAMPLE PREPROCESSING
        #
        # ----------------------------
        run_sample_preprocessing(ctx)

        # ----------------------------
        # STEP 4: VARIANT CONFIRMATION
        #
        # ----------------------------
        if "variant_confirmation" in ctx.modes:
            run_variant_confirmation(ctx)


        # ----------------------------
        # STEP 5: VARIANT EVIDENCE PREPARATION (GENEBE AND PharmCAT)
        #
        # ----------------------------
        run_variant_evidence_preparation(ctx)



        # ----------------------------
        # STEP 6: VARIANT COLLECTION (GENEBE AND/OR CLINVAR, STRs and SMN1-copy, pharmCAT)
        #
        # ----------------------------
        run_variant_collection(ctx)


        # ----------------------------
        # STEP 7: VARIANT SELECTION from the set of VARIANT COLLECTION
        #
        # ----------------------------
        run_variant_selection(ctx)

        if debug_dump_ctx:
            debug_ctx_path = Path(outdir) / "ctx_backup_RR_CFTR_couple_screening.pkl"
            with open(debug_ctx_path, 'wb') as f:
                pickle.dump(ctx, f)

    else:
        debug_ctx_path = Path(outdir) / "ctx_backup_RR_CFTR_couple_screening.pkl"
        with open(debug_ctx_path, "rb") as f:
            ctx = pickle.load(f)
        # ----------------------------
        # STEP 8: REPORT GENERATION
        #
        # ----------------------------
        run_report_generation(ctx)
