"""
Implementation of ``sftool resources setup``.
"""

from pathlib import Path
from sftool.utils.resource_utils import load_bundled_resources

import click


SETUP_HELP = """
Prepare the datasets and generated resources required by SFtool.

The command creates a resource directory containing ClinVar databases,
catalog BED/JSON files, HPO data, PharmCAT positions, and optionally
reference genomes.

\b
Example:
  sftool resources setup \\
    --output-dir /data/sftool_resources \\
    --clinvar-evidence 1 \\
    --resource-version bundled \\
    --download-reference-genomes
"""


@click.command(
    help=SETUP_HELP,
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(
        path_type=Path,
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
    help="Directory where SFtool resources will be installed.",
)
@click.option(
    "--clinvar-evidence",
    type=click.IntRange(min=1, max=5),
    default=1,
    show_default=True,
    help=(
            "Minimum ClinVar evidence level used to generate the "
            "processed ClinVar databases."
    ),
)
@click.option(
    "--resource-version",
    type=click.Choice(
        ["bundled", "latest"],
        case_sensitive=False,
    ),
    default="bundled",
    show_default=True,
    help=(
            "Resource version policy: use versions bundled with SFtool "
            "or resolve the latest available versions."
    ),
)
@click.option(
    "--download-reference-genomes",
    is_flag=True,
    default=False,
    help="Download the GRCh37 and GRCh38 reference genomes.",
)
def setup(
        output_dir: Path,
        clinvar_evidence: int,
        resource_version: str,
        download_reference_genomes: bool,
) -> None:
    """
    Prepare SFtool resources.

    Resource installation is implemented incrementally in the following
    resource-management tasks.
    """

    resource_version = resource_version.lower()

    if resource_version == "bundled":
        bundled_resources = load_bundled_resources()

    click.echo(
        "Bundled resource specification loaded "
        f"(schema version "
        f"{bundled_resources['schema_version']})."
    )

    click.echo("SFtool resource setup")
    click.echo(f"Output directory: {output_dir}")
    click.echo(f"ClinVar evidence: {clinvar_evidence}")
    click.echo(f"Resource version: {resource_version}")
    click.echo(
        "Download reference genomes: "
        f"{'yes' if download_reference_genomes else 'no'}"
    )

    click.echo()
    click.echo(
        "Resource preparation is not implemented yet.",
        err=True,
    )