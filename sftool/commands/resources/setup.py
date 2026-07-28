"""
Implementation of ``sftool resources setup``.
"""

from pathlib import Path
from sftool.utils.catalog_utils import (
    CatalogGenerationError,
    prepare_catalog_resources,
)

from sftool.utils.clinvar_utils import (
    ClinVarProcessingError,
    build_clinvar_databases,
    download_bundled_clinvar_snapshot,
)

from sftool.utils.resource_utils import (
    load_bundled_resources,
    ResourceOperationError,
    ResourceSpecificationError,
    download_versioned_resource
)

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
    type=click.IntRange(min=1, max=4),
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
    help=(
            "Download the GRCh37 and GRCh38 reference genomes."
            "Reserved for future reference genome download support. "
            "Currently, users must provide their own reference genome "
            "through the SFtool configuration."
    ),
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

    if download_reference_genomes:
        click.echo(
            "Reference genome download is not implemented yet. "
            "Users must currently provide their own reference genome "
            "through the SFtool configuration. "
            "This option is reserved for a future release.",
            err=True,
        )


    resource_version = resource_version.lower()

    if resource_version != "bundled":
        raise click.ClickException(
            f"Only bundled versions are implemented at the moment"
        )

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
    click.echo("Downloading bundled ClinVar snapshot...")

    try:
        clinvar_resources = download_bundled_clinvar_snapshot(
            output_root=output_dir,
            clinvar_specification=bundled_resources["clinvar"],
        )
    except (
            ResourceOperationError,
            ResourceSpecificationError,
    ) as error:
        raise click.ClickException(str(error)) from error

    click.echo(
        "  Version: "
        f"{clinvar_resources['version']}"
    )
    click.echo(
        "  Variant summary: "
        f"{clinvar_resources['variant_summary']}"
    )
    click.echo(
        "  Submission summary: "
        f"{clinvar_resources['submission_summary']}"
    )
    click.echo(
        "ClinVar source snapshot downloaded successfully."
    )

    click.echo()
    click.echo("Generating ClinVar assembly databases...")

    try:
        clinvar_databases = build_clinvar_databases(
            variant_summary_path=Path(
                clinvar_resources["variant_summary"]
            ),
            output_root=output_dir,
            version=clinvar_resources["version"],
        )
    except ClinVarProcessingError as error:
        raise click.ClickException(str(error)) from error

    for assembly, database_path in clinvar_databases.items():
        click.echo(
            f"  {assembly}: {database_path}"
        )

    click.echo(
        "ClinVar assembly databases generated successfully."
    )


    click.echo()
    click.echo("Preparing catalog resources...")

    try:
        catalog_resources = prepare_catalog_resources(
            output_root=output_dir,
        )
    except (
            CatalogGenerationError,
            ResourceOperationError,
    ) as error:
        raise click.ClickException(str(error)) from error

    for assembly, catalogs in (
            catalog_resources["assemblies"].items()
    ):
        click.echo(f"  {assembly}")

        for category, resources in catalogs.items():
            click.echo(
                f"    {category}: "
                f"{resources['bed']}, "
                f"{resources['chr_bed']}, "
                f"{resources['json']}"
            )

    click.echo(
        "  RR_STR: "
        f"{catalog_resources['RR_STR']['csv']}"
    )

    click.echo("Catalog resources prepared successfully.")

    click.echo()
    click.echo("Downloading HPO resource...")

    try:
        hpo_resource = download_versioned_resource(
            output_root=output_dir,
            resource_name="hpo",
            specification=bundled_resources["hpo"],
            validator=validate_hpo_gene_to_phenotype,
        )
    except (
            ResourceOperationError,
            ResourceSpecificationError,
    ) as error:
        raise click.ClickException(str(error)) from error

    click.echo(f"  Version: {hpo_resource['version']}")
    click.echo(f"  File: {hpo_resource['path']}")
    click.echo("HPO resource downloaded successfully.")

    click.echo()
    click.echo("Downloading PharmCAT positions resource...")

    try:
        pharmcat_resource = download_versioned_resource(
            output_root=output_dir,
            resource_name="pharmcat",
            specification=bundled_resources["pharmcat"],
            validator=validate_vcf_resource,
        )
    except (
            ResourceOperationError,
            ResourceSpecificationError,
    ) as error:
        raise click.ClickException(str(error)) from error

    click.echo(f"  Version: {pharmcat_resource['version']}")
    click.echo(f"  File: {pharmcat_resource['path']}")
    click.echo("PharmCAT positions resource downloaded successfully.")