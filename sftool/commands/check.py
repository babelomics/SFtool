import click


@click.command(
    help=(
            "Check whether SFtool inputs, configuration files and external "
            "dependencies are correctly available."
    )
)
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the SFtool configuration JSON file."
)
@click.option(
    "--samples",
    "samples_path",
    required=False,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Optional samples_info JSON file to validate sample-specific inputs."
)
def check(config_path, samples_path):
    """
    Check SFtool environment and input files.
    """
    # TO BE DONE