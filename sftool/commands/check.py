import click


CHECK_HELP = """
Validate the SFtool execution environment and configuration.

\b
Examples:
  sftool check --config config.json 

  sftool check --config config.json --samples samples_info.json
"""

@click.command(
    help=CHECK_HELP
)
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the SFtool configuration JSON file (see template in examples/config.schema.json)."
)
@click.option(
    "--samples",
    "samples_path",
    required=False,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Optional samples_info JSON file to validate sample-specific inputs (see template in examples/samples_info.schema.json)."
)
def check(config_path, samples_path):
    """
    Check SFtool environment and input files.
    """
    # TO BE DONE