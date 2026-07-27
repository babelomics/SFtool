"""
SFtool resource-management commands.

This package defines the ``sftool resources`` command group and its
subcommands.
"""

import click

from sftool.commands.resources.setup import setup


RESOURCES_HELP = """
Manage SFtool datasets and shared resources.

Resources must be prepared before running an SFtool analysis.

\b
Available operations:
  setup   Download and generate the resources required by SFtool.

\b
Example:
  sftool resources setup \\
    --output-dir /data/sftool_resources
"""


@click.group(
    help=RESOURCES_HELP,
)
def resources() -> None:
    """
    Manage SFtool datasets and shared resources.
    """


resources.add_command(setup)