#!/usr/bin/python3
# -*- coding: utf-8 -*-

import click

from sftool.commands.run import run
from sftool.commands.check import check

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"]
}

@click.group(
    context_settings=CONTEXT_SETTINGS,
    invoke_without_command=True
)
@click.version_option()
@click.pass_context

def main(ctx):
    """
    SFtool: manage personal, reproductive and pharmacogenomic secondary findings.
    """
    if ctx.invoked_subcommand is None:
        click.secho("Error: no command specified.\n", fg="red", bold=True)
        click.echo("Available commands:")
        click.echo("  run     Execute the complete SFtool workflow.")
        click.echo("  check   Validate the execution environment and configuration.")
        click.echo()
        click.echo("Run 'sftool --help' for more information.")
        ctx.exit(1)


main.add_command(run)
main.add_command(check)

if __name__ == "__main__":
    main()
