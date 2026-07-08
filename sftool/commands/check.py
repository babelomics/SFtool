from pathlib import Path
from sftool.core.bootstrap import (
    load_json,
    validate_samples_info,
    validate_config,
    validate_execution_context,
)
from sftool.core.config import Config
from sftool.core.context import ExecutionContext, SampleContext
from sftool.utils.errors import RuntimeDependencyError, ValidationError
from sftool.utils.runtime import check_runtime_dependencies
import click


def _ok(message: str):
    click.secho(f"✓ {message}", fg="green")


def _fail(message: str):
    click.secho(f"✗ {message}", fg="red")


def _warn(message: str):
    click.secho(f"⚠ {message}", fg="yellow")


def _build_context_for_check(
        samples_info: dict,
        config: Config,
        output_dir: Path,
) -> ExecutionContext:
    ctx = ExecutionContext(
        execution_meta=samples_info["execution"],
        config=config,
        output_dir=output_dir,
    )

    for sample_data in samples_info["samples"]:
        sample_ctx = SampleContext(sample_data=sample_data, exec_ctx=ctx)
        ctx.add_sample(sample_ctx)

    return ctx


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
    Check SFtool environment and input files without running the workflow
    """
    click.echo()
    click.secho("SFtool environment check", bold=True)
    click.echo()

    try:
        config_path = Path(config_path)
        config_data = load_json(config_path)
        _ok(f"Config JSON loaded: {config_path}")

        samples_info = None

        if samples_path:
            samples_path = Path(samples_path)
            samples_info = load_json(samples_path)
            _ok(f"Samples JSON loaded: {samples_path}")

            validate_samples_info(samples_info)
            _ok("samples_info structure is valid")

            validate_config(config_data, samples_info)
            _ok("config structure is valid for requested categories")
        else:
            validate_config(config_data)
            _ok("config structure is valid")

            _warn(
                "No samples file provided. Sample-specific checks and "
                "category-dependent checks will be skipped."
            )

        config = Config(config_data)
        _ok("Config object created")

        # Check versions
        versions = check_runtime_dependencies(config.paths)
        for tool in versions.values():
            _ok(
                f"{tool.name}: {tool.version} "
                f"(minimum {tool.minimum})"
            )

        if samples_info:
            check_output_dir = Path.cwd() / ".sftool_check_tmp"

            ctx = _build_context_for_check(
                samples_info=samples_info,
                config=config,
                output_dir=check_output_dir,
            )

            validate_execution_context(ctx)
            _ok("Execution context is valid")

        click.echo()
        click.secho("SFtool check completed successfully.", fg="green", bold=True)

    except (
            ValidationError,
            RuntimeDependencyError,
            FileNotFoundError,
            ValueError,
            KeyError,
    ) as e:
        click.echo()
        _fail(str(e))
        click.echo()
        click.secho("SFtool check failed.", fg="red", bold=True)
        raise click.Abort()