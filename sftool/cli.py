#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Herramienta para el manejo automático de hallazgos secundarios.

Esta herramienta permite a los usuarios analizar archivos VCF para el manejo automático de hallazgos secundarios relacionados con riesgo personal, riesgo reproductivo y farmacogenético.


@Usage python3.10 -m sftool.cli
    --samples <samples_info.json>
    --config <config_info.json>
     --outdir <root_output_dir>
     --force

@Author Javier Perez FLorido, Edurne Urrutia Lafuente
@Date 2023/08/01
@email javier.perez.florido.sspa@juntadeandalucia.es, edurlaf@gmail.com
@github https://github.com/babelomics/SFtool
"""

import click

from sftool.commands.run import run
from sftool.commands.check import check


@click.group(
    context_settings={"help_option_names": ["-h", "--help"]}
)
@click.version_option()
def main():
    """
    SFtool: manage personal, reproductive and pharmacogenomic secondary findings.
    """
    pass


main.add_command(run)
main.add_command(check)

if __name__ == "__main__":
    main()
