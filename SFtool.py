#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Herramienta para el manejo automático de hallazgos secundarios.

Esta herramienta permite a los usuarios analizar archivos VCF para el manejo automático de hallazgos secundarios relacionados con riesgo personal, riesgo reproductivo y farmacogenético.

@Dependencies InterVar and AnnoVar 
@Usage python3 SFtool.py input_file.vcf --mode <Option: 'basic' or 'advanced'> --evidence <integer> --assembly <Option: '37' or '38'>
@Arguments:
    -vcf (str): Ruta al archivo VCF de entrada.
    -outpath (str): Ruta al directorio donde se guardarán los resultados.
    -mode (str): Modo de análisis (básico o avanzado).
    -evidence (int): Nivel de evidencia de ClinVar para el modo avanzado (1-4).#comprobar que lo he puesto de memoria
    -assembly (str): Ensamblaje genómico a utilizar (GRCh37 o GRCh38).

@Author Javier Perez FLorido, Edurne Urrutia Lafuente
@Date 2023/08/01
@email javier.perez.florido.sspa@juntadeandalucia.es, edurlaf@gmail.com
@github https://github.com/babelomics/secondaryfindings
"""

import os
import sys

from modules.misc.errors import (
    BootstrapError
)

from modules.misc.arguments import parse_arguments
from modules.bootstrap import bootstrap_execution
from modules.FG.run_fg_module import run_pharmacogenomic_risk_module
from modules.misc.report_utils import generate_report
from steps.catalog_generation import run as run_catalog_generation
from steps.clinvar_setup import run as run_clinvar_setup
from steps.sample_preprocessing import run as run_sample_preprocessing
from steps.variant_evidence_preparation import run as run_variant_evidence_preparation
from steps.variant_collection import run as run_variant_collection
from steps.variant_selection import run as run_variant_selection


def main():

    # --------------------------
    # Parse CLI arguments
    # --------------------------
    args = parse_arguments()
    outdir = args.outdir

    # ----------------------------
    # STEP 1: SFtool bootstraping: JSON inputs validation, create execution context object and validate run dependencies
    # ----------------------------

    try:
        ctx = bootstrap_execution(args.samples, args.config, outdir)
    except BootstrapError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    # ----------------------------
    # STEP 2: JSON and BED files catalog generation
    # ----------------------------
    run_catalog_generation(ctx)

    # ----------------------------
    # STEP 3: CLINVAR DDBB MANAGEMENT
    #           Only if profile is advanced and for PR and RR categories
    # ----------------------------
    profile = ctx.profile
    # Get unique list of categories for all samples
    categories = sorted({
        c for s in ctx.samples for c in s.categories
    })

    if profile == 'advanced' and ("PR" in categories or "RR" in categories):
        run_clinvar_setup(ctx)


    # ----------------------------
    # STEP 4: SAMPLE PREPROCESSING
    #
    # ----------------------------
    run_sample_preprocessing(ctx)


    # ----------------------------
    # STEP 5: VARIANT EVIDENCE PREPARATION (GENEBE AND/OR CLINVAR)
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
    #           Only for PR and RR categories
    # ----------------------------
    if "PR" in categories or "RR" in categories:
        run_variant_selection(ctx)



    # ------------------------------------------------------------
    # Sample-level (same semantics as samples_data["samples"][0])
    # ------------------------------------------------------------
    sample = ctx.samples[0]
    vcf_file = str(sample.vcf)
    assembly = ctx.assembly

    """
    Execute modules selected by the user according to categories. For reproductive risk (rr) category, take into account results from SMAca or STRipy if available
    """
    # Run modules selected by user

    pr_results = None
    rr_results = None

    haplot_results = None
    pharmCAT_report_file = None



    if "PGx" in categories: # Run Pharmacogenetic (FG) module - pharmCAT
        if assembly == "GRCh38": # pharmCAT is only allowed for GRCh38 assembly
            python_path = ctx.config.paths.python
            pharmCAT_path = ctx.config.paths.pharmCAT
            htslib_path = ctx.config.paths.htslib
            java_path = ctx.config.paths.java
            bcftools_path = ctx.config.paths.bcftools
            out_path = ctx.base_output_dir
            [pharmCAT_report_file, haplot_results] = run_pharmacogenomic_risk_module(vcf_file, python_path, pharmCAT_path, bcftools_path, htslib_path, java_path, out_path)
        else:
            print("Farmacogenomic module (pharmCAT) is available only for GRCh38 human assembly")

    """
    Create report
    """
    out_path = outdir
    generate_report(pr_results, rr_results, haplot_results, pharmCAT_report_file, ctx.config, args, ctx.outputs["clinvar"]["clinvar_db"], categories, out_path)


    
        
if __name__ == "__main__":
    main()
