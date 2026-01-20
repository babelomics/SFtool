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
from modules.PR_RR.run_pers_repro_risk_module import run_pers_repro_risk_module
from modules.misc.vcf_utils import normalize_vcf, intersect_vcf_with_bed
from modules.misc.report_utils import generate_report
from modules.STRipy.STRipy_collection import STRipy_collection
from modules.SMAca.parse_SMAca_output import parse_SMAca_output

from steps.catalog_generation import run as run_catalog_generation
from steps.clinvar_setup import run as run_clinvar_setup
from steps.sample_preprocessing import run as run_sample_preprocessing
from steps.variant_evidence_preparation import run as run_variant_evidence_preparation


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
    #           Only for PR and RR categories
    # ----------------------------
    if "PR" in categories or "RR" in categories:
        run_sample_preprocessing(ctx)


    # ----------------------------
    # STEP 5: VARIANT EVIDENCE PREPARATION (GENEBE AND/OR CLINVAR)
    #           Only for PR and RR categories
    # ----------------------------
    if "PR" in categories or "RR" in categories:
        run_variant_evidence_preparation(ctx)


    catalogs_cfg = ctx.config.catalogs

    # ------------------------------------------------------------
    # Sample-level (same semantics as samples_data["samples"][0])
    # ------------------------------------------------------------
    sample = ctx.samples[0]
    vcf_file = str(sample.vcf)

    # ------------------------------------------------------------
    # Run-level
    # ------------------------------------------------------------
    assembly = ctx.assembly

    """
    Execute modules selected by the user according to categories. For reproductive risk (rr) category, take into account results from SMAca or STRipy if available
    """
    # Run modules selected by user

    pr_results = None
    rr_results = None
    STRipy_results_rr = None
    SMAca_results_rr = None
    haplot_results = None
    pharmCAT_report_file = None
    genebe_path = ctx.config.paths.genebe
    java_path = ctx.config.paths.java


    genebe_apikey = ctx.config.genebe_credentials.api_key
    genebe_username = ctx.config.genebe_credentials.username

    clinvar_evidence = ctx.clinvar_evidence

    if "PR" in categories:
        # Run Personal Risk (PR) module. P/LP variants from GeneBe and/or CLINVAR in Genes related to pr category
        pr_results = run_pers_repro_risk_module(sample.vcf_outputs["intersected"]['PR'], assembly, profile, clinvar_evidence, ctx.outputs["clinvar"]["clinvar_db"], ctx.outputs["clinvar"]["clinvar_summary_db"], 'pr', ctx.config.catalogs.personal_risk_geneset, genebe_path, java_path, genebe_apikey, genebe_username)
    if "RR" in categories:
        # Run Reproductive Risk (RR) module. P/LP variants from GeneBe and/or CLINVAR in Genes related to rr category
        rr_results = run_pers_repro_risk_module(sample.vcf_outputs["intersected"]['RR'], assembly, profile, clinvar_evidence, ctx.outputs["clinvar"]["clinvar_db"], ctx.outputs["clinvar"]["clinvar_summary_db"], 'rr', ctx.config.catalogs.reproductive_risk_geneset, genebe_path, java_path, genebe_apikey, genebe_username)
        # Parse STRipy JSON file (if provided)
        STRipy_output = ctx.samples[0].stripy_path
        if STRipy_output != "None":
            reproductive_risk_geneset_STR_file = (
                catalogs_cfg.reproductive_risk_geneset_STR
            )
            STRipy_results_rr = STRipy_collection(reproductive_risk_geneset_STR_file, STRipy_output)
        # Parse SMAca CSV file (if provided)
        SMAca_output = sample.smaca_path
        if SMAca_output != "None" and "RR" in categories:

            smaca_cv_fail_threshold = ctx.config.smaca_thresholds.cv_fail
            smaca_cv_warn_threshold = ctx.config.smaca_thresholds.cv_warn
            smaca_low_cov_abs = ctx.config.smaca_thresholds.low_cov_absolute
            smaca_low_cov_rel = ctx.config.smaca_thresholds.low_cov_relative
            SMAca_results_rr = parse_SMAca_output(SMAca_output, smaca_cv_fail_threshold, smaca_cv_warn_threshold, smaca_low_cov_abs, smaca_low_cov_rel)
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
