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
    ValidationError,
    RuntimeDependencyError
)

from modules.misc.arguments import parse_arguments
from modules.misc.validation import validate_all

from modules.misc.build_json_bed_files import build_json_bed_files
from modules.misc.clinvar_utils import clinvar_manager
from modules.FG.run_fg_module import run_pharmacogenomic_risk_module
from modules.PR_RR.run_pers_repro_risk_module import run_pers_repro_risk_module
from modules.misc.runtime import check_runtime_dependencies
from modules.misc.vcf_utils import normalize_vcf, intersect_vcf_with_bed
from modules.misc.report_utils import generate_report
from modules.STRipy.parse_STRipy_output import parse_STRipy_output
from modules.SMAca.parse_SMAca_output import parse_SMAca_output

from modules.context import ExecutionContext
from modules.config import (
    CatalogConfig,
    ClinVarConfig,
    GeneBeConfig,
    SMAcaConfig,
    ReferenceDataConfig,
    PathsConfig
)


def main():

    # --------------------------
    # Parse CLI arguments
    # --------------------------
    args = parse_arguments()

    samples_path = args.samples
    config_path = args.config
    outdir = args.outdir

    try:
        # --------------------------
        # Validate JSON inputs
        # --------------------------
        samples_data, config_data = validate_all(samples_path, config_path)

        # --------------------------
        # Prepare output and temporal directory
        # --------------------------
        if os.path.exists(outdir):
            if not args.force:
                print(f"ERROR: Output directory '{outdir}' already exists. Use --force to overwrite.")
                sys.exit(1)
        else:
            os.makedirs(outdir, exist_ok=True)

    except ValidationError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    try:
        # --------------------------
        # Check dependencies
        # --------------------------
        check_runtime_dependencies(config_data)
    except RuntimeDependencyError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    # Create execution context
    ctx = ExecutionContext(samples_data, outdir)

    """
    1. Generate JSON and BED files for PR or RR categories
    """
    catalogs_cfg = CatalogConfig(config_data["catalogs"])
    reference_cfg = ReferenceDataConfig(config_data["references"])

    # ------------------------------------------------------------
    # Sample-level (same semantics as samples_data["samples"][0])
    # ------------------------------------------------------------
    sample = ctx.samples[0]
    categories = sample.categories
    vcf_file = str(sample.vcf)

    # ------------------------------------------------------------
    # Run-level
    # ------------------------------------------------------------
    assembly = ctx.assembly
    paths_cfg = PathsConfig(config_data["paths"])

    categories_path = paths_cfg.categories

    # ------------------------------------------------------------
    # Reference genome (via ReferenceDataConfig wrapper)
    # ------------------------------------------------------------
    reference_genome = reference_cfg.cfg["genomes"][assembly]

    # ------------------------------------------------------------
    # Personal Risk (PR)
    # ------------------------------------------------------------
    if "PR" in categories:
        personal_risk_geneset_file = (
            catalogs_cfg.cfg["personal_risk_geneset"]
        )

        bed_path = f"{categories_path}/PR/PR_risk_genes_{assembly}.bed"
        if not os.path.exists(bed_path):
            build_json_bed_files(
                "pr",
                assembly,
                categories_path,
                personal_risk_geneset_file,
                vcf_file,
            )

    # ------------------------------------------------------------
    # Reproductive Risk (RR)
    # ------------------------------------------------------------
    if "RR" in categories:
        reproductive_risk_geneset_file = (
            catalogs_cfg.cfg["reproductive_risk_geneset"]
        )

        bed_path = f"{categories_path}/RR/RR_risk_genes_{assembly}.bed"
        if not os.path.exists(bed_path):
            build_json_bed_files(
                "rr",
                assembly,
                categories_path,
                reproductive_risk_geneset_file,
                vcf_file,
            )


    """
    In advanced mode, check/update clinVar database
    """
    # If "advanced" mode, check whether Clinvar Database exists
    profile = ctx.profile
    clinvar_cfg = ClinVarConfig(config_data["clinvar"])

    if profile == 'advanced' and ("PR" in categories or "RR" in categories):
        [clinvar_db, clinvar_submission] = clinvar_manager(
            clinvar_cfg.db_path,
            clinvar_cfg.version,
            assembly,
        )
    else:
        clinvar_db = None
        clinvar_submission = None

    """
    VCF normalization: only of PR or RR cateogry (FG has its own normalization procedure)
    """

    bcftools_path = paths_cfg.bcftools
    if "PR" in categories or "RR" in categories:
        temp_path = ctx.tmp_dir
        norm_vcf_file = normalize_vcf(vcf_file, temp_path, bcftools_path, reference_genome)

    """
    Normalized VCF and BED intersection for each category
    """
    input_vcf_files = {}
    for category in categories:
        if category == "PR" or category == "RR":
            category_bed_file = os.path.join(categories_path + category.upper(), category + '_risk_genes_' + assembly + '.bed')
            generated_vcf_file = intersect_vcf_with_bed(norm_vcf_file, category_bed_file, temp_path, category)
            input_vcf_files[category] = generated_vcf_file

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
    genebe_path = paths_cfg.genebe
    java_path = paths_cfg.java

    genebe_cfg = GeneBeConfig(config_data["genebe_credentials"])
    genebe_apikey = genebe_cfg.api_key
    genebe_username = genebe_cfg.username

    clinvar_evidence = ctx.clinvar_evidence

    if "PR" in categories:
        # Run Personal Risk (PR) module. P/LP variants from GeneBe and/or CLINVAR in Genes related to pr category
        pr_results = run_pers_repro_risk_module(input_vcf_files['PR'], assembly, profile, clinvar_evidence, clinvar_db, clinvar_submission, 'pr', personal_risk_geneset_file, genebe_path, java_path, genebe_apikey, genebe_username)
    if "RR" in categories:
        # Run Reproductive Risk (RR) module. P/LP variants from GeneBe and/or CLINVAR in Genes related to rr category
        rr_results = run_pers_repro_risk_module(input_vcf_files['RR'], assembly, profile, clinvar_evidence, clinvar_db, clinvar_submission, 'rr', reproductive_risk_geneset_file, genebe_path, java_path, genebe_apikey, genebe_username)
        # Parse STRipy JSON file (if provided)
        STRipy_output = samples_data.get("samples")[0].get("stripy_path")
        if STRipy_output != "None":
            reproductive_risk_geneset_STR_file = (
                catalogs_cfg.cfg["reproductive_risk_geneset_STR"]
            )
            STRipy_results_rr = parse_STRipy_output(reproductive_risk_geneset_STR_file, STRipy_output)
        # Parse SMAca CSV file (if provided)
        SMAca_output = sample.smaca_path
        if SMAca_output != "None" and "RR" in categories:
            smaca_cfg = SMAcaConfig(config_data.get("smaca_thresholds", {}))
            smaca_cv_fail_threshold = smaca_cfg.cv_fail
            smaca_cv_warn_threshold = smaca_cfg.cv_warn
            smaca_low_cov_abs = smaca_cfg.low_cov_absolute
            smaca_low_cov_rel = smaca_cfg.low_cov_relative
            SMAca_results_rr = parse_SMAca_output(SMAca_output, smaca_cv_fail_threshold, smaca_cv_warn_threshold, smaca_low_cov_abs, smaca_low_cov_rel)
    if "PGx" in categories: # Run Pharmacogenetic (FG) module - pharmCAT
        if assembly == "GRCh38": # pharmCAT is only allowed for GRCh38 assembly
            python_path = config_data.get("paths").get("python")
            pharmCAT_path = paths_cfg.pharmCAT
            htslib_path = paths_cfg.htslib
            java_path = paths_cfg.java
            bcftools_path = paths_cfg.bcftools
            out_path = outdir
            [pharmCAT_report_file, haplot_results] = run_pharmacogenomic_risk_module(vcf_file, python_path, pharmCAT_path, bcftools_path, htslib_path, java_path, out_path)
        else:
            print("Farmacogenomic module (pharmCAT) is available only for GRCh38 human assembly")

    """
    Create report
    """
    out_path = outdir
    generate_report(pr_results, rr_results, haplot_results, pharmCAT_report_file, config_data, args, clinvar_db, categories, out_path)


    
        
if __name__ == "__main__":
    main()
