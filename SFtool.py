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

import sys
import pickle
from pathlib import Path

from sftool.utils.errors import (
    BootstrapError
)

from sftool.utils.arguments import parse_arguments
from sftool.core.bootstrap import bootstrap_execution
from sftool.steps.catalog_generation import run as run_catalog_generation
from sftool.steps.clinvar_setup import run as run_clinvar_setup
from sftool.steps.sample_preprocessing import run as run_sample_preprocessing
from sftool.steps.variant_evidence_preparation import run as run_variant_evidence_preparation
from sftool.steps.variant_collection import run as run_variant_collection
from sftool.steps.variant_selection import run as run_variant_selection
from sftool.steps.report_generation import run as run_report_generation


def main():

    # --------------------------
    # Parse CLI arguments
    # --------------------------
    args = parse_arguments()
    outdir = args.outdir

    if not args.debug_load_ctx:

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
        #           Only if 'clinvar' is present in variant_classification_sources
        # ----------------------------
        variant_classification_sources = ctx.variant_classification_sources
        # Get unique list of categories for all samples
        categories = sorted({
            c for s in ctx.samples for c in s.categories
        })

        if 'clinvar' in variant_classification_sources and ("PR" in categories or "RR" in categories):
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
        #
        # ----------------------------
        run_variant_selection(ctx)

        if args.debug_dump_ctx:
            debug_ctx_path = Path(outdir) / "ctx_backup_rr_screening_CFTR.pkl"
            with open(debug_ctx_path, 'wb') as f:
                pickle.dump(ctx, f)

    else:
        debug_ctx_path = Path(outdir) / "ctx_backup_rr_screening_CFTR.pkl"
        with open(debug_ctx_path, "rb") as f:
            ctx = pickle.load(f)
        # ----------------------------
        # STEP 6: REPORT GENERATION
        #
        # ----------------------------
        run_report_generation(ctx)





if __name__ == "__main__":
    main()
