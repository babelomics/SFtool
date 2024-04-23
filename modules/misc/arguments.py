# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 01:04:02 2023

@author: kindi
"""

import argparse
import sys

def arguments():
    """
    Get the arguments
    """
    
    parser = argparse.ArgumentParser(description="Secondary findings analysis tool. \n")
    
    # Argumento para el archivo VCF
    parser.add_argument('vcf_file', metavar='VCF_FILE', type=str, help='VCF input file (vcf.gz format)')
    
    # Argumento para el modo de análisis (básico o avanzado)
    parser.add_argument('--mode', choices=['basic', 'advanced'], default='basic', help='Analysis mode (basic or advanced). Basic by default')
    
    # Argumento para el nivel de evidencia (solo en modo avanzado)
    parser.add_argument('--evidence', type=int, choices=range(1, 5), default=1, help='Evidence level (1-4). Only in advanced mode. 1 by default')
    
    # Argumento para el genoma de referencia
    parser.add_argument('--assembly', type=int, choices=[37, 38], default=37, help='Reference genome. GRCh37 by default')
    
    # Argumento para el archivo de texto de HPOs
    parser.add_argument("--hpos_file", default=None, help="TXT file with a list of HPOs for the patient (optional)")

    # Argument for config file
    parser.add_argument("--config_file", default=None, help="Configuration file", required=True)

    # Categories to be run: PR, RR, FG (by default all categories are selected)
    parser.add_argument("--categories", type=str, default="PR,RR,FG", help="PR,RR,FG or any combination of those categories or a single one. All categories by default")
    
    try:
        args = parser.parse_args()
        tmp_categories = args.categories.replace(" " ,"").strip() # Remove spaces
        tmp_categories = tmp_categories.split(",")  # To a list

        if (len(tmp_categories) == 3 and 'PR' in tmp_categories and 'RR' in tmp_categories and 'FG' in tmp_categories) or \
                (len(tmp_categories) == 2 and (('PR' in tmp_categories and 'RR' in tmp_categories) or ('PR' in tmp_categories and 'FG' in tmp_categories) or ('RR' in tmp_categories and 'FG' in tmp_categories))) or \
                (len(tmp_categories) == 1 and ('PR' in tmp_categories or 'RR' in tmp_categories or 'FG' in tmp_categories)):
            return [args, list(map(lambda x: x.lower(), tmp_categories))]
        else:
            raise argparse.ArgumentError(args.categories, "Should be a comma separated string with any (or all) of the following categories: PR,RR,FG")
    
    except:
        print("\nPlease, introduce required arguments")
        parser.print_help(sys.stderr)
        sys.exit()




