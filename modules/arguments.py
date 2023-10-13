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
    
    parser = argparse.ArgumentParser(description="Herramienta de análisis de hallazgos secundarios. \n\
                                     \npython3 secondary_findings.py input_file.vcf -o <out-path-dir> --mode <Option: 'basic' or 'advanced'> --evidence <integer> --assembly <Option: '37' or '38'>\n \
                                     \nEXAMPLE: python3 secondary_findings.py example.vcf -o results --mode basic\n")
    
    # Argumento para el archivo VCF
    parser.add_argument('vcf_file', metavar='VCF_FILE', type=str, help='Archivo VCF de entrada')
    
    # Argumento para el arhivo de salida
    parser.add_argument('-o', '--outpath', dest='out_path',
                   action='store', required=False, default='',
                   help='Output path. Default current directory')
    
    # Argumento para el modo de análisis (básico o avanzado)
    parser.add_argument('--mode', choices=['basic', 'advanced'], default='basic', help='Modo de análisis (basic o advanced)')
    
    # Argumento para el nivel de evidencia (solo en modo avanzado)
    parser.add_argument('--evidence', type=int, choices=range(1, 5), default=1, help='Nivel de evidencia (1-4) en modo avanzado')
    
    # Argumento para el genoma de referencia
    parser.add_argument('--assembly', type=int, choices=['37', '38'], default='37', help='Genoma de referencia')
    
    try:
        args = parser.parse_args()      
        return args
    
    except:
        print("\nPor favor, introduzca los argumentos requeridos.")
        sys.exit()


# se podrían manejar por separado los diferentes errores (ver hcatgpo: flinotfounderror, argumenterror, valueerror)

