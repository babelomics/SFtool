# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 01:04:02 2023

@author: jpflorido
"""

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="SFtool: Secondary Findings Analysis Tool"
    )

    # Mandatory JSON input files
    parser.add_argument(
        "--samples",
        required=True,
        help="Path to samples_info.json"
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to config.json"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for SFtool results"
    )

    # Optional runtime option: overwrite output directory
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output directory if it already exists"
    )

    return parser.parse_args()
