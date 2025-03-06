

import os
import configparser

def check_dependencies(genebe_path, bcftools_path):

    print("Checking dependencies...")

    # Check whether InterVar is installed
    if not os.path.exists(os.path.join(genebe_path, "GeneBeClient.jar")):
        print("GeneBe is not installed. Exiting")
        exit(1)

    # Check whether bcftools is installed
    if not os.path.exists(os.path.join(bcftools_path, "bcftools")):
        print("bcftools is not installed. Exiting")
        exit(1)

