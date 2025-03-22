

import os
import configparser

def check_dependencies(genebe_path, bcftools_path, java_path, pharmCAT_path, python_path):

    print("Checking dependencies...")

    # Check whether GeneBe client is installed
    if not os.path.exists(os.path.join(genebe_path, "GeneBeClient.jar")):
        print("GeneBe is not installed. Exiting")
        exit(1)

    # Check whether bcftools is installed
    if not os.path.exists(os.path.join(bcftools_path, "bcftools")):
        print("bcftools is not installed. Exiting")
        exit(1)

    # Check whether Java is installed
    if not os.path.exists(java_path):
        print("Java is not installed. Exiting")
        exit(1)

    # Check whether pharmCAT is installed
    if not os.path.exists(os.path.join(pharmCAT_path, "pharmcat.jar")):
        print("pharmCAT is not installed. Exiting")
        exit(1)

    # Check whether python3 is installed
    if not os.path.exists(python_path):
        print("Python3 is not installed. Exiting")
        exit(1)

