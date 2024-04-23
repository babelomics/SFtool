

import os
import configparser

def check_dependencies(intervar_path, bcftools_path):

    print("Checking dependencies...")

    # Check whether InterVar is installed
    if not os.path.exists(os.path.join(intervar_path, "Intervar.py")):
        print("InterVar is not installed. Exiting")
        exit(1)

    # Check whether bcftools is installed
    if not os.path.exists(os.path.join(bcftools_path, "bcftools")):
        print("bcftools is not installed. Exiting")
        exit(1)

    # Check whether annovar is installed. Annovar must be installed as a dependency of Intervar, and annovar paths must added to Intervar's config.ini file through the variables
    # convert2annovar, table_annovar and annotate_variation

    config = configparser.ConfigParser()

    config.read(os.path.join(intervar_path, "config.ini"))

    if not os.path.exists(config.get('Annovar', 'convert2annovar')) or not os.path.exists(config.get('Annovar','table_annovar')) or not os.path.exists(config.get('Annovar','annotate_variation')):
        print("Annovar files (convert2annovar.pl or table_annovar.pl or annotate_variation.pl) do not exist. Exiting")
        exit(1)

