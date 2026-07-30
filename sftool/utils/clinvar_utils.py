# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 15:12:18 2023

@author: Javier Perez Florido, Edurne Urrutia
"""

from __future__ import annotations

import os
import gzip
import csv
import urllib.request
from datetime import datetime
import shutil
from sftool.utils.catalog_utils import read_catalog_csv
from collections import Counter
from pathlib import Path
from typing import Any
import json
from collections.abc import Mapping
from sftool.utils.resource_utils import write_json

from sftool.utils.resource_utils import (
    ResourceOperationError,
    ResourceSpecificationError,
    download_file,
    ensure_directory,
    render_resource_filename,
    resolve_resource_url,
)

SUPPORTED_CLINVAR_ASSEMBLIES = (
    "GRCh37",
    "GRCh38",
)

CLINVAR_DATABASE_COLUMNS = (
    "Type",
    "Name",
    "GeneSymbol",
    "ClinicalSignificance",
    "ClinSigSimple",
    "RS# (dbSNP)",
    "VariationID",
    "PhenotypeIDS",
    "PhenotypeList",
    "Assembly",
    "Chromosome",
    "Start",
    "Stop",
    "ReviewStatus",
    "SubmitterCategories",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
)

SUPPORTED_CLINVAR_CATALOGS = (
    "PR",
    "RR",
)

SUPPORTED_CLINVAR_EVIDENCE_LEVELS = (
    1,
    2,
    3,
    4,
)

ALLOWED_CLINVAR_VARIANT_TYPES = {
    "deletion",
    "duplication",
    "insertion",
    "indel",
    "single nucleotide variant",
    "microsatellite",
    "variation",
}

class ClinVarProcessingError(RuntimeError):
    """Raised when a ClinVar source file cannot be processed."""

def build_clinvar_assembly_database(
        *,
        variant_summary_path: Path,
        output_root: Path,
        assembly: str,
        version: str,
        overwrite: bool = False,
) -> Path:
    """
    Generate the processed ClinVar database for one genome assembly.

    Parameters
    ----------
    variant_summary_path
        Path to the downloaded ClinVar ``variant_summary`` gzip file.
    output_root
        Root directory of the installed SFtool resources.
    assembly
        Genome assembly to retain: ``GRCh37`` or ``GRCh38``.
    version
        ClinVar version recorded in the bundled resource specification.
    overwrite
        Replace an existing generated database when true.

    Returns
    -------
    Path
        Generated assembly-specific ClinVar database.
    """
    variant_summary_path = Path(variant_summary_path)
    output_root = Path(output_root)

    if assembly not in SUPPORTED_CLINVAR_ASSEMBLIES:
        raise ClinVarProcessingError(
            f"Unsupported ClinVar assembly: {assembly}"
        )

    if not variant_summary_path.is_file():
        raise ClinVarProcessingError(
            "ClinVar variant summary does not exist or is not a file: "
            f"{variant_summary_path}"
        )

    output_directory = ensure_directory(
        output_root / "clinvar" / assembly
    )

    output_path = (
            output_directory
            / f"clinvar_database_{assembly}_{version}.txt"
    )

    if output_path.exists() and not overwrite:
        return output_path

    temporary_path = output_path.with_name(
        f"{output_path.name}.tmp"
    )

    try:
        with gzip.open(
                variant_summary_path,
                mode="rt",
                encoding="utf-8",
                newline="",
        ) as source_handle:
            reader = csv.reader(
                source_handle,
                delimiter="\t",
            )

            try:
                header = next(reader)
            except StopIteration as error:
                raise ClinVarProcessingError(
                    "ClinVar variant summary is empty: "
                    f"{variant_summary_path}"
                ) from error

            missing_columns = [
                column
                for column in CLINVAR_DATABASE_COLUMNS
                if column not in header
            ]

            if missing_columns:
                raise ClinVarProcessingError(
                    "ClinVar variant summary is missing required "
                    f"columns: {', '.join(missing_columns)}"
                )

            column_positions = [
                header.index(column)
                for column in CLINVAR_DATABASE_COLUMNS
            ]
            assembly_position = header.index("Assembly")

            with temporary_path.open(
                    mode="w",
                    encoding="utf-8",
                    newline="",
            ) as output_handle:
                writer = csv.writer(
                    output_handle,
                    delimiter="\t",
                    lineterminator="\n",
                )

                writer.writerow(CLINVAR_DATABASE_COLUMNS)

                for line_number, row in enumerate(
                        reader,
                        start=2,
                ):
                    if len(row) != len(header):
                        raise ClinVarProcessingError(
                            "Malformed ClinVar row at line "
                            f"{line_number}: expected {len(header)} "
                            f"columns, found {len(row)}"
                        )

                    if row[assembly_position] != assembly:
                        continue

                    writer.writerow(
                        row[position]
                        for position in column_positions
                    )

        temporary_path.replace(output_path)

    except ClinVarProcessingError:
        temporary_path.unlink(missing_ok=True)
        raise

    except (OSError, EOFError, gzip.BadGzipFile) as error:
        temporary_path.unlink(missing_ok=True)

        raise ClinVarProcessingError(
            "Could not generate ClinVar database for "
            f"{assembly} from {variant_summary_path}"
        ) from error

    return output_path

def build_clinvar_databases(
        *,
        variant_summary_path: Path,
        output_root: Path,
        version: str,
        assemblies: Sequence[str] = SUPPORTED_CLINVAR_ASSEMBLIES,
        overwrite: bool = False,
) -> dict[str, str]:
    """
    Generate the processed ClinVar databases for the requested assemblies.
    """
    databases: dict[str, str] = {}

    for assembly in assemblies:
        database_path = build_clinvar_assembly_database(
            variant_summary_path=variant_summary_path,
            output_root=output_root,
            assembly=assembly,
            version=version,
            overwrite=overwrite,
        )

        databases[assembly] = str(database_path)

    return databases

def validate_clinvar_catalog(category: str) -> str:
    normalized_category = category.upper()

    if normalized_category not in SUPPORTED_CLINVAR_CATALOGS:
        raise ClinVarProcessingError(
            f"Unsupported ClinVar catalog: {category}. "
            "Expected PR or RR."
        )

    return normalized_category


def validate_clinvar_evidence_level(
        evidence_level: int,
) -> int:
    try:
        normalized_level = int(evidence_level)
    except (TypeError, ValueError) as error:
        raise ClinVarProcessingError(
            f"Invalid ClinVar evidence level: {evidence_level}"
        ) from error

    if normalized_level not in SUPPORTED_CLINVAR_EVIDENCE_LEVELS:
        raise ClinVarProcessingError(
            "Unsupported ClinVar evidence level: "
            f"{normalized_level}. Expected one of "
            f"{SUPPORTED_CLINVAR_EVIDENCE_LEVELS}."
        )

    return normalized_level

def load_clinvar_catalog_genes(
        category: str,
        category_geneset_file: Path,
) -> set[str]:
    normalized_category = validate_clinvar_catalog(category)

    category_geneset_file = Path(category_geneset_file)

    if not category_geneset_file.is_file():
        raise ClinVarProcessingError(
            "ClinVar catalog file does not exist: "
            f"{category_geneset_file}"
        )

    try:
        _, genes = read_catalog_csv(
            category_geneset_file,
            normalized_category,
        )
    except Exception as error:
        raise ClinVarProcessingError(
            "Could not read ClinVar catalog "
            f"{normalized_category}: {category_geneset_file}"
        ) from error

    return {
        str(gene).strip()
        for gene in genes
        if str(gene).strip()
    }

def load_clinvar_submission_summaries(
        submission_summary_path: Path,
) -> dict[str, str]:
    submission_summary_path = Path(submission_summary_path)

    if not submission_summary_path.is_file():
        raise ClinVarProcessingError(
            "ClinVar submission summary does not exist: "
            f"{submission_summary_path}"
        )

    clinical_significance_data: dict[str, Counter] = {}

    try:
        with gzip.open(
                submission_summary_path,
                mode="rt",
                encoding="utf-8",
                newline="",
        ) as source_handle:
            header = None

            for line in source_handle:
                if line.startswith("#"):
                    header = (
                        line.lstrip("#")
                        .rstrip("\n")
                        .split("\t")
                    )
                    continue

                if header is None:
                    continue

                row = dict(
                    zip(
                        header,
                        line.rstrip("\n").split("\t"),
                    )
                )

                variation_id = row.get("VariationID", "")
                clinical_significance = (
                    row.get("ClinicalSignificance", "").strip()
                )
                contributes = row.get(
                    "ContributesToAggregateClassification",
                    "",
                )

                if (
                        not variation_id
                        or not clinical_significance
                        or contributes != "yes"
                ):
                    continue

                counter = clinical_significance_data.setdefault(
                    variation_id,
                    Counter(),
                )
                counter[clinical_significance] += 1

    except (
            OSError,
            EOFError,
            gzip.BadGzipFile,
    ) as error:
        raise ClinVarProcessingError(
            "Could not read ClinVar submission summary: "
            f"{submission_summary_path}"
        ) from error

    return {
        variation_id: "; ".join(
            f"{label} ({count})"
            for label, count in counter.items()
        )
        for variation_id, counter
        in clinical_significance_data.items()
    }

def build_clinvar_variant_entry(
        fields: list[str],
        submission_summaries: Mapping[str, str],
) -> tuple[str, dict[str, Any]] | None:
    variant_type = fields[0]

    if variant_type.lower() not in ALLOWED_CLINVAR_VARIANT_TYPES:
        return None

    position_vcf = fields[15]

    try:
        if int(position_vcf) == -1:
            return None
    except ValueError:
        return None

    variation_id = fields[6]

    variant_key = (
        f"{fields[10]}:"
        f"{position_vcf}:"
        f"{fields[16]}:"
        f"{fields[17]}"
    )

    stars = map_review_status(fields[13])

    entry = {
        "VariantName": fields[1],
        "Gene": fields[2],
        "ClinicalSignificance": fields[3],
        "ClinSigSimple": fields[4],
        "rs": (
            f"rs{fields[5]}"
            if fields[5]
            else ""
        ),
        "ReviewStatus": (
            f"({stars}) {fields[13]}"
        ),
        "ClinvarID": variation_id,
        "PhenotypeIDS": fields[7],
        "ClinSigSummary": submission_summaries.get(
            variation_id,
            "",
        ),
    }

    return variant_key, entry

def build_catalog_clinvar_database(
        *,
        clinvar_database_path: Path,
        submission_summaries: Mapping[str, str],
        category: str,
        category_geneset_file: Path,
        assembly: str,
        evidence_level: int,
        output_root: Path,
        overwrite: bool = False,
) -> Path:
    normalized_category = validate_clinvar_catalog(category)
    normalized_evidence = validate_clinvar_evidence_level(
        evidence_level
    )

    if assembly not in SUPPORTED_CLINVAR_ASSEMBLIES:
        raise ClinVarProcessingError(
            f"Unsupported ClinVar assembly: {assembly}"
        )

    clinvar_database_path = Path(clinvar_database_path)

    if not clinvar_database_path.is_file():
        raise ClinVarProcessingError(
            "Processed ClinVar assembly database does not exist: "
            f"{clinvar_database_path}"
        )

    catalog_genes = load_clinvar_catalog_genes(
        normalized_category,
        category_geneset_file,
    )

    output_directory = ensure_directory(
        Path(output_root)
        / "clinvar"
        / assembly
        / normalized_category
    )

    output_path = (
            output_directory
            / (
                f"clinvar_{assembly}_"
                f"{normalized_category}_"
                f"evidence_{normalized_evidence}.json"
            )
    )

    if output_path.exists() and not overwrite:
        return output_path

    variants: dict[str, dict[str, Any]] = {}

    try:
        with clinvar_database_path.open(
                mode="r",
                encoding="utf-8",
                newline="",
        ) as source_handle:
            reader = csv.reader(
                source_handle,
                delimiter="\t",
            )

            try:
                header = next(reader)
            except StopIteration as error:
                raise ClinVarProcessingError(
                    "Processed ClinVar database is empty: "
                    f"{clinvar_database_path}"
                ) from error

            if tuple(header) != CLINVAR_DATABASE_COLUMNS:
                raise ClinVarProcessingError(
                    "Unexpected processed ClinVar database header: "
                    f"{clinvar_database_path}"
                )

            for line_number, fields in enumerate(
                    reader,
                    start=2,
            ):
                if len(fields) != len(CLINVAR_DATABASE_COLUMNS):
                    raise ClinVarProcessingError(
                        "Malformed processed ClinVar row at line "
                        f"{line_number}: {clinvar_database_path}"
                    )

                row_genes = {
                    gene.strip()
                    for gene in fields[2].split(";")
                    if gene.strip()
                }

                if catalog_genes.isdisjoint(row_genes):
                    continue

                stars = map_review_status(fields[13])

                if stars < normalized_evidence:
                    continue

                result = build_clinvar_variant_entry(
                    fields,
                    submission_summaries,
                )

                if result is None:
                    continue

                variant_key, entry = result
                variants[variant_key] = entry

        write_json(
            data=variants,
            destination=output_path,
        )

    except ClinVarProcessingError:
        raise
    except OSError as error:
        raise ClinVarProcessingError(
            "Could not generate catalog ClinVar database: "
            f"{output_path}"
        ) from error

    return output_path

def build_catalog_clinvar_databases(
        *,
        assembly_databases: Mapping[str, str],
        submission_summary_path: Path,
        catalog_files: Mapping[str, Path],
        output_root: Path,
        evidence_levels: Sequence[int] = (
                SUPPORTED_CLINVAR_EVIDENCE_LEVELS
        ),
        overwrite: bool = False,
) -> dict[str, dict[str, dict[str, str]]]:
    missing_assemblies = [
        assembly
        for assembly in SUPPORTED_CLINVAR_ASSEMBLIES
        if assembly not in assembly_databases
    ]

    if missing_assemblies:
        raise ClinVarProcessingError(
            "Missing processed ClinVar databases for: "
            f"{', '.join(missing_assemblies)}"
        )

    missing_catalogs = [
        category
        for category in SUPPORTED_CLINVAR_CATALOGS
        if category not in catalog_files
    ]

    if missing_catalogs:
        raise ClinVarProcessingError(
            "Missing ClinVar catalog files for: "
            f"{', '.join(missing_catalogs)}"
        )

    normalized_evidence_levels = [
        validate_clinvar_evidence_level(level)
        for level in evidence_levels
    ]

    submission_summaries = load_clinvar_submission_summaries(
        submission_summary_path
    )

    generated: dict[
        str,
        dict[str, dict[str, str]],
    ] = {}

    for assembly in SUPPORTED_CLINVAR_ASSEMBLIES:
        generated[assembly] = {}

        for category in SUPPORTED_CLINVAR_CATALOGS:
            generated[assembly][category] = {}

            for evidence_level in normalized_evidence_levels:
                output_path = build_catalog_clinvar_database(
                    clinvar_database_path=Path(
                        assembly_databases[assembly]
                    ),
                    submission_summaries=submission_summaries,
                    category=category,
                    category_geneset_file=Path(
                        catalog_files[category]
                    ),
                    assembly=assembly,
                    evidence_level=evidence_level,
                    output_root=output_root,
                    overwrite=overwrite,
                )

                generated[assembly][category][
                    str(evidence_level)
                ] = str(output_path)

    return generated

def map_review_status(review_status):
    """
    Mapea el estado de revisión (review status) a un nivel de evidencia (evidence level) equivalente.

    Args:
        review_status (str): Estado de revisión en texto.

    Returns:
        int: Nivel de evidencia equivalente en número.
    """
    # Mapear "Review status" a número de estrellas
    mapping = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, conflicting classifications": 1,
        "criteria provided, single submitter": 1,
        "no classification for the single variant": 0,
        "no classifications from unflagged records": 0,
        "no assertion criteria provided": 0,
        "no classification provided": 0
    }
    return mapping.get(review_status.lower(), 0)  # Valor predeterminado es 0 si no se encuentra en el mapeo

def run_clinvar(evidence_level, clinvar_db, clinvar_submission, category, category_geneset_file, output_clinvar_file):
    """
    Run clinvar using the database according to an evidence level

    Args:
        evidence_level (int): Evidence level for variants
        clinvar_db (str): Path to CLINVAR database with variants
        clinvar_submission (str): Path to CLINVAR submission summary
        assembly (str): Reference genome version
        category (str): either pr or rr
        category_geneset_file (str): Path to CSV file for the given category
        output_clinvar_file (str): Output JSON file where clinvar results are saved

    Returns:
        dict:  A path to file that that contains variants from Clinvar and their related information.

    Raises:
        Exception: An exception occurs when an error arises in Clinvar
    """
    try:
        # Select genes for current category
        # Read CSV and store it in a dictionary
        genes_dct, genes_lst = read_catalog_csv(Path(category_geneset_file), category)

        # Read Clinvar database
        clinvar_dct = {}  #  dictionary to store information from CLINVAR

        all_clinvar_id = []

        allowed_types = {"deletion", "duplication", "insertion", "indel", "single nucleotide variant", "microsatellite", "variation"} # Only check for variants that corresponds to SNVs and small indels (these are the variants expected in the VCF file)

        with open(clinvar_db, "r") as db_file:
            for line in db_file:
                line = line.rstrip()
                if line == "":
                    continue
                fields = line.strip().split("\t")
                gene = fields[2]
                pos = fields[15]
                variant_type = fields[0]
                if any(g in genes_lst for g in gene.split(';')) and variant_type.lower() in allowed_types and int(pos) != -1:
                    # Only parse entries for the corresponding category (a given entry in clinvar can contain a set of genes separated by ;), belongs to the allowed types or has a position defined
                    variant = f"{fields[10]}:{pos}:{fields[16]}:{fields[17]}"
                    variant_name = fields[1]
                    clinical_significance = fields[3]
                    clinsigsimple = fields[4]
                    rs_id = fields[5]
                    review_status = fields[13]
                    stars = map_review_status(review_status)
                    phenotypeIDS = fields[7]
                    if stars >= int(evidence_level):
                        clinvar_id = fields[6]
                        clinvar_dct[variant] = {
                            "VariantName": variant_name,
                            "Gene": gene,
                            "ClinicalSignificance": clinical_significance,
                            "ClinSigSimple": clinsigsimple,
                            "rs": 'rs'+ rs_id,
                            "ReviewStatus": '(' + str(stars) + ') ' + review_status,
                            "ClinvarID": clinvar_id,
                            "PhenotypeIDS": phenotypeIDS
                        }
                        all_clinvar_id.append(clinvar_id)

        # For a given clinvar entry from clinvar_dct, include an aggregated results of Clinical Significance for each entry

        clinvar_ids = set(map(str, all_clinvar_id))
        clinical_significance_data = {vid: Counter() for vid in clinvar_ids}

        # Read submission summary file and get data
        with gzip.open(clinvar_submission, 'rt', encoding='utf-8') as f:
            lines = f.readlines()

        # Get the last header line which has column names
        header_line = [line for line in lines if line.startswith("#")][-1]
        header = header_line.lstrip("#").strip().split("\t")


        # Move through data (lines not starting with #)
        data_lines = [line for line in lines if not line.startswith("#")]

        # Build a CSV
        reader = csv.DictReader(data_lines, fieldnames=header, delimiter='\t')

        for row in reader:
            var_id = row['VariationID']
            if var_id in clinvar_ids:
                cs = row['ClinicalSignificance'].strip()
                if cs and row["ContributesToAggregateClassification"] == "yes": # Only get entries which contribute to the aggregate classification
                    clinical_significance_data[var_id][cs] += 1

        # Asign statistics to each entry
        for entry in clinvar_dct.values():
            var_id = str(entry.get("ClinvarID"))
            if var_id in clinical_significance_data:
                counter = clinical_significance_data[var_id]
                summary = "; ".join(f"{label} ({count})" for label, count in counter.items()) if counter else "No data"
                entry['ClinSigSummary'] = summary
            else:
                entry['ClinSigSummary'] = ""

        # Save dictionary to JSON file
        with output_clinvar_file.open("w", encoding="utf-8") as fh:
            json.dump(
                clinvar_dct,
                fh,
                indent=2,
                sort_keys=True,
                ensure_ascii=False
            )


    except Exception as e:
        print(f"Error when filtering variants: {e}")


def download_bundled_clinvar_snapshot(
        *,
        output_root: Path,
        clinvar_specification: dict[str, Any],
        overwrite: bool = False,
) -> dict[str, str]:
    """
    Download the ClinVar files pinned by the bundled resource specification.

    The source files are stored under::

        <output_root>/clinvar/source/

    Parameters
    ----------
    output_root
        Root directory of the installed SFtool resources.
    clinvar_specification
        ``clinvar`` section loaded from ``bundled_resources.json``.
    overwrite
        Replace previously downloaded files when true.

    Returns
    -------
    dict[str, str]
        ClinVar version, source URLs, and downloaded paths.

    Raises
    ------
    ResourceSpecificationError
        If the bundled ClinVar definition is incomplete or invalid.
    ResourceOperationError
        If a source file cannot be downloaded.
    """
    output_root = Path(output_root)
    source_directory = ensure_directory(
        output_root / "clinvar" / "source"
    )

    try:
        version = clinvar_specification["version"]
        archive_base_url = clinvar_specification[
            "archive_base_url"
        ]
        variant_summary_template = clinvar_specification[
            "variant_summary_filename"
        ]
        submission_summary_template = clinvar_specification[
            "submission_summary_filename"
        ]
    except KeyError as error:
        raise ResourceSpecificationError(
            "The bundled ClinVar specification is missing "
            f"the required field: {error.args[0]}"
        ) from error

    variant_summary_filename = render_resource_filename(
        variant_summary_template,
        version=version,
    )
    submission_summary_filename = render_resource_filename(
        submission_summary_template,
        version=version,
    )

    variant_summary_url = resolve_resource_url(
        archive_base_url,
        variant_summary_filename,
    )
    submission_summary_url = resolve_resource_url(
        archive_base_url,
        submission_summary_filename,
    )

    variant_summary_path = download_file(
        variant_summary_url,
        source_directory / variant_summary_filename,
        overwrite=overwrite,
        )
    submission_summary_path = download_file(
        submission_summary_url,
        source_directory / submission_summary_filename,
        overwrite=overwrite,
        )

    return {
        "version": version,
        "variant_summary": str(variant_summary_path),
        "submission_summary": str(submission_summary_path),
        "variant_summary_url": variant_summary_url,
        "submission_summary_url": submission_summary_url,
    }

def process_clinvar_data(assembly, release_date, clinvar_path):
    """
    Download and process CLINVAR database for a given assembly version
    
    Args:
        assembly (str): Assembly version. Either GRCh37 or GRCh38
        release_date (datetime.datetime): Clinvar release date
        clinvar_path: Path to clinvar directory
    
    Returns:
        str: Output file with processed data
    
    Raises:
        Exception: When an error occurs

    """
    # Columns of interest
    columns_of_interest_names = ["Type", "Name", "GeneSymbol",
                                 "ClinicalSignificance", "ClinSigSimple", "RS# (dbSNP)", "VariationID",
                                 "PhenotypeIDS", "PhenotypeList", "Assembly", 
                                 "Chromosome", "Start", "Stop", "ReviewStatus", 
                                 "SubmitterCategories", "PositionVCF", 
                                 "ReferenceAlleleVCF", "AlternateAlleleVCF"]
    
    # Output file
    output_file = f"{clinvar_path}clinvar_database_{assembly}_{release_date.strftime('%Y%m%d')}.txt"
    
    # Process CLINVAR file
    with gzip.open(f"{clinvar_path}variant_summary.txt.gz", "rt") as gz_file, open(output_file, "w") as output:
        csv_writer = csv.writer(output, delimiter="\t")
        header_line = gz_file.readline().strip()
        header_fields = header_line.split("\t")
        columns_of_interest_positions = [header_fields.index(col) for col in columns_of_interest_names]
        
        # Find the index of the "Assembly" column
        for idx, field in enumerate(header_fields):
            if field == "Assembly":
                assembly_column_index = idx
                break
    
        # Add columns of interest to the output file
        csv_writer.writerow(columns_of_interest_names)
    
        for line in gz_file:
            row = line.strip().split("\t")
            if row[assembly_column_index] == assembly:  # Filter according to genome version (GRCh37 o GRCh38)
                relevant_fields = [row[pos] for pos in columns_of_interest_positions]
                csv_writer.writerow(relevant_fields)
    
    return output_file


def get_clinvar(clinvar_path, assembly):
    """
    Download and process Clinvar database: include variant file and submission summary
    
    Args:
        clinvar_path: Path to CLINVAR directory database
    """
    try:        
        # CLINVAR URL: variant's file
        clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        
        # Open URL
        response = urllib.request.urlopen(clinvar_url)

        # Check whether response is OK (HTTP 200 code)
        if response.status != 200:
            print(f"Error downloading CLINVAR variants file. HTTP code: {response.status}")
            exit(1)
        
        # Open a local file for writing in binary mode
        out_rawfile = f"{clinvar_path}variant_summary.txt.gz"
        with open(out_rawfile, 'wb') as output_file:
            # Copy the response content to the local file
            shutil.copyfileobj(response, output_file)
        print(f"File downloaded to {out_rawfile}")
        
        # Get Clinvar release date
        last_modified = response.headers['Last-Modified']
        if last_modified is None:
            print("Clinvar release date cannot be obtained")
            exit(1)
        
        release_date = datetime.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
        
        # Process CLINVAR file for the assembly
        if assembly == "GRCh37":
            clinvar_variant_output_file = process_clinvar_data("GRCh37", release_date, clinvar_path)
            print(f"CLINVAR GRCh37 file is downloaded and processed. Version: {release_date.strftime('%Y%m%d')}")
        else:  # Assembly 38
            clinvar_variant_output_file = process_clinvar_data("GRCh38", release_date, clinvar_path)
            print(f"CLINVAR GRCh38 file is downloaded and processed. Version: {release_date.strftime('%Y%m%d')}")
        
        # Remove donwloaded file
        os.remove(f"{clinvar_path}variant_summary.txt.gz")

        # Download summary file with summaries for each clinvar entry
        clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"

        # Open URL
        response = urllib.request.urlopen(clinvar_url)

        # Check whether response is OK (HTTP 200 code)
        if response.status != 200:
            print(f"Error downloading CLINVAR submission summary. HTTP code: {response.status}")
            exit(1)

        # Open a local file for writing in binary mode
        clinvar_variant_summary_output_file = f"{clinvar_path}clinvar_submission_" + str(release_date.strftime('%Y%m%d')) +".txt.gz"
        with open(clinvar_variant_summary_output_file, 'wb') as output_file:
            # Copy the response content to the local file
            shutil.copyfileobj(response, output_file)
        print(f"File downloaded to {clinvar_variant_summary_output_file}")

        return [clinvar_variant_output_file, clinvar_variant_summary_output_file]
    
    except Exception as e:
        print(f"Error found: {str(e)}")




def clinvar_manager(clinvar_path, clinvar_ddbb_version, assembly):
    """

    Manage clinvar database: get version contains in the config file or download the latest version

    :param clinvar_path: directory path where Clinvar dabatase is stored
    :param clinvar_ddbb_version: clinvar version. Either a date in the form of YYYYMMDD or 'latest' for downloading the latest one
    :param assembly: reference genome version, either 37 or 38
    :return:
    """


    if clinvar_ddbb_version == "latest": # Download the latest version
        print("Downloading the latest version of Clinvar...")
        [clinvar_db, clinvar_summary_db] = get_clinvar(clinvar_path, assembly)
    else:  # Use the version contained in the config file
        print("Using existing Clinvar database (version " + clinvar_ddbb_version +")...")
        clinvar_file = os.path.join(clinvar_path, "clinvar_database_" + str(assembly) + "_" + clinvar_ddbb_version + ".txt")
        clinvar_summary_file = os.path.join(clinvar_path, "clinvar_submission" + "_" + clinvar_ddbb_version + ".txt.gz")
        if os.path.exists(clinvar_file) and os.path.exists(clinvar_summary_file):
            clinvar_db = clinvar_file
            clinvar_summary_db = clinvar_summary_file
        else:
            print(clinvar_file + "and/or" + clinvar_summary_file + "do not exist")
            exit(1)
    return [clinvar_db, clinvar_summary_db]
