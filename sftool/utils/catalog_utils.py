
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 19:07:52 2023

@author: Javier Perez Florido, Edurne Urrutia
"""
import csv
from pathlib import Path
import shutil
from natsort import natsorted
import requests
from typing import NamedTuple
from importlib.resources import as_file
from sftool.utils.resource_utils import (
    ensure_directory,
    write_json,
)

GENERATED_CATALOGS = ("PR", "RR")

CATALOG_DESCRIPTIONS = {
    "PR": "Secondary findings of personal risk",
    "RR": "Secondary findings of reproductive risk",
}

GRCH37_GENE_ALIASES = {
    "MMUT": "MUT",
    "ELP1": "IKBKAP",
    "G6PC1": "G6PC",
    "GBA1": "GBA",
}

ENSEMBL_SERVERS = {
    "GRCh37": "https://grch37.rest.ensembl.org",
    "GRCh38": "https://rest.ensembl.org",
}

CATALOG_SOURCE_FILES = {
    "PR": (
        "PR",
        "PR_risk_genes_ACMG_SF_v3.1.csv",
    ),
    "RR": (
        "RR",
        "RR_risk_genes_ACMG_CS_v2021.csv",
    ),
    "RR_STR": (
        "RR",
        "RR_risk_genes_STR_ACMG_CS_v2021.csv",
    ),
}

class CatalogGenerationError(RuntimeError):
    """Raised when an SFtool catalog cannot be generated."""
class GeneCoordinate(NamedTuple):
    chromosome: str
    start: int
    end: int
    gene_symbol: str

def read_catalog_csv(
        source_csv: Path,
        category: str,
) -> tuple[dict, list[str]]:
    try:
        description = CATALOG_DESCRIPTIONS[category]
    except KeyError as error:
        raise CatalogGenerationError(
            f"Unsupported generated catalog: {category}"
        ) from error

    catalog = {
        "category": description,
        "genes": [],
    }

    genes: list[str] = []

    try:
        with source_csv.open(
                "r",
                encoding="latin1",
                newline="",
        ) as handle:
            reader = csv.DictReader(handle)

            if not reader.fieldnames or "Gene" not in reader.fieldnames:
                raise CatalogGenerationError(
                    f"Catalog CSV does not contain a 'Gene' column: "
                    f"{source_csv}"
                )

            for row in reader:
                gene_symbol = row["Gene"].strip()

                if not gene_symbol:
                    raise CatalogGenerationError(
                        f"Empty gene symbol found in {source_csv}"
                    )

                catalog["genes"].append(
                    {
                        "gene_symbol": gene_symbol,
                        "phenotype": row.get("Phenotype", ""),
                        "ACMG_version": row.get(
                            "ACMG SF List Version",
                            "",
                        ),
                        "OMIM_disorder": row.get(
                            "OMIM Disorder",
                            "",
                        ),
                        "inheritance": row.get(
                            "Inheritance",
                            "",
                        ),
                        "variants_to_report": row.get(
                            "Variants to Report",
                            "",
                        ),
                        "specific_variant_GRCh38": row.get(
                            "Specific variant GRCh38",
                            "",
                        ),
                        "specific_variant_GRCh37": row.get(
                            "Specific variant GRCh37",
                            "",
                        ),
                        "specific_consequence": row.get(
                            "Specific consequence",
                            "",
                        ),
                    }
                )

                genes.append(gene_symbol)

    except OSError as error:
        raise CatalogGenerationError(
            f"Could not read catalog CSV: {source_csv}"
        ) from error

    return catalog, genes


def get_gene_location_ensembl(
        gene_symbol: str,
        assembly: str,
        *,
        timeout: float = 30.0,
) -> dict[str, object]:
    try:
        server = ENSEMBL_SERVERS[assembly]
    except KeyError as error:
        raise CatalogGenerationError(
            f"Unsupported genome assembly: {assembly}"
        ) from error

    url = (
        f"{server}/lookup/symbol/human/{gene_symbol}"
        "?content-type=application/json"
    )

    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        data = response.json()
    except requests.RequestException as error:
        raise CatalogGenerationError(
            f"Could not retrieve coordinates for gene "
            f"{gene_symbol!r} using {assembly}"
        ) from error
    except ValueError as error:
        raise CatalogGenerationError(
            f"Invalid JSON returned by Ensembl for gene "
            f"{gene_symbol!r} using {assembly}"
        ) from error

    required_fields = {
        "seq_region_name",
        "start",
        "end",
    }

    if not isinstance(data, dict) or not required_fields.issubset(data):
        raise CatalogGenerationError(
            f"Incomplete Ensembl response for gene "
            f"{gene_symbol!r} using {assembly}"
        )

    return {
        "Gene_symbol": gene_symbol,
        "Chromosome": data["seq_region_name"],
        "Start": data["start"],
        "End": data["end"],
    }

def collect_gene_coordinates(
        genes: list[str],
        assembly: str,
) -> list[GeneCoordinate]:
    coordinates: list[GeneCoordinate] = []

    for gene_symbol in genes:
        query_symbol = gene_symbol

        if assembly == "GRCh37":
            query_symbol = GRCH37_GENE_ALIASES.get(
                gene_symbol,
                gene_symbol,
            )

        location = get_gene_location_ensembl(
            query_symbol,
            assembly,
        )

        coordinates.append(
            GeneCoordinate(
                chromosome=str(location["Chromosome"]),
                start=int(location["Start"]),
                end=int(location["End"]),
                gene_symbol=gene_symbol,
            )
        )

    return natsorted(
        coordinates,
        key=lambda coordinate: (
            coordinate.chromosome,
            coordinate.start,
            coordinate.end,
            coordinate.gene_symbol,
        ),
    )

def format_chromosome(
        chromosome: str,
        *,
        chr_prefix: bool,
) -> str:
    if not chr_prefix:
        return chromosome

    if chromosome == "MT":
        return "chrM"

    return f"chr{chromosome}"

def write_bed_file(
        coordinates: list[GeneCoordinate],
        destination: Path,
        *,
        chr_prefix: bool,
) -> Path:
    ensure_directory(destination.parent)

    try:
        with destination.open(
                "w",
                encoding="utf-8",
        ) as handle:
            for coordinate in coordinates:
                chromosome = format_chromosome(
                    coordinate.chromosome,
                    chr_prefix=chr_prefix,
                )

                handle.write(
                    f"{chromosome}\t"
                    f"{coordinate.start}\t"
                    f"{coordinate.end}\t"
                    f"{coordinate.gene_symbol}\n"
                )
    except OSError as error:
        raise CatalogGenerationError(
            f"Could not write BED file: {destination}"
        ) from error

    return destination


def get_bundled_catalog_resource(category: str):
    try:
        subdirectory, filename = CATALOG_SOURCE_FILES[category]
    except KeyError as error:
        raise CatalogGenerationError(
            f"Unsupported catalog: {category}"
        ) from error

    return (
        files("sftool.data.categories")
        .joinpath(subdirectory)
        .joinpath(filename)
    )

def copy_rr_str_catalog(
        output_root: Path,
        *,
        overwrite: bool = False,
) -> Path:
    source_resource = get_bundled_catalog_resource(
        "RR_STR"
    )

    destination_dir = ensure_directory(
        output_root / "catalogs" / "RR_STR"
    )

    destination = (
            destination_dir
            / "RR_risk_genes_STR_ACMG_CS_v2021.csv"
    )

    try:
        with as_file(source_resource) as source_path:
            shutil.copy2(
                source_path,
                destination,
            )
    except OSError as error:
        raise CatalogGenerationError(
            f"Could not copy RR_STR catalog to {destination}"
        ) from error

    return destination

def build_catalog_resources(
        category: str,
        assembly: str,
        source_csv: Path,
        output_dir: Path,
) -> dict[str, Path]:
    if category not in GENERATED_CATALOGS:
        raise CatalogGenerationError(
            f"Catalog {category!r} cannot be generated as BED/JSON"
        )

    if assembly not in ENSEMBL_SERVERS:
        raise CatalogGenerationError(
            f"Unsupported genome assembly: {assembly}"
        )

    output_dir = ensure_directory(output_dir)

    catalog_data, genes = read_catalog_csv(
        source_csv,
        category,
    )

    coordinates = collect_gene_coordinates(
        genes,
        assembly,
    )

    bed_path = output_dir / f"{category}.bed"
    chr_bed_path = output_dir / f"{category}.chr.bed"
    json_path = output_dir / f"{category}.json"

    write_bed_file(
        coordinates,
        bed_path,
        chr_prefix=False,
    )

    write_bed_file(
        coordinates,
        chr_bed_path,
        chr_prefix=True,
    )

    write_json(
        catalog_data,
        json_path,
    )

    return {
        "bed": bed_path,
        "chr_bed": chr_bed_path,
        "json": json_path,
    }

def prepare_catalog_resources(
        output_root: Path,
) -> dict[str, object]:
    output_root = ensure_directory(output_root)

    generated_catalogs: dict[str, object] = {
        "assemblies": {},
        "RR_STR": {},
    }

    for assembly in ("GRCh37", "GRCh38"):
        assembly_output_dir = ensure_directory(
            output_root / "catalogs" / assembly
        )

        generated_catalogs["assemblies"][assembly] = {}

        for category in GENERATED_CATALOGS:
            source_resource = get_bundled_catalog_resource(
                category
            )

            with as_file(source_resource) as source_csv:
                outputs = build_catalog_resources(
                    category=category,
                    assembly=assembly,
                    source_csv=source_csv,
                    output_dir=assembly_output_dir,
                )

            generated_catalogs["assemblies"][assembly][
                category
            ] = {
                name: str(path)
                for name, path in outputs.items()
            }

    rr_str_path = copy_rr_str_catalog(output_root)

    generated_catalogs["RR_STR"] = {
        "csv": str(rr_str_path),
    }

    return generated_catalogs