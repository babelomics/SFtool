
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 19:07:52 2023

@author: Javier Perez Florido, Edurne Urrutia
"""
import csv
import time
from pathlib import Path
import shutil
from natsort import natsorted
import requests
import json
from typing import NamedTuple
from importlib.resources import as_file, files
from sftool.utils.resource_utils import (
    ensure_directory,
    write_json,
)

GENERATED_CATALOGS = ("PR", "RR")

CATALOG_VERSIONS = { "PR": "ACMG_SF_3.1", "RR": "ACMG_CS_2021", "RR_STR": "ACMG_CS_STR_2021", }

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

RETRYABLE_STATUS_CODES = {
    429,
    500,
    502,
    503,
    504,
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



def read_catalog_json(
        catalog_json: Path | str,
        category: str,
) -> tuple[dict, list[str]]:
    """
    Load an installed PR or RR catalog JSON resource.

    The JSON resource is generated during ``sftool resources setup``
    and contains the catalog metadata and its list of genes.

    Parameters
    ----------
    catalog_json
        Path to the installed catalog JSON resource.
    category
        Supported catalog category: PR or RR.

    Returns
    -------
    tuple[dict, list[str]]
        The complete catalog dictionary and the ordered list of gene
        symbols.
    """
    if category not in GENERATED_CATALOGS:
        raise CatalogGenerationError(
            f"Unsupported generated catalog: {category}"
        )

    catalog_path = Path(catalog_json)

    try:
        with catalog_path.open(
                "r",
                encoding="utf-8",
        ) as handle:
            catalog = json.load(handle)

    except OSError as error:
        raise CatalogGenerationError(
            f"Could not read catalog JSON: {catalog_path}"
        ) from error

    except json.JSONDecodeError as error:
        raise CatalogGenerationError(
            f"Invalid catalog JSON: {catalog_path}: {error}"
        ) from error

    if not isinstance(catalog, dict):
        raise CatalogGenerationError(
            "Catalog JSON must contain an object: "
            f"{catalog_path}"
        )

    genes_data = catalog.get("genes")

    if not isinstance(genes_data, list):
        raise CatalogGenerationError(
            "Catalog JSON field 'genes' must contain a list: "
            f"{catalog_path}"
        )

    genes: list[str] = []

    for index, gene_entry in enumerate(genes_data):
        if not isinstance(gene_entry, dict):
            raise CatalogGenerationError(
                "Catalog JSON gene entry must be an object at "
                f"index {index}: {catalog_path}"
            )

        gene_symbol = gene_entry.get("gene_symbol")

        if (
                not isinstance(gene_symbol, str)
                or not gene_symbol.strip()
        ):
            raise CatalogGenerationError(
                "Catalog JSON gene entry contains an invalid "
                f"'gene_symbol' at index {index}: {catalog_path}"
            )

        genes.append(gene_symbol.strip())

    return catalog, genes

def _get_ensembl_response(
        url: str,
        *,
        timeout: float,
        max_attempts: int = 5,
) -> requests.Response:
    """
    Send a GET request to Ensembl, retrying transient network and HTTP errors.

    Retries are performed for:

    - Connection errors
    - Timeouts
    - HTTP 429
    - HTTP 500
    - HTTP 502
    - HTTP 503
    - HTTP 504

    Non-transient HTTP errors, such as 400 or 404, are raised immediately.
    """
    if max_attempts < 1:
        raise ValueError("max_attempts must be greater than or equal to 1")

    last_error: requests.RequestException | None = None

    for attempt in range(1, max_attempts + 1):
        response: requests.Response | None = None

        try:
            response = requests.get(
                url,
                timeout=timeout,
            )

            if response.status_code not in RETRYABLE_STATUS_CODES:
                response.raise_for_status()
                return response

            response_body = response.text.strip()[:500]

            last_error = requests.HTTPError(
                (
                        f"Ensembl returned HTTP {response.status_code}"
                        + (
                            f": {response_body}"
                            if response_body
                            else ""
                        )
                ),
                response=response,
            )

        except (
                requests.Timeout,
                requests.ConnectionError,
        ) as error:
            last_error = error

        except requests.RequestException:
            # Do not retry non-transient errors such as HTTP 400 or 404.
            raise

        if attempt == max_attempts:
            break

        retry_after = None

        if response is not None:
            retry_after_header = response.headers.get("Retry-After")

            if retry_after_header:
                try:
                    retry_after = float(retry_after_header)
                except ValueError:
                    retry_after = None

        wait_seconds = (
            retry_after
            if retry_after is not None
            else 2 ** (attempt - 1)
        )

        time.sleep(wait_seconds)

    if last_error is None:
        raise requests.RequestException(
            f"Ensembl request failed for an unknown reason. URL: {url}"
        )

    raise last_error


def get_gene_location_ensembl(
        gene_symbol: str,
        assembly: str,
        *,
        timeout: float = 30.0,
        max_attempts: int = 5,
) -> dict[str, object]:
    """
    Retrieve genomic coordinates for a gene symbol from Ensembl.

    Transient network and server errors are retried automatically.
    """
    try:
        server = ENSEMBL_SERVERS[assembly]
    except KeyError as error:
        raise CatalogGenerationError(
            f"Unsupported genome assembly: {assembly!r}"
        ) from error

    url = (
        f"{server}/lookup/symbol/human/{gene_symbol}"
        "?content-type=application/json"
    )

    try:
        response = _get_ensembl_response(
            url,
            timeout=timeout,
            max_attempts=max_attempts,
        )

    except requests.Timeout as error:
        raise CatalogGenerationError(
            f"Ensembl request timed out while retrieving coordinates "
            f"for gene {gene_symbol!r} using {assembly}. "
            f"Attempts: {max_attempts}. "
            f"Timeout per attempt: {timeout} seconds. "
            f"URL: {url}. "
            f"Underlying error: {error}"
        ) from error

    except requests.ConnectionError as error:
        raise CatalogGenerationError(
            f"Could not connect to Ensembl while retrieving coordinates "
            f"for gene {gene_symbol!r} using {assembly}. "
            f"Attempts: {max_attempts}. "
            f"Check the network connection and Ensembl server availability. "
            f"URL: {url}. "
            f"Underlying error: {error}"
        ) from error

    except requests.HTTPError as error:
        error_response = error.response

        status_code = (
            error_response.status_code
            if error_response is not None
            else None
        )

        response_body = (
            error_response.text.strip()[:500]
            if error_response is not None
            else ""
        )

        if status_code == 400:
            reason = "Ensembl rejected the request"
        elif status_code == 404:
            reason = (
                f"Gene {gene_symbol!r} was not found by Ensembl "
                f"for {assembly}"
            )
        elif status_code == 429:
            reason = (
                "Ensembl rate limit exceeded after all retry attempts"
            )
        elif status_code in {500, 502, 503, 504}:
            reason = (
                "Ensembl returned a temporary server or gateway error "
                "after all retry attempts"
            )
        else:
            reason = "Ensembl returned an HTTP error"

        message = (
            f"{reason} while retrieving coordinates for gene "
            f"{gene_symbol!r} using {assembly}. "
            f"Attempts: {max_attempts}. "
            f"URL: {url}"
        )

        if status_code is not None:
            message += f". HTTP status: {status_code}"

        if response_body:
            message += f". Response: {response_body!r}"

        raise CatalogGenerationError(message) from error

    except requests.RequestException as error:
        raise CatalogGenerationError(
            f"Unexpected network error while retrieving coordinates "
            f"for gene {gene_symbol!r} using {assembly}. "
            f"Attempts: {max_attempts}. "
            f"URL: {url}. "
            f"Error type: {type(error).__name__}. "
            f"Underlying error: {error}"
        ) from error

    try:
        data = response.json()
    except ValueError as error:
        response_body = response.text.strip()[:500]

        message = (
            f"Invalid JSON returned by Ensembl for gene "
            f"{gene_symbol!r} using {assembly}. "
            f"HTTP status: {response.status_code}. "
            f"URL: {url}"
        )

        if response_body:
            message += f". Response: {response_body!r}"

        raise CatalogGenerationError(message) from error

    required_fields = {
        "seq_region_name",
        "start",
        "end",
    }

    if not isinstance(data, dict):
        raise CatalogGenerationError(
            f"Unexpected Ensembl response type for gene "
            f"{gene_symbol!r} using {assembly}: expected a JSON object, "
            f"received {type(data).__name__}. "
            f"URL: {url}"
        )

    missing_fields = required_fields.difference(data)

    if missing_fields:
        available_fields = ", ".join(
            sorted(map(str, data.keys()))
        )

        raise CatalogGenerationError(
            f"Incomplete Ensembl response for gene "
            f"{gene_symbol!r} using {assembly}. "
            f"Missing fields: {', '.join(sorted(missing_fields))}. "
            f"Available fields: {available_fields or 'none'}. "
            f"URL: {url}"
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
    chromosome = chromosome.strip()

    if chromosome.startswith("chr"):
        chromosome = chromosome[3:]

    if chromosome in {"M", "MT"}:
        return "chrM" if chr_prefix else "MT"

    return f"chr{chromosome}" if chr_prefix else chromosome

def write_bed_file(
        coordinates: list[GeneCoordinate],
        destination: Path,
        *,
        chr_prefix: bool,
) -> Path:
    temporary_path = destination.with_name(
        f"{destination.name}.tmp"
    )

    try:
        with temporary_path.open(
                "w",
                encoding="utf-8",
                newline="",
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

        temporary_path.replace(destination)

    except OSError as error:
        temporary_path.unlink(missing_ok=True)

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

def write_bed_variants(
        coordinates: list[GeneCoordinate],
        output_dir: Path,
        category: str,
) -> dict[str, Path]:
    """
    Write both BED chromosome conventions for one catalog.

    Returns paths for:
    - standard chromosome names: 1, 2, X, Y, MT
    - chr-prefixed names: chr1, chr2, chrX, chrY, chrM
    """
    bed_path = output_dir / f"{category}.bed"
    chr_bed_path = output_dir / f"{category}.chr.bed"

    write_bed_file(
        coordinates=coordinates,
        destination=bed_path,
        chr_prefix=False,
    )

    write_bed_file(
        coordinates=coordinates,
        destination=chr_bed_path,
        chr_prefix=True,
    )

    return {
        "bed": bed_path,
        "chr_bed": chr_bed_path,
    }
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

    bed_outputs = write_bed_variants(
        coordinates=coordinates,
        output_dir=output_dir,
        category=category,
    )

    json_path = output_dir / f"{category}.json"

    write_json(
        catalog_data,
        json_path,
    )

    return {
        **bed_outputs,
        "json": json_path,
    }

def prepare_catalog_resources(
        output_root: Path,
) -> dict[str, object]:
    output_root = ensure_directory(output_root)

    generated_catalogs: dict[str, object] = {
        "assemblies": {},
        "RR_STR": {
            "version": CATALOG_VERSIONS["RR_STR"],
        },
    }

    for assembly in ("GRCh37", "GRCh38"):
        assembly_output_dir = ensure_directory(
            output_root / "catalogs" / assembly
        )

        assembly_catalogs: dict[str, dict[str, object]] = {}

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

            assembly_catalogs[category] = {
                "version": CATALOG_VERSIONS[category],
                **{
                    name: str(path)
                    for name, path in outputs.items()
                },
            }

        generated_catalogs["assemblies"][assembly] = (
            assembly_catalogs
        )

    rr_str_path = copy_rr_str_catalog(
        output_root
    )

    generated_catalogs["RR_STR"]["csv"] = str(rr_str_path)

    return generated_catalogs